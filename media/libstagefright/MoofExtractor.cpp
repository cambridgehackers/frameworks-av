/*
 * Copyright (C) 2009 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//#define LOG_NDEBUG 0
#define LOG_TAG "MoofExtractor"
#include <utils/Log.h>

#include "include/MoofExtractor.h"
#include "include/SampleTable.h"
#include "include/ESDS.h"

#include <arpa/inet.h>

#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>
#include <semaphore.h>

#include <media/stagefright/foundation/ABitReader.h>
#include <media/stagefright/foundation/ADebug.h>
#include <media/stagefright/foundation/AMessage.h>
#include <media/stagefright/DataSource.h>
#include <media/stagefright/MediaBuffer.h>
#include <media/stagefright/MediaBufferGroup.h>
#include <media/stagefright/MediaDefs.h>
#include <media/stagefright/MediaSource.h>
#include <media/stagefright/MetaData.h>
#include <media/stagefright/TimeSource.h>
#include <media/stagefright/Utils.h>
#include <utils/String8.h>

namespace android {

    static void hexdump(const char *prefix, const void *_data, size_t size, bool ascii=false) {
    const uint8_t *data = (const uint8_t *)_data;
    size_t offset = 0;
    char line[512];
    int pos = 0;
    while (offset < size) {
        pos += snprintf(line+pos, sizeof(line)-pos, "%s", prefix);
        pos += snprintf(line+pos, sizeof(line)-pos, "0x%04x  ", offset);

        size_t n = size - offset;
        if (n > 32) {
            n = 32;
        }

        for (size_t i = 0; i < 32; ++i) {
            if ((i % 4) == 0) {
                pos += snprintf(line+pos, sizeof(line)-pos, " ");
            }

            if (offset + i < size) {
                pos += snprintf(line+pos, sizeof(line)-pos, "%02x", data[offset + i]);
            } else {
                pos += snprintf(line+pos, sizeof(line)-pos, "  ");
            }
        }

        pos += snprintf(line+pos, sizeof(line)-pos, " ");

        if (ascii) {
            for (size_t i = 0; i < n; ++i) {
                if (isprint(data[offset + i])) {
                    pos += snprintf(line+pos, sizeof(line)-pos, "%c", data[offset + i]);
                } else {
                    pos += snprintf(line+pos, sizeof(line)-pos, ".");
                }
            }
        }

        ALOGV("%s", line);
        pos = 0;

        offset += 32;
    }
    ALOGV("%s hexdump offset %d(0x%x) size %d", prefix, offset, offset, size);
}

static void MakeFourCCString(uint32_t x, char *s) {
    s[0] = x >> 24;
    s[1] = (x >> 16) & 0xff;
    s[2] = (x >> 8) & 0xff;
    s[3] = x & 0xff;
    s[4] = '\0';
}

size_t parseNALSize(uint32_t nalLengthSize, const uint8_t *data) {
    switch (nalLengthSize) {
        case 1:
            return *data;
        case 2:
            return U16_AT(data);
        case 3:
            return ((size_t)data[0] << 16) | U16_AT(&data[1]);
        case 4:
            return U32_AT(data);
    }

    // This cannot happen, mNALLengthSize springs to life by adding 1 to
    // a 2-bit integer.
    CHECK(!"Should not be here.");

    return 0;
}

// This custom data source wraps an existing one and satisfies requests
// falling entirely within a cached range from the cache while forwarding
// all remaining requests to the wrapped datasource.
// This is used to cache the full sampletable metadata for a single track,
// possibly wrapping multiple times to cover all tracks, i.e.
// Each MoofDataSource caches the sampletable metadata for a single track.

struct MoofDataSource : public DataSource {
    MoofDataSource(const sp<DataSource> &source);

    virtual status_t initCheck() const;
    virtual ssize_t readAt(off64_t offset, void *data, size_t size);
    virtual status_t getSize(off64_t *size);
    virtual uint32_t flags();

    status_t setCachedRange(off64_t offset, size_t size);
    void clearCache();

protected:
    virtual ~MoofDataSource();

private:
    Mutex mLock;

    sp<DataSource> mSource;
    off64_t mCachedOffset;
    size_t mCachedSize;
    uint8_t *mCache;

    MoofDataSource(const MoofDataSource &);
    MoofDataSource &operator=(const MoofDataSource &);
};

MoofDataSource::MoofDataSource(const sp<DataSource> &source)
    : mSource(source),
      mCachedOffset(0),
      mCachedSize(0),
      mCache(NULL) {
}

MoofDataSource::~MoofDataSource() {
    clearCache();
}

void MoofDataSource::clearCache() {
    if (mCache) {
        free(mCache);
        mCache = NULL;
    }

    mCachedOffset = 0;
    mCachedSize = 0;
}

status_t MoofDataSource::initCheck() const {
    return mSource->initCheck();
}

ssize_t MoofDataSource::readAt(off64_t offset, void *data, size_t size) {
    Mutex::Autolock autoLock(mLock);

    bool cacheHit = (offset >= mCachedOffset
                     && (offset + size) <= (mCachedOffset + mCachedSize));
    ALOGV("MoofDataSource::readAt %lld size %d cached[%lld:%lld] %s",
          offset, size, mCachedOffset, mCachedOffset+mCachedSize,
          cacheHit ? "" : "miss");

    if (cacheHit) {
        memcpy(data, &mCache[offset - mCachedOffset], size);
        return size;
    }

    return mSource->readAt(offset, data, size);
}

status_t MoofDataSource::getSize(off64_t *size) {
    ALOGV("MoofDataSource::getSize");
    return mSource->getSize(size);
}

uint32_t MoofDataSource::flags() {
    return mSource->flags();
}

status_t MoofDataSource::setCachedRange(off64_t offset, size_t size) {
    Mutex::Autolock autoLock(mLock);

    ALOGV("MoofDataSource::setCachedRange %lld size %d", offset, size);
    clearCache();

    mCache = (uint8_t *)malloc(size);

    if (mCache == NULL) {
        return -ENOMEM;
    }

    mCachedOffset = offset;
    mCachedSize = size;

    ssize_t err = mSource->readAt(mCachedOffset, mCache, mCachedSize);

    if (err < (ssize_t)size) {
        clearCache();

        return ERROR_IO;
    }

    return OK;
}

struct MoofSample {
    MoofSample() : sampleIndex(0), sampleSize(0), sampleOffset(0), sampleDuration(0), sampleTime(0),
                   buffer(0) {}
    void setBuffer(MediaBuffer *newBuffer) {
        if (buffer)
            buffer->release();
        buffer = newBuffer;
    }
    bool hasBuffer() const { return buffer != NULL; }
    MediaBuffer *takeBuffer() {
        MediaBuffer *result = buffer;
        buffer = NULL;
        return result;
    }

    uint32_t sampleIndex;
    uint32_t sampleSize;
    uint64_t sampleOffset;
    uint64_t sampleDuration; // ticks
    uint64_t sampleTime;     // microseconds
    uint32_t sampleFlags;
    uint32_t cts;
    double   mfts;
    sp<MetaData> meta;
private:
    MediaBuffer *buffer;
};

class MoofSampleTable : public RefBase, MediaBufferObserver {
public:
    MoofSampleTable(uint64_t offset, sp<MoofDataSource> &dataSource);
    ~MoofSampleTable();
    int getSample(size_t trackIndex, size_t sampleIndex, MoofSample *sample);
    int setTimeScale(size_t trackIndex, uint32_t timeScale);
    void setTrackMetaData(size_t trackIndex, sp<MetaData> &meta) { mMetaData[trackIndex] = meta; }
    void onMetaDataReady();

    // MediaBufferObserver
    virtual void signalBufferReturned(MediaBuffer *buffer);
private:
    enum {
        kMfts = 'mfts'
    };
    status_t acquireBuffer(MediaBuffer **buffer, int track, uint64_t size);
    status_t readSample(const sp<MetaData> &format, MoofSample *sample, MediaBuffer *buffer);
    int readSamples();

    static void *threadWrapper(void *me);
    void run();

private:
    Mutex mSampleLock[2];
    pthread_t mThread;
    sem_t mProducerSemaphore;
    sem_t mConsumerSemaphore;
    bool  mDone;
    double mSampleTableTimeStamp;
    double mFirstMfts;

    List<MoofSample> mSamples[2];
    size_t mTrackIndexes[2];
    uint64_t mOffset;
    uint32_t mUsPerTick[2];
    uint64_t mLastDecodeTime[2];
    sp<MoofDataSource> mDataSource;
    MediaBufferGroup *mGroup[2];
    uint32_t mSampleIndex[2];
    sp<MetaData> mMetaData[2];
};

MoofSampleTable::MoofSampleTable(uint64_t offset, sp<MoofDataSource> &dataSource)
    : mOffset(offset)
    , mDataSource(dataSource)
{
    sem_init(&mConsumerSemaphore, 0, 1);
    sem_init(&mProducerSemaphore, 0, 0);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    ALOGV("starting thread");
    mDone = false;
    pthread_create(&mThread, &attr, threadWrapper, this);
    pthread_attr_destroy(&attr);

    memset(mTrackIndexes, 0, sizeof(mTrackIndexes));
    memset(mLastDecodeTime, 0, sizeof(mLastDecodeTime));
    memset(mUsPerTick, 0, sizeof(mUsPerTick));
    memset(mSampleIndex, 0, sizeof(mSampleIndex));
    for (int i = 0; i < 2; i++) {
        mGroup[i] = new MediaBufferGroup;
        // prime the pump
        for (int j = 0; j < 100; j++) {
            MediaBuffer *buffer = new MediaBuffer(1000);
            mGroup[i]->add_buffer(buffer);
        }
    }
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    mSampleTableTimeStamp = tv.tv_sec + tv.tv_usec / 1.0e6;
}

MoofSampleTable::~MoofSampleTable()
{
    mDone = true;
    void *ret;
    ALOGV("~MoofSampleTable: pthread_joining ...");
    pthread_join(mThread, &ret);
    ALOGV("~MoofSampleTable: thread done");

    sem_destroy(&mConsumerSemaphore);
    sem_destroy(&mProducerSemaphore);
}

void MoofSampleTable::onMetaDataReady()
{
    sem_post(&mProducerSemaphore);
}

void MoofSampleTable::signalBufferReturned(MediaBuffer *buffer)
{
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    double bufferMfts = 0;
    buffer->meta_data()->findDouble(kMfts, &bufferMfts);
    ALOGD("buffer mfts %.3f returned after %.3f seconds ",
          bufferMfts, (tv.tv_sec + tv.tv_usec / 1.0e6) - bufferMfts);

    //FIXME: should reuse the buffer
    buffer->setObserver(NULL);
    buffer->release();
}

void *MoofSampleTable::threadWrapper(void *self)
{
    MoofSampleTable *sampleTable = static_cast<MoofSampleTable *>(self);
    sampleTable->run();
    return NULL;
}

void MoofSampleTable::run()
{
    char skipbuf[1024];
    mDataSource->readAt(0, skipbuf, mOffset);
    while (1) {
        ALOGV("MoofSampleTable::run() waiting for producer semaphore");
        sem_wait(&mProducerSemaphore);
        ALOGV("MoofSampleTable::run() producer semaphore proceeding");
        readSamples();
        sem_post(&mConsumerSemaphore);
    }
}

int MoofSampleTable::acquireBuffer(MediaBuffer **pbuffer, int track, uint64_t size)
{
    MediaBuffer *buffer = new MediaBuffer(size);
    buffer->add_ref();
    buffer->setObserver(this);
    ALOGV("buffer %p size() %d size %lld sample lock %d", buffer, buffer->size(), size, track);
    CHECK(buffer != NULL);
    CHECK(buffer->size() >= size);

    *pbuffer = buffer;
    return OK;
}

status_t MoofSampleTable::readSample(const sp<MetaData> &format,
                                     MoofSample *sample, MediaBuffer *buffer)
{
    ALOGV("readSample format=%p sample=%p buffer=%p\n", &format, sample, buffer);
    sample->setBuffer(buffer);
    uint64_t offset = sample->sampleOffset;
    uint64_t size = sample->sampleSize;
    bool isSyncSample = ((sample->sampleFlags & 0x00010000) == 0);

    const char *mime;
    bool success = format->findCString(kKeyMIMEType, &mime);
    CHECK(success);

    bool isAVC = !strcasecmp(mime, MEDIA_MIMETYPE_VIDEO_AVC);
    int32_t drm = 0;
    bool usesDRM = (format->findInt32(kKeyIsDRM, &drm) && drm != 0);
    int prefixLen = 0;
    if (isAVC && !usesDRM)
        prefixLen = 4;

    ssize_t num_bytes_read =
        mDataSource->readAt(offset, (uint8_t *)buffer->data(), size);
    ALOGV("readSample offset %lld size %lld read %ld", offset, size, num_bytes_read);

    if (num_bytes_read < (ssize_t)size) {
        ALOGV("ERROR_IO num_bytes_read=%ld", num_bytes_read);
        buffer->release();
        buffer = NULL;
        return ERROR_IO;
    }


    if (!isAVC) {

        buffer->set_range(0, size);
        buffer->meta_data()->clear();
        buffer->meta_data()->setInt64(kKeyTime, sample->sampleTime);
        buffer->meta_data()->setDouble(kMfts, sample->mfts);

        uint64_t targetSampleTimeUs = sample->sampleTime;
        buffer->meta_data()->setInt64(
            kKeyTargetTime, targetSampleTimeUs);

        if (isSyncSample) {
            buffer->meta_data()->setInt32(kKeyIsSyncFrame, 1);
        }

        return OK;
    } else {
        // Whole NAL units are returned but each fragment is prefixed by
        // the start code (0x00 00 00 01).

        uint32_t nalLengthSize = 0;
        {
            uint32_t type;
            const void *data;
            size_t size;
            CHECK(format->findData(kKeyAVCC, &type, &data, &size));

            const uint8_t *ptr = (const uint8_t *)data;

            CHECK(size >= 7);
            CHECK_EQ((unsigned)ptr[0], 1u);  // configurationVersion == 1

            // The number of bytes used to encode the length of a NAL unit.
            nalLengthSize = 1 + (ptr[4] & 3);
        }

        if (usesDRM) {
            CHECK(buffer != NULL);
            buffer->set_range(0, size);

        } else {
            const uint8_t *srcBuffer = (const uint8_t *)buffer->data();
            size_t srcOffset = 0;
            size_t numNalUnits = 0;

            bool isMalFormed = (srcOffset + nalLengthSize > size);
            size_t nalLength = 0;
            if (!isMalFormed) {
                nalLength = parseNALSize(nalLengthSize, &srcBuffer[srcOffset]);
                srcOffset += nalLengthSize;
                isMalFormed = srcOffset + nalLength > size;
            }
            ALOGV("%s:%d isMalFormed %d nalLengthSize %d nalLength %d srcOffset %d size %lld", __FILE__, __LINE__,
                  isMalFormed, nalLengthSize, nalLength, srcOffset, size);

            if (isMalFormed) {
                ALOGE("Video is malformed");
                //buffer->release();
                //buffer = NULL;
                //return ERROR_MALFORMED;
            }

            //CHECK(nalLength != 0);
            CHECK(nalLengthSize == prefixLen);
            // overwrite NAL length with start code 00 00 00 01
            uint8_t *dstData = (uint8_t *)buffer->data(); //+prefixLen-nalLengthSize;
            size_t dstOffset = 0;
            CHECK(dstOffset + 4 <= buffer->size());

            dstData[dstOffset++] = 0;
            dstData[dstOffset++] = 0;
            dstData[dstOffset++] = 0;
            dstData[dstOffset++] = 1;

            numNalUnits++;

            //ALOGV("copied %d nal units %lld bytes", numNalUnits, size-nalLengthSize+prefixLen);
            //CHECK_EQ(srcOffset, size);
            CHECK(buffer != NULL);
            CHECK(size <= buffer->size());
            buffer->set_range(0, size);
            if (1)
            ALOGV("set_range prefixLen %d nalLengthSize %d offset %d size %lld buffer size %d",
                  prefixLen, nalLengthSize, 0, size, buffer->size());
            //hexdump(buffer->data(), size+prefixLen);
        }

        buffer->meta_data()->clear();
        buffer->meta_data()->setInt64(kKeyTime, sample->sampleTime);
        buffer->meta_data()->setDouble(kMfts, sample->mfts);

        // FIXME
        int64_t targetSampleTimeUs = sample->sampleTime;
        if (targetSampleTimeUs >= 0) {
            buffer->meta_data()->setInt64(kKeyTargetTime, targetSampleTimeUs);
        }

        if (isSyncSample) {
            buffer->meta_data()->setInt32(kKeyIsSyncFrame, 1);
        }

        return OK;
    }
}

int MoofSampleTable::readSamples()
{
    
    bool done = false;
    uint32_t trackIndex = 0;

    uint32_t sampleCount = 0;
    uint32_t defaultSampleDuration = 0;
    uint32_t defaultSampleSize = 0;
    uint32_t defaultSampleFlags = 0;
    uint64_t decodeTimeUs;

    uint64_t sampleDataOffset = 0;
    double mfts_source_time = 0;

    ALOGV("readSamples mOffset=%lld", mOffset);
    while (!done) {
        uint32_t hdr[2];
        size_t bytesRead = 0;
        if ((bytesRead = mDataSource->readAt(mOffset, hdr, 8)) < 8) {
            ALOGE("short read @ %lld / %d", mOffset, bytesRead);
            return ERROR_IO;
        }
        uint64_t chunk_size = ntohl(hdr[0]);
        uint32_t chunk_type = ntohl(hdr[1]);
        off64_t data_offset = mOffset + 8;
        uint64_t chunk_data_size = chunk_size - 8;
        if (chunk_size == 1) {
            if (mDataSource->readAt(mOffset + 8, &chunk_size, 8) < 8) {
                ALOGE("%s:%d: short read @ %lld", __FILE__, __LINE__, mOffset);
                return ERROR_IO;
            }
            chunk_size = ntoh64(chunk_size);
            data_offset += 8;
            chunk_data_size -= 8;

            if (chunk_size < 16) {
                // The smallest valid chunk is 16 bytes long in this case.
                ALOGV("malformed offset %lld chunk_size %lld", mOffset, chunk_size);
                return ERROR_MALFORMED;
            }
        } else if (chunk_size < 8) {
            // The smallest valid chunk is 8 bytes long.
            ALOGV("malformed offset %lld chunk_size %lld", mOffset, chunk_size);
            return ERROR_MALFORMED;
        }
        char chunk[5];
        MakeFourCCString(chunk_type, chunk);
        ALOGV("%s @ %lld size %lld readSamples",
              chunk, mOffset, chunk_size);

        switch(chunk_type) {
        case FOURCC('m','o','o','f'):
            mDataSource->setCachedRange(data_offset, chunk_data_size);
            mOffset = data_offset;
            break;
        case FOURCC('t','r','a','f'):
            // parse the contents of this chunk
            mOffset = data_offset;
            break;
        case FOURCC('m','d','a','t'): {
            mOffset += chunk_size;
            mDataSource->setCachedRange(data_offset, chunk_data_size);
            // now read the samples
            int numSamplesRead = 0;
            for (int trackIndex = 0; trackIndex < 2; trackIndex++) {
                //ALOGV("track %d samples %d", trackIndex, mSamples[trackIndex].size());

                ALOGV("%s:%d: locking sample lock %d", __FILE__, __LINE__, trackIndex);
		mSampleLock[trackIndex].lock();

                for (List<MoofSample>::iterator it = mSamples[trackIndex].begin();
                     it != mSamples[trackIndex].end(); ++it) {
                    MoofSample &sample = *it;
                    uint64_t sampleOffset = sample.sampleOffset;
                    if (sampleOffset >= data_offset
                        && sampleOffset < data_offset + chunk_data_size) {
                        MediaBuffer *buffer;

                        ALOGV("sample lock %d acquiring buffer", trackIndex);
                        acquireBuffer(&buffer, trackIndex, chunk_data_size);

                        ALOGV("sample lock %d read sample %d", trackIndex, sample.sampleIndex);
                        readSample(mMetaData[trackIndex], &sample, buffer);

                        numSamplesRead++;
                    }
                }

                ALOGV("%s:%d: unlocking sample lock %d", __FILE__, __LINE__, trackIndex);
                mSampleLock[trackIndex].unlock();

            }
            mDataSource->clearCache();
            if (!numSamplesRead)
                ALOGE("read %d samples at offset %lld", numSamplesRead, data_offset);
            return OK;
        } break;
        case FOURCC('m','f','t','s'): {
            uint32_t box[2];
            if (mDataSource->readAt(data_offset, box, chunk_size-8) < (chunk_size-8)) {
                ALOGE("%s:%d readAt short read", __FILE__, __LINE__);
                return ERROR_IO;
            }
            struct timeval tv;
            struct timezone tz;
            gettimeofday(&tv, &tz);
            box[0] = ntohl(box[0]);
            box[1] = ntohl(box[1]);
            mfts_source_time = box[0] + box[1] / 1.0e6;
            double local_time = tv.tv_sec + tv.tv_usec / 1.0e6;
            ALOGD("mfts delta %.3g source %x.%x local %lx.%lx",
                  local_time - mfts_source_time, box[0], box[1], tv.tv_sec, tv.tv_usec);
            mOffset += chunk_size;
        } break;
        case FOURCC('t','f','h','d'):
        case FOURCC('t','r','u','n'):
        case FOURCC('t','f','d','t'): {
            uint32_t box[256];
            if (mDataSource->readAt(data_offset, box, chunk_size-8) < (chunk_size-8)) {
                ALOGE("%s:%d readAt short read", __FILE__, __LINE__);
                return ERROR_IO;
            }
            uint32_t flags = ntohl(box[0]);
            switch (chunk_type) {
            case FOURCC('t','f','h','d'): {
                trackIndex = ntohl(box[1])-1;
                if (trackIndex >= 2)
                    ALOGE("trackIndex %d out of range", trackIndex);
                CHECK(trackIndex < 2);
                decodeTimeUs = mLastDecodeTime[trackIndex];
                enum tfhd_flags {
                    tfhd_default                         = 0x000000,
                    tfhd_base_data_offset_present        = 0x000001,
                    tfhd_default_sample_duration_present = 0x000008,
                    tfhd_default_sample_size_present     = 0x000010,
                    tfhd_default_sample_flags_present    = 0x000020,
                };

                int infoOffset = 2;
                if (flags & tfhd_base_data_offset_present) {
                    memcpy((char*)&sampleDataOffset, box+infoOffset, sizeof(sampleDataOffset));
                    infoOffset += 2;
                    sampleDataOffset = ntoh64(sampleDataOffset);
                }
                if (flags & tfhd_default_sample_duration_present) {
                    defaultSampleDuration = ntohl(box[infoOffset]);
                    infoOffset++;
                }
                if (flags & tfhd_default_sample_size_present) {
                    defaultSampleSize = ntohl(box[infoOffset]);
                    infoOffset++;
                }
                if (flags & tfhd_default_sample_size_present) {
                    defaultSampleFlags = ntohl(box[infoOffset]);
                    infoOffset++;
                }
                if (1)
                ALOGV("tfhd %lld flags %d trackIndex=%x dataOffset=%lld",
                      chunk_size, flags, trackIndex, sampleDataOffset);
            } break;
            case FOURCC('t','r','u','n'): {
                sampleCount = ntohl(box[1]);
                if (0)
                ALOGV("trun %lld %d sampleCount %d durationTicks %d durationSize %d",
                      chunk_size, flags, sampleCount, defaultSampleDuration, defaultSampleSize);
                enum trun_flags {
                    trun_default                    = 0x000000,
                    trun_data_offset_present        = 0x000001,
                    trun_first_sample_flags_present = 0x000004,
                    trun_sample_duration_present    = 0x000100,
                    trun_sample_size_present        = 0x000200,
                    trun_sample_flags_present       = 0x000400,
                    trun_sample_composition_time_present = 0x000800,
                };
                int wordNumber = 2; // first word of trun sample info
                MoofSample sample;
                sample.sampleTime = decodeTimeUs;
                sample.mfts = mfts_source_time;
                sample.sampleDuration = defaultSampleDuration;
                sample.sampleSize = defaultSampleSize;
                sample.sampleOffset = sampleDataOffset;
                if (flags & trun_data_offset_present) {
                    ALOGV("trun data offset present %d", ntohl(box[wordNumber]));
                    sample.sampleOffset = ntohl(box[wordNumber]);
                    wordNumber++;
                }
                if (flags & trun_first_sample_flags_present) {
                    sample.sampleFlags = ntohl(box[wordNumber]);
                    //FIXME
                    defaultSampleFlags = 0x00010000; // default is nonsync if first sample flags provided 
                    wordNumber++;
                }
                if (trackIndex >= 2)
                    ALOGE("trackIndex %d out of range", trackIndex);
                CHECK(trackIndex < 2);
                for (size_t i = 0; i < sampleCount; i++) {
                    sample.sampleIndex = ++mSampleIndex[trackIndex];
                    if (i != 0)
                        sample.sampleFlags = defaultSampleFlags;
                    if (flags & trun_sample_duration_present) {
                        sample.sampleDuration = ntohl(box[wordNumber]);
                        wordNumber++;
                    }
                    if (flags & trun_sample_size_present) {
                        sample.sampleSize = ntohl(box[wordNumber]);
                        wordNumber++;
                    }
                    if (flags & trun_sample_flags_present) {
                        sample.sampleFlags = ntohl(box[wordNumber]);
                        wordNumber++;
                    }
                    ALOGV("track %d sample %d ticks %llx dur %lld size %d flags %x offset %lld time %lld timescale %d",
                          trackIndex,
                          sample.sampleIndex,
                          sample.sampleDuration,
                          (uint64_t)(sample.sampleDuration * mUsPerTick[trackIndex]),
                          sample.sampleSize, sample.sampleFlags,
                          sample.sampleOffset, sample.sampleTime,
                          mUsPerTick[trackIndex]);

                    ALOGV("%s:%d: locking sample lock %d", __FILE__, __LINE__, trackIndex);
		    mSampleLock[trackIndex].lock();

                    mSamples[trackIndex].push_back(sample);

                    ALOGV("%s:%d: unlocking sample lock %d", __FILE__, __LINE__, trackIndex);
		    mSampleLock[trackIndex].unlock();

                    sample.sampleOffset += sample.sampleSize;
                    if (mUsPerTick[trackIndex] == 0)
                        ALOGE("zero mUsPerTick for track %d", trackIndex);
                    //CHECK(mUsPerTick[trackIndex] != 0);
                    sample.sampleTime += sample.sampleDuration * mUsPerTick[trackIndex];
                    mLastDecodeTime[trackIndex] = sample.sampleTime;
                }
            } break;
            case FOURCC('t','f','d','t'): {
                decodeTimeUs = ntohl(box[1]) * mUsPerTick[trackIndex];
                if (1)
                ALOGD("tfdt ticks decodetime %lld (us) delta %lld (us) ", 
                      decodeTimeUs,
                      decodeTimeUs - mLastDecodeTime[trackIndex]);
                mLastDecodeTime[trackIndex] = decodeTimeUs;
            } break;
            }
            mOffset += chunk_size;
        } break;
        default:
            mOffset += chunk_size;
        }
    }
    return OK;
}

int MoofSampleTable::getSample(size_t trackIndex, size_t sampleIndex, MoofSample *sample)
{
    ALOGV("%s:%d: locking sample lock %d", __FILE__, __LINE__, trackIndex);
    Mutex::Autolock autoLock(mSampleLock[trackIndex]);

    if (trackIndex >= 2)
        ALOGE("getSample trackIndex %d", trackIndex);
    CHECK(trackIndex < 2);

    int numTries = 0;
    do {
        for (List<MoofSample>::iterator it = mSamples[trackIndex].begin();
             it != mSamples[trackIndex].end(); ++it) {
            if (it->sampleIndex == sampleIndex) {
                *sample = *it;
                ALOGV("got track %d sample %d buffer %d", trackIndex, sampleIndex, sample->hasBuffer());
                if (!sample->hasBuffer())
                    break;
                mSamples[trackIndex].erase(it);
                ALOGV("%s:%d: unlocking sample lock %d", __FILE__, __LINE__, trackIndex);
                return OK;
            }
            if (it->sampleIndex < sampleIndex) {
                ALOGD("discarding track %d sample %d", trackIndex, it->sampleIndex);
                MediaBuffer *mb = it->takeBuffer();
                if (mb)
                    mb->release();
                mSamples[trackIndex].erase(it);
            }
        }

        ALOGV("%s:%d: unlocking sample lock %d", __FILE__, __LINE__, trackIndex);
        mSampleLock[trackIndex].unlock();

        // tell it to read more samples
        sem_post(&mProducerSemaphore);

        // and wait for the response
        ALOGV("getSample waiting for consumer semaphore");
        sem_wait(&mConsumerSemaphore);
        ALOGV("getSample consumer semaphore proceeding");

        ALOGV("%s:%d: locking sample lock %d", __FILE__, __LINE__, trackIndex);
        mSampleLock[trackIndex].lock();

    } while (numTries++ < 50); // totally arbitrary 50
    ALOGE("track %d sample %d out of range ", trackIndex, sampleIndex);
    ALOGV("%s:%d: unlocking sample lock %d", __FILE__, __LINE__, trackIndex);
    return -1;    
}

int MoofSampleTable::setTimeScale(size_t trackIndex, uint32_t timeScale)
{

    CHECK(trackIndex < 2);
    mUsPerTick[trackIndex] = 1.0e6 / timeScale;
    ALOGV("setTimeScale track %d timeScale %d us/tick %d",
          trackIndex, timeScale, mUsPerTick[trackIndex]);
    return OK;
}

class MoofSource : public MediaSource {
public:
    // Caller retains ownership of both "dataSource" and "sampleTable".
    MoofSource(size_t trackIndex,
               const sp<MetaData> &format,
               const sp<DataSource> &dataSource,
               uint32_t timescale,
               sp<MoofSampleTable> &moofSampleTable);

    virtual status_t start(MetaData *params = NULL);
    virtual status_t stop();

    virtual sp<MetaData> getFormat();

    virtual status_t read(
            MediaBuffer **buffer, const ReadOptions *options = NULL);

protected:
    virtual ~MoofSource();

private:
    Mutex mLock;

    sp<MetaData> mFormat;
    sp<DataSource> mDataSource;
    int32_t mTimescale;
    uint32_t mCurrentSampleIndex;
    SystemTimeSource mTimeSource;
    uint64_t mLastTimeReadUs;
    size_t mTrackIndex;
    sp<MoofSampleTable> mSampleTable;

    bool mIsAVC;
    size_t mNALLengthSize;

    bool mStarted;

    MediaBufferGroup *mGroup;

    MediaBuffer *mBuffer;

    bool mWantsNALFragments;

    uint8_t *mSrcBuffer;
    size_t   mSrcBufferSize;

    MoofSource(const MoofSource &);
    MoofSource &operator=(const MoofSource &);
};

////////////////////////////////////////////////////////////////////////////////

static const char *FourCC2MIME(uint32_t fourcc) {
    switch (fourcc) {
        case FOURCC('m', 'p', '4', 'a'):
            return MEDIA_MIMETYPE_AUDIO_AAC;

        case FOURCC('s', 'a', 'm', 'r'):
            return MEDIA_MIMETYPE_AUDIO_AMR_NB;

        case FOURCC('s', 'a', 'w', 'b'):
            return MEDIA_MIMETYPE_AUDIO_AMR_WB;

        case FOURCC('m', 'p', '4', 'v'):
            return MEDIA_MIMETYPE_VIDEO_MPEG4;

        case FOURCC('s', '2', '6', '3'):
        case FOURCC('h', '2', '6', '3'):
        case FOURCC('H', '2', '6', '3'):
            return MEDIA_MIMETYPE_VIDEO_H263;

        case FOURCC('a', 'v', 'c', '1'):
            return MEDIA_MIMETYPE_VIDEO_AVC;

        default:
            CHECK(!"should not be here.");
            return NULL;
    }
}

MoofExtractor::MoofExtractor(const sp<DataSource> &source)
    : mDataSource(new MoofDataSource(source)),
      mInitCheck(NO_INIT),
      mHasVideo(false),
      mNumTracks(0),
      mFirstTrack(NULL),
      mLastTrack(NULL),
      mFileMetaData(new MetaData),
      mFirstSINF(NULL),
      mIsDrm(false) {
    ALOGV("MoofExtractor");
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    mStartupTimeStamp = tv.tv_sec + tv.tv_usec / 1.0e6;    
}

MoofExtractor::~MoofExtractor() {
    Track *track = mFirstTrack;
    while (track) {
        Track *next = track->next;

        delete track;
        track = next;
    }
    mFirstTrack = mLastTrack = NULL;

    SINF *sinf = mFirstSINF;
    while (sinf) {
        SINF *next = sinf->next;
        delete sinf->IPMPData;
        delete sinf;
        sinf = next;
    }
    mFirstSINF = NULL;
}

sp<MetaData> MoofExtractor::getMetaData() {
    status_t err;
    if ((err = readMetaData()) != OK) {
        return new MetaData;
    }

    ALOGV("getMetaData");

    return mFileMetaData;
}

uint32_t MoofExtractor::flags() const {
    return 0;
}

size_t MoofExtractor::countTracks() {
    status_t err;
    if ((err = readMetaData()) != OK) {
        return 0;
    }

    size_t n = 0;
    Track *track = mFirstTrack;
    while (track) {
        ++n;
        track = track->next;
    }

    //FIXME
    n = 1;
    ALOGV("countTracks %d", n);

    return n;
}

sp<MetaData> MoofExtractor::getTrackMetaData(
        size_t index, uint32_t flags) {
    status_t err;
    if ((err = readMetaData()) != OK) {
        return NULL;
    }

    ALOGV("getTrackMetaData %d %x", index, flags);
    Track *track = mFirstTrack;
    while (index > 0) {
        if (track == NULL) {
            return NULL;
        }

        track = track->next;
        --index;
    }

    if (track == NULL) {
        return NULL;
    }

    return track->meta;
}

status_t MoofExtractor::readMetaData()
{
    if (mInitCheck != NO_INIT) {
        return mInitCheck;
    }

    ALOGV("readMetaData");

    off64_t offset = 0;
    status_t err = 0;
    while ((err = parseChunk(&offset, 0)) == OK) {
    }
    ALOGD("metadata read offset=%lld/%llx err=%d", offset, offset, err);

    if (mInitCheck == OK) {
        if (mHasVideo) {
            mFileMetaData->setCString(
                    kKeyMIMEType, MEDIA_MIMETYPE_CONTAINER_DASH);
        } else {
            mFileMetaData->setCString(kKeyMIMEType, "audio/mp4");
        }

        mInitCheck = OK;
    } else {
        ALOGV("%s:%d: mInitCheck=%d", __FILE__, __LINE__, err);
        mInitCheck = err;
    }

    if (mInitCheck == OK)
        mOffset = offset;
    else
        mOffset = 0;

    mMoofSampleTable = new MoofSampleTable(mOffset, mDataSource);
    ALOGD("new mMoofSampleTable");
    Track *track = mFirstTrack;
    size_t trackIndex = 0;
    while (track) {
        mMoofSampleTable->setTimeScale(trackIndex, track->timescale);
        mMoofSampleTable->setTrackMetaData(trackIndex, track->meta);
        trackIndex++;
        track = track->next;
    }

    mMoofSampleTable->onMetaDataReady();

    CHECK_NE(err, (status_t)NO_INIT);
    return mInitCheck;
}

char* MoofExtractor::getDrmTrackInfo(size_t trackID, int *len) {
    if (mFirstSINF == NULL) {
        return NULL;
    }

    SINF *sinf = mFirstSINF;
    while (sinf && (trackID != sinf->trackID)) {
        sinf = sinf->next;
    }

    if (sinf == NULL) {
        return NULL;
    }

    *len = sinf->len;
    return sinf->IPMPData;
}

// Reads an encoded integer 7 bits at a time until it encounters the high bit clear.
int32_t MoofExtractor::readSize(off64_t offset,
                                const sp<DataSource> DataSource, uint8_t *numOfBytes) {
    uint32_t size = 0;
    uint8_t data;
    bool moreData = true;
    *numOfBytes = 0;

    while (moreData) {
        if (DataSource->readAt(offset, &data, 1) < 1) {
            return -1;
        }
        offset ++;
        moreData = (data >= 128) ? true : false;
        size = (size << 7) | (data & 0x7f); // Take last 7 bits
        (*numOfBytes) ++;
    }

    return size;
}

status_t MoofExtractor::parseDrmSINF(off64_t *offset, off64_t data_offset)
{
    uint8_t updateIdTag;
    if (mDataSource->readAt(data_offset, &updateIdTag, 1) < 1) {
        ALOGE("%s:%d: short read @ %lld", __FILE__, __LINE__, data_offset);
        return ERROR_IO;
    }
    data_offset ++;

    if (0x01/*OBJECT_DESCRIPTOR_UPDATE_ID_TAG*/ != updateIdTag) {
        return ERROR_MALFORMED;
    }

    uint8_t numOfBytes;
    int32_t size = readSize(data_offset, mDataSource, &numOfBytes);
    if (size < 0) {
        return ERROR_IO;
    }
    int32_t classSize = size;
    data_offset += numOfBytes;

    while(size >= 11 ) {
        uint8_t descriptorTag;
        if (mDataSource->readAt(data_offset, &descriptorTag, 1) < 1) {
            return ERROR_IO;
        }
        data_offset ++;

        if (0x11/*OBJECT_DESCRIPTOR_ID_TAG*/ != descriptorTag) {
            return ERROR_MALFORMED;
        }

        uint8_t buffer[8];
        //ObjectDescriptorID and ObjectDescriptor url flag
        if (mDataSource->readAt(data_offset, buffer, 2) < 2) {
            return ERROR_IO;
        }
        data_offset += 2;

        if ((buffer[1] >> 5) & 0x0001) { //url flag is set
            return ERROR_MALFORMED;
        }

        if (mDataSource->readAt(data_offset, buffer, 8) < 8) {
            return ERROR_IO;
        }
        data_offset += 8;

        if ((0x0F/*ES_ID_REF_TAG*/ != buffer[1])
                || ( 0x0A/*IPMP_DESCRIPTOR_POINTER_ID_TAG*/ != buffer[5])) {
            return ERROR_MALFORMED;
        }

        SINF *sinf = new SINF;
        sinf->trackID = U16_AT(&buffer[3]);
        sinf->IPMPDescriptorID = buffer[7];
        sinf->next = mFirstSINF;
        mFirstSINF = sinf;

        size -= (8 + 2 + 1);
    }

    if (size != 0) {
        return ERROR_MALFORMED;
    }

    if (mDataSource->readAt(data_offset, &updateIdTag, 1) < 1) {
        return ERROR_IO;
    }
    data_offset ++;

    if(0x05/*IPMP_DESCRIPTOR_UPDATE_ID_TAG*/ != updateIdTag) {
        return ERROR_MALFORMED;
    }

    size = readSize(data_offset, mDataSource, &numOfBytes);
    if (size < 0) {
        return ERROR_IO;
    }
    classSize = size;
    data_offset += numOfBytes;

    while (size > 0) {
        uint8_t tag;
        int32_t dataLen;
        if (mDataSource->readAt(data_offset, &tag, 1) < 1) {
            return ERROR_IO;
        }
        data_offset ++;

        if (0x0B/*IPMP_DESCRIPTOR_ID_TAG*/ == tag) {
            uint8_t id;
            dataLen = readSize(data_offset, mDataSource, &numOfBytes);
            if (dataLen < 0) {
                return ERROR_IO;
            } else if (dataLen < 4) {
                return ERROR_MALFORMED;
            }
            data_offset += numOfBytes;

            if (mDataSource->readAt(data_offset, &id, 1) < 1) {
                return ERROR_IO;
            }
            data_offset ++;

            SINF *sinf = mFirstSINF;
            while (sinf && (sinf->IPMPDescriptorID != id)) {
                sinf = sinf->next;
            }
            if (sinf == NULL) {
                return ERROR_MALFORMED;
            }
            sinf->len = dataLen - 3;
            sinf->IPMPData = new char[sinf->len];

            if (mDataSource->readAt(data_offset + 2, sinf->IPMPData, sinf->len) < sinf->len) {
                return ERROR_IO;
            }
            data_offset += sinf->len;

            size -= (dataLen + numOfBytes + 1);
        }
    }

    if (size != 0) {
        return ERROR_MALFORMED;
    }

    return UNKNOWN_ERROR;  // Return a dummy error.
}

struct PathAdder {
    PathAdder(Vector<uint32_t> *path, uint32_t chunkType)
        : mPath(path) {
        mPath->push(chunkType);
    }

    ~PathAdder() {
        mPath->pop();
    }

private:
    Vector<uint32_t> *mPath;

    PathAdder(const PathAdder &);
    PathAdder &operator=(const PathAdder &);
};

static bool underMetaDataPath(const Vector<uint32_t> &path) {
    return path.size() >= 5
        && path[0] == FOURCC('m', 'o', 'o', 'v')
        && path[1] == FOURCC('u', 'd', 't', 'a')
        && path[2] == FOURCC('m', 'e', 't', 'a')
        && path[3] == FOURCC('i', 'l', 's', 't');
}

// Given a time in seconds since Jan 1 1904, produce a human-readable string.
static void convertTimeToDate(int64_t time_1904, String8 *s) {
    time_t time_1970 = time_1904 - (((66 * 365 + 17) * 24) * 3600);

    char tmp[32];
    strftime(tmp, sizeof(tmp), "%Y%m%dT%H%M%S.000Z", gmtime(&time_1970));

    s->setTo(tmp);
}

status_t MoofExtractor::parseChunk(off64_t *offset, int depth) {
    ALOGV("entering parseChunk %lld/%d", *offset, depth);
    uint32_t hdr[2];
    if (mDataSource->readAt(*offset, hdr, 8) < 8) {
        return ERROR_IO;
    }
    uint64_t chunk_size = ntohl(hdr[0]);
    uint32_t chunk_type = ntohl(hdr[1]);
    off64_t data_offset = *offset + 8;

    if (chunk_size == 1) {
        if (mDataSource->readAt(*offset + 8, &chunk_size, 8) < 8) {
            return ERROR_IO;
        }
        chunk_size = ntoh64(chunk_size);
        data_offset += 8;

        if (chunk_size < 16) {
            // The smallest valid chunk is 16 bytes long in this case.
            return ERROR_MALFORMED;
        }
    } else if (chunk_size < 8) {
        // The smallest valid chunk is 8 bytes long.
        return ERROR_MALFORMED;
    }

    char chunk[5];
    MakeFourCCString(chunk_type, chunk);
    ALOGV("chunk: %s @ %lld", chunk, *offset);

    PathAdder autoAdder(&mPath, chunk_type);

    off64_t chunk_data_size = *offset + chunk_size - data_offset;

    if (chunk_type != FOURCC('c', 'p', 'r', 't')
            && chunk_type != FOURCC('c', 'o', 'v', 'r')
            && mPath.size() == 5 && underMetaDataPath(mPath)) {
        off64_t stop_offset = *offset + chunk_size;
        *offset = data_offset;
        while (*offset < stop_offset) {
            status_t err = parseChunk(offset, depth + 1);
            if (err != OK) {
                return err;
            }
        }

        if (*offset != stop_offset) {
            return ERROR_MALFORMED;
        }

        return OK;
    }

    switch(chunk_type) {
        case FOURCC('f', 't', 'y', 'p'):
        case FOURCC('f', 'r', 'e', 'e'):
            mDataSource->setCachedRange(data_offset, chunk_size-data_offset);
            *offset += chunk_size;
        break;
        case FOURCC('m', 'o', 'o', 'v'):
            mDataSource->setCachedRange(data_offset, chunk_size-data_offset);
            // fll through
        case FOURCC('t', 'r', 'a', 'k'):
        case FOURCC('m', 'd', 'i', 'a'):
        case FOURCC('m', 'i', 'n', 'f'):
        case FOURCC('d', 'i', 'n', 'f'):
        case FOURCC('s', 't', 'b', 'l'):
        case FOURCC('m', 'v', 'e', 'x'):
        case FOURCC('t', 'r', 'a', 'f'):
        case FOURCC('m', 'f', 'r', 'a'):
        case FOURCC('u', 'd', 't', 'a'):
        case FOURCC('i', 'l', 's', 't'):
        {
            if (chunk_type == FOURCC('s', 't', 'b', 'l')) {
                ALOGV("sampleTable chunk is %d bytes long.", (size_t)chunk_size);

                if (mDataSource->flags()
                        & (DataSource::kWantsPrefetching
                            | DataSource::kIsCachingDataSource)) {

                    mDataSource->setCachedRange(*offset, chunk_size);
                }

            }

            bool isTrack = false;
            if (chunk_type == FOURCC('t', 'r', 'a', 'k')) {
                isTrack = true;

                Track *track = new Track;
                track->next = NULL;
                if (mLastTrack) {
                    mLastTrack->next = track;
                } else {
                    mFirstTrack = track;
                }
                mLastTrack = track;
                mNumTracks++;

                track->meta = new MetaData;
                track->includes_expensive_metadata = false;
                track->skipTrack = false;
                track->timescale = 0;
                track->meta->setCString(kKeyMIMEType, "application/octet-stream");

            }

            off64_t stop_offset = *offset + chunk_size;
            *offset = data_offset;
            while (*offset < stop_offset) {
                status_t err = parseChunk(offset, depth + 1);
                if (err != OK) {
                    return err;
                }
            }

            if (*offset != stop_offset) {
                return ERROR_MALFORMED;
            }

            if (isTrack) {
                if (mLastTrack->skipTrack) {
                    Track *cur = mFirstTrack;

                    if (cur == mLastTrack) {
                        delete cur;
                        mFirstTrack = mLastTrack = NULL;
                    } else {
                        while (cur && cur->next != mLastTrack) {
                            cur = cur->next;
                        }
                        cur->next = NULL;
                        delete mLastTrack;
                        mLastTrack = cur;
                    }

                    return OK;
                }

                status_t err = verifyTrack(mLastTrack);

                if (err != OK) {
                    return err;
                }
            } else if (chunk_type == FOURCC('m', 'o', 'o', 'v')) {
                mInitCheck = OK;

                if (!mIsDrm) {
                    return UNKNOWN_ERROR;  // Return a dummy error.
                } else {
                    return OK;
                }
            }
            break;
        }

        case FOURCC('t', 'k', 'h', 'd'):
        {
            status_t err;
            if ((err = parseTrackHeader(data_offset, chunk_data_size)) != OK) {
                return err;
            }

            *offset += chunk_size;
            break;
        }

        case FOURCC('m', 'd', 'h', 'd'):
        {
            if (chunk_data_size < 4) {
                return ERROR_MALFORMED;
            }

            uint8_t version;
            if (mDataSource->readAt(
                        data_offset, &version, sizeof(version))
                    < (ssize_t)sizeof(version)) {
                return ERROR_IO;
            }

            off64_t timescale_offset;

            if (version == 1) {
                timescale_offset = data_offset + 4 + 16;
            } else if (version == 0) {
                timescale_offset = data_offset + 4 + 8;
            } else {
                return ERROR_IO;
            }

            uint32_t timescale;
            if (mDataSource->readAt(
                        timescale_offset, &timescale, sizeof(timescale))
                    < (ssize_t)sizeof(timescale)) {
                return ERROR_IO;
            }

            mLastTrack->timescale = ntohl(timescale);
            ALOGV("track %p timescale %d", mLastTrack, mLastTrack->timescale);

            int64_t duration;
            if (version == 1) {
                if (mDataSource->readAt(
                            timescale_offset + 4, &duration, sizeof(duration))
                        < (ssize_t)sizeof(duration)) {
                    return ERROR_IO;
                }
                duration = ntoh64(duration);
            } else {
                int32_t duration32;
                if (mDataSource->readAt(
                            timescale_offset + 4, &duration32, sizeof(duration32))
                        < (ssize_t)sizeof(duration32)) {
                    return ERROR_IO;
                }
                duration = ntohl(duration32);
            }
            mLastTrack->meta->setInt64(
                    kKeyDuration, (duration * 1000000) / mLastTrack->timescale);

            uint8_t lang[2];
            off64_t lang_offset;
            if (version == 1) {
                lang_offset = timescale_offset + 4 + 8;
            } else if (version == 0) {
                lang_offset = timescale_offset + 4 + 4;
            } else {
                return ERROR_IO;
            }

            if (mDataSource->readAt(lang_offset, &lang, sizeof(lang))
                    < (ssize_t)sizeof(lang)) {
                return ERROR_IO;
            }

            // To get the ISO-639-2/T three character language code
            // 1 bit pad followed by 3 5-bits characters. Each character
            // is packed as the difference between its ASCII value and 0x60.
            char lang_code[4];
            lang_code[0] = ((lang[0] >> 2) & 0x1f) + 0x60;
            lang_code[1] = ((lang[0] & 0x3) << 3 | (lang[1] >> 5)) + 0x60;
            lang_code[2] = (lang[1] & 0x1f) + 0x60;
            lang_code[3] = '\0';

            mLastTrack->meta->setCString(
                    kKeyMediaLanguage, lang_code);

            *offset += chunk_size;
            break;
        }

        case FOURCC('s', 't', 's', 'd'):
        {
            if (chunk_data_size < 8) {
                return ERROR_MALFORMED;
            }

            uint8_t buffer[8];
            if (chunk_data_size < (off64_t)sizeof(buffer)) {
                return ERROR_MALFORMED;
            }

            if (mDataSource->readAt(
                        data_offset, buffer, 8) < 8) {
                return ERROR_IO;
            }

            if (U32_AT(buffer) != 0) {
                // Should be version 0, flags 0.
                return ERROR_MALFORMED;
            }

            uint32_t entry_count = U32_AT(&buffer[4]);

            if (entry_count > 1) {
                // For 3GPP timed text, there could be multiple tx3g boxes contain
                // multiple text display formats. These formats will be used to
                // display the timed text.
                const char *mime;
                CHECK(mLastTrack->meta->findCString(kKeyMIMEType, &mime));
                if (strcasecmp(mime, MEDIA_MIMETYPE_TEXT_3GPP)) {
                    // For now we only support a single type of media per track.
                    mLastTrack->skipTrack = true;
                    *offset += chunk_size;
                    break;
                }
            }

            off64_t stop_offset = *offset + chunk_size;
            *offset = data_offset + 8;
            for (uint32_t i = 0; i < entry_count; ++i) {
                status_t err = parseChunk(offset, depth + 1);
                if (err != OK) {
                    return err;
                }
            }

            if (*offset != stop_offset) {
                return ERROR_MALFORMED;
            }
            break;
        }

        case FOURCC('m', 'p', '4', 'a'):
        case FOURCC('s', 'a', 'm', 'r'):
        case FOURCC('s', 'a', 'w', 'b'):
        {
            uint8_t buffer[8 + 20];
            if (chunk_data_size < (ssize_t)sizeof(buffer)) {
                // Basic AudioSampleEntry size.
                return ERROR_MALFORMED;
            }

            if (mDataSource->readAt(
                        data_offset, buffer, sizeof(buffer)) < (ssize_t)sizeof(buffer)) {
                return ERROR_IO;
            }

            uint16_t data_ref_index = U16_AT(&buffer[6]);
            uint16_t num_channels = U16_AT(&buffer[16]);

            uint16_t sample_size = U16_AT(&buffer[18]);
            uint32_t sample_rate = U32_AT(&buffer[24]) >> 16;

            if (!strcasecmp(MEDIA_MIMETYPE_AUDIO_AMR_NB,
                            FourCC2MIME(chunk_type))) {
                // AMR NB audio is always mono, 8kHz
                num_channels = 1;
                sample_rate = 8000;
            } else if (!strcasecmp(MEDIA_MIMETYPE_AUDIO_AMR_WB,
                               FourCC2MIME(chunk_type))) {
                // AMR WB audio is always mono, 16kHz
                num_channels = 1;
                sample_rate = 16000;
            }

#if 0
            printf("*** coding='%s' %d channels, size %d, rate %d\n",
                   chunk, num_channels, sample_size, sample_rate);
#endif

            mLastTrack->meta->setCString(kKeyMIMEType, FourCC2MIME(chunk_type));
            mLastTrack->meta->setInt32(kKeyChannelCount, num_channels);
            mLastTrack->meta->setInt32(kKeySampleRate, sample_rate);

            off64_t stop_offset = *offset + chunk_size;
            *offset = data_offset + sizeof(buffer);
            while (*offset < stop_offset) {
                status_t err = parseChunk(offset, depth + 1);
                if (err != OK) {
                    return err;
                }
            }

            if (*offset != stop_offset) {
                return ERROR_MALFORMED;
            }
            break;
        }

        case FOURCC('m', 'p', '4', 'v'):
        case FOURCC('s', '2', '6', '3'):
        case FOURCC('H', '2', '6', '3'):
        case FOURCC('h', '2', '6', '3'):
        case FOURCC('a', 'v', 'c', '1'):
        {
            mHasVideo = true;

            uint8_t buffer[78];
            if (chunk_data_size < (ssize_t)sizeof(buffer)) {
                // Basic VideoSampleEntry size.
                return ERROR_MALFORMED;
            }

            if (mDataSource->readAt(
                        data_offset, buffer, sizeof(buffer)) < (ssize_t)sizeof(buffer)) {
                return ERROR_IO;
            }

            uint16_t data_ref_index = U16_AT(&buffer[6]);
            uint16_t width = U16_AT(&buffer[6 + 18]);
            uint16_t height = U16_AT(&buffer[6 + 20]);

            // The video sample is not stand-compliant if it has invalid dimension.
            // Use some default width and height value, and
            // let the decoder figure out the actual width and height (and thus
            // be prepared for INFO_FOMRAT_CHANGED event).
            if (width == 0)  width  = 352;
            if (height == 0) height = 288;

            ALOGV("*** coding='%s' width=%d height=%d\n",
                  chunk, width, height);

            mLastTrack->meta->setCString(kKeyMIMEType, FourCC2MIME(chunk_type));
            mLastTrack->meta->setInt32(kKeyWidth, width);
            mLastTrack->meta->setInt32(kKeyHeight, height);

            off64_t stop_offset = *offset + chunk_size;
            *offset = data_offset + sizeof(buffer);
            while (*offset < stop_offset) {
                status_t err = parseChunk(offset, depth + 1);
                if (err != OK) {
                    return err;
                }
            }

            if (*offset != stop_offset) {
                return ERROR_MALFORMED;
            }
            break;
        }

        case FOURCC('s', 't', 'c', 'o'):
        case FOURCC('c', 'o', '6', '4'):
        {
            // this box should be empty -- per fragment metadata in a moof box
            *offset += chunk_size;
            break;
        }

        case FOURCC('s', 't', 's', 'c'):
        {
            // this box should be empty -- per fragment metadata in a moof box

            *offset += chunk_size;
            break;
        }

        case FOURCC('s', 't', 's', 'z'):
        case FOURCC('s', 't', 'z', '2'):
        {
            // this box should be empty -- per fragment metadata in a moof box

            size_t max_size = 200000;
            mLastTrack->meta->setInt32(kKeyMaxInputSize, max_size);
            *offset += chunk_size;

            // Calculate average frame rate.
            const char *mime;
            CHECK(mLastTrack->meta->findCString(kKeyMIMEType, &mime));
            if (!strncasecmp("video/", mime, 6)) {
                size_t nSamples = 0;
                int64_t durationUs;
                if (mLastTrack->meta->findInt64(kKeyDuration, &durationUs)) {
                    if (durationUs > 0) {
                        int32_t frameRate = (nSamples * 1000000LL +
                                    (durationUs >> 1)) / durationUs;
                        mLastTrack->meta->setInt32(kKeyFrameRate, frameRate);
                    }
                }
            }

            break;
        }

        case FOURCC('s', 't', 't', 's'):
        {
            // this box should be empty -- per fragment metadata in a moof box

            *offset += chunk_size;
            break;
        }

        case FOURCC('c', 't', 't', 's'):
        {
            // this box should be empty -- per fragment metadata in a moof box

            *offset += chunk_size;
            break;
        }

        case FOURCC('s', 't', 's', 's'):
        {
            // this box should be empty -- per fragment metadata in a moof box

            *offset += chunk_size;
            break;
        }

        // @xyz
        case FOURCC('\xA9', 'x', 'y', 'z'):
        {
            // Best case the total data length inside "@xyz" box
            // would be 8, for instance "@xyz" + "\x00\x04\x15\xc7" + "0+0/",
            // where "\x00\x04" is the text string length with value = 4,
            // "\0x15\xc7" is the language code = en, and "0+0" is a
            // location (string) value with longitude = 0 and latitude = 0.
            if (chunk_data_size < 8) {
                return ERROR_MALFORMED;
            }

            // Worst case the location string length would be 18,
            // for instance +90.0000-180.0000, without the trailing "/" and
            // the string length + language code.
            char buffer[18];

            // Substracting 5 from the data size is because the text string length +
            // language code takes 4 bytes, and the trailing slash "/" takes 1 byte.
            off64_t location_length = chunk_data_size - 5;
            if (location_length >= (off64_t) sizeof(buffer)) {
                return ERROR_MALFORMED;
            }

            if (mDataSource->readAt(
                        data_offset + 4, buffer, location_length) < location_length) {
                return ERROR_IO;
            }

            buffer[location_length] = '\0';
            mFileMetaData->setCString(kKeyLocation, buffer);
            *offset += chunk_size;
            break;
        }

        case FOURCC('e', 's', 'd', 's'):
        {
            if (chunk_data_size < 4) {
                return ERROR_MALFORMED;
            }

            uint8_t buffer[256];
            if (chunk_data_size > (off64_t)sizeof(buffer)) {
                return ERROR_BUFFER_TOO_SMALL;
            }

            if (mDataSource->readAt(
                        data_offset, buffer, chunk_data_size) < chunk_data_size) {
                return ERROR_IO;
            }

            if (U32_AT(buffer) != 0) {
                // Should be version 0, flags 0.
                return ERROR_MALFORMED;
            }

            ALOGV("esds data size %lld", chunk_data_size - 4);
            mLastTrack->meta->setData(
                    kKeyESDS, kTypeESDS, &buffer[4], chunk_data_size - 4);

            if (mPath.size() >= 2
                    && mPath[mPath.size() - 2] == FOURCC('m', 'p', '4', 'a')) {
                // Information from the ESDS must be relied on for proper
                // setup of sample rate and channel count for Moof Audio.
                // The generic header appears to only contain generic
                // information...

                status_t err = updateAudioTrackInfoFromESDS_MPEG4Audio(
                        &buffer[4], chunk_data_size - 4);

                if (err != OK) {
                    return err;
                }
            }

            *offset += chunk_size;
            break;
        }

        case FOURCC('a', 'v', 'c', 'C'):
        {
            char buffer[256];
            if (chunk_data_size > (off64_t)sizeof(buffer)) {
                return ERROR_BUFFER_TOO_SMALL;
            }

            if (mDataSource->readAt(
                        data_offset, buffer, chunk_data_size) < chunk_data_size) {
                return ERROR_IO;
            }

            mLastTrack->meta->setData(
                    kKeyAVCC, kTypeAVCC, buffer, chunk_data_size);

            *offset += chunk_size;
            break;
        }

        case FOURCC('d', '2', '6', '3'):
        {
            /*
             * d263 contains a fixed 7 bytes part:
             *   vendor - 4 bytes
             *   version - 1 byte
             *   level - 1 byte
             *   profile - 1 byte
             * optionally, "d263" box itself may contain a 16-byte
             * bit rate box (bitr)
             *   average bit rate - 4 bytes
             *   max bit rate - 4 bytes
             */
            char buffer[23];
            if (chunk_data_size != 7 &&
                chunk_data_size != 23) {
                ALOGE("Incorrect D263 box size %lld", chunk_data_size);
                return ERROR_MALFORMED;
            }

            if (mDataSource->readAt(
                    data_offset, buffer, chunk_data_size) < chunk_data_size) {
                return ERROR_IO;
            }

            mLastTrack->meta->setData(kKeyD263, kTypeD263, buffer, chunk_data_size);

            *offset += chunk_size;
            break;
        }

        case FOURCC('m', 'e', 't', 'a'):
        {
            uint8_t buffer[4];
            if (chunk_data_size < (off64_t)sizeof(buffer)) {
                return ERROR_MALFORMED;
            }

            if (mDataSource->readAt(
                        data_offset, buffer, 4) < 4) {
                return ERROR_IO;
            }

            if (U32_AT(buffer) != 0) {
                // Should be version 0, flags 0.

                // If it's not, let's assume this is one of those
                // apparently malformed chunks that don't have flags
                // and completely different semantics than what's
                // in the Moof specs and skip it.
                *offset += chunk_size;
                return OK;
            }

            off64_t stop_offset = *offset + chunk_size;
            *offset = data_offset + sizeof(buffer);
            while (*offset < stop_offset) {
                status_t err = parseChunk(offset, depth + 1);
                if (err != OK) {
                    return err;
                }
            }

            if (*offset != stop_offset) {
                return ERROR_MALFORMED;
            }
            break;
        }

        case FOURCC('m', 'e', 'a', 'n'):
        case FOURCC('n', 'a', 'm', 'e'):
        case FOURCC('d', 'a', 't', 'a'):
        {
            if (mPath.size() == 6 && underMetaDataPath(mPath)) {
                status_t err = parseMetaData(data_offset, chunk_data_size);

                if (err != OK) {
                    return err;
                }
            }

            *offset += chunk_size;
            break;
        }

        case FOURCC('m', 'v', 'h', 'd'):
        {
            if (chunk_data_size < 12) {
                return ERROR_MALFORMED;
            }

            uint8_t header[12];
            if (mDataSource->readAt(
                        data_offset, header, sizeof(header))
                    < (ssize_t)sizeof(header)) {
                return ERROR_IO;
            }

            int64_t creationTime;
            if (header[0] == 1) {
                creationTime = U64_AT(&header[4]);
            } else if (header[0] != 0) {
                return ERROR_MALFORMED;
            } else {
                creationTime = U32_AT(&header[4]);
            }

            String8 s;
            convertTimeToDate(creationTime, &s);

            mFileMetaData->setCString(kKeyDate, s.string());

            *offset += chunk_size;
            break;
        }

        case FOURCC('m', 'd', 'a', 't'):
        {
            if (!mIsDrm) {
                *offset += chunk_size;
                break;
            }

            if (chunk_size < 8) {
                return ERROR_MALFORMED;
            }

            return parseDrmSINF(offset, data_offset);
        }

        case FOURCC('h', 'd', 'l', 'r'):
        {
            uint32_t buffer;
            if (mDataSource->readAt(
                        data_offset + 8, &buffer, 4) < 4) {
                return ERROR_IO;
            }

            uint32_t type = ntohl(buffer);
            // For the 3GPP file format, the handler-type within the 'hdlr' box
            // shall be 'text'. We also want to support 'sbtl' handler type
            // for a practical reason as various Moof containers use it.
            if (type == FOURCC('t', 'e', 'x', 't') || type == FOURCC('s', 'b', 't', 'l')) {
                mLastTrack->meta->setCString(kKeyMIMEType, MEDIA_MIMETYPE_TEXT_3GPP);
            }

            *offset += chunk_size;
            break;
        }

        case FOURCC('t', 'x', '3', 'g'):
        {
            uint32_t type;
            const void *data;
            size_t size = 0;
            if (!mLastTrack->meta->findData(
                    kKeyTextFormatData, &type, &data, &size)) {
                size = 0;
            }

            uint8_t *buffer = new uint8_t[size + chunk_size];

            if (size > 0) {
                memcpy(buffer, data, size);
            }

            if ((size_t)(mDataSource->readAt(*offset, buffer + size, chunk_size))
                    < chunk_size) {
                delete[] buffer;
                buffer = NULL;

                return ERROR_IO;
            }

            mLastTrack->meta->setData(
                    kKeyTextFormatData, 0, buffer, size + chunk_size);

            delete[] buffer;

            *offset += chunk_size;
            break;
        }

        case FOURCC('c', 'o', 'v', 'r'):
        {
            if (mFileMetaData != NULL) {
                ALOGV("chunk_data_size = %lld and data_offset = %lld",
                        chunk_data_size, data_offset);
                uint8_t *buffer = new uint8_t[chunk_data_size + 1];
                if (mDataSource->readAt(
                    data_offset, buffer, chunk_data_size) != (ssize_t)chunk_data_size) {
                    delete[] buffer;
                    buffer = NULL;

                    return ERROR_IO;
                }
                const int kSkipBytesOfDataBox = 16;
                mFileMetaData->setData(
                    kKeyAlbumArt, MetaData::TYPE_NONE,
                    buffer + kSkipBytesOfDataBox, chunk_data_size - kSkipBytesOfDataBox);
            }

            *offset += chunk_size;
            break;
        }

        case FOURCC('-', '-', '-', '-'):
        {
            mLastCommentMean.clear();
            mLastCommentName.clear();
            mLastCommentData.clear();
            *offset += chunk_size;
            break;
        }

        default:
        {
            ALOGV("default chunk size %x %lld", chunk_type, chunk_size);
            *offset += chunk_size;
            break;
        }
    }

    return OK;
}

status_t MoofExtractor::parseTrackHeader(
        off64_t data_offset, off64_t data_size) {
    if (data_size < 4) {
        return ERROR_MALFORMED;
    }

    ALOGV("parseTrackHeader");
    uint8_t version;
    if (mDataSource->readAt(data_offset, &version, 1) < 1) {
        return ERROR_IO;
    }

    size_t dynSize = (version == 1) ? 36 : 24;

    uint8_t buffer[36 + 60];

    if (data_size != (off64_t)dynSize + 60) {
        return ERROR_MALFORMED;
    }

    if (mDataSource->readAt(
                data_offset, buffer, data_size) < (ssize_t)data_size) {
        return ERROR_IO;
    }

    uint64_t ctime, mtime, duration;
    int32_t id;

    if (version == 1) {
        ctime = U64_AT(&buffer[4]);
        mtime = U64_AT(&buffer[12]);
        id = U32_AT(&buffer[20]);
        duration = U64_AT(&buffer[28]);
    } else {
        CHECK_EQ((unsigned)version, 0u);

        ctime = U32_AT(&buffer[4]);
        mtime = U32_AT(&buffer[8]);
        id = U32_AT(&buffer[12]);
        duration = U32_AT(&buffer[20]);
    }

    mLastTrack->meta->setInt32(kKeyTrackID, id);

    size_t matrixOffset = dynSize + 16;
    int32_t a00 = U32_AT(&buffer[matrixOffset]);
    int32_t a01 = U32_AT(&buffer[matrixOffset + 4]);
    int32_t dx = U32_AT(&buffer[matrixOffset + 8]);
    int32_t a10 = U32_AT(&buffer[matrixOffset + 12]);
    int32_t a11 = U32_AT(&buffer[matrixOffset + 16]);
    int32_t dy = U32_AT(&buffer[matrixOffset + 20]);

#if 0
    ALOGI("x' = %.2f * x + %.2f * y + %.2f",
         a00 / 65536.0f, a01 / 65536.0f, dx / 65536.0f);
    ALOGI("y' = %.2f * x + %.2f * y + %.2f",
         a10 / 65536.0f, a11 / 65536.0f, dy / 65536.0f);
#endif

    uint32_t rotationDegrees;

    static const int32_t kFixedOne = 0x10000;
    if (a00 == kFixedOne && a01 == 0 && a10 == 0 && a11 == kFixedOne) {
        // Identity, no rotation
        rotationDegrees = 0;
    } else if (a00 == 0 && a01 == kFixedOne && a10 == -kFixedOne && a11 == 0) {
        rotationDegrees = 90;
    } else if (a00 == 0 && a01 == -kFixedOne && a10 == kFixedOne && a11 == 0) {
        rotationDegrees = 270;
    } else if (a00 == -kFixedOne && a01 == 0 && a10 == 0 && a11 == -kFixedOne) {
        rotationDegrees = 180;
    } else {
        ALOGW("We only support 0,90,180,270 degree rotation matrices");
        rotationDegrees = 0;
    }

    if (rotationDegrees != 0) {
        mLastTrack->meta->setInt32(kKeyRotation, rotationDegrees);
    }

    // Handle presentation display size, which could be different
    // from the image size indicated by kKeyWidth and kKeyHeight.
    uint32_t width = U32_AT(&buffer[dynSize + 52]);
    uint32_t height = U32_AT(&buffer[dynSize + 56]);
    mLastTrack->meta->setInt32(kKeyDisplayWidth, width >> 16);
    mLastTrack->meta->setInt32(kKeyDisplayHeight, height >> 16);

    return OK;
}

status_t MoofExtractor::parseMetaData(off64_t offset, size_t size) {
    if (size < 4) {
        return ERROR_MALFORMED;
    }

    ALOGV("parseMetaData");
    uint8_t *buffer = new uint8_t[size + 1];
    if (mDataSource->readAt(
                offset, buffer, size) != (ssize_t)size) {
        delete[] buffer;
        buffer = NULL;

        return ERROR_IO;
    }

    uint32_t flags = U32_AT(buffer);

    uint32_t metadataKey = 0;
    char chunk[5];
    MakeFourCCString(mPath[4], chunk);
    ALOGV("meta: %s @ %lld", chunk, offset);
    switch (mPath[4]) {
        case FOURCC(0xa9, 'a', 'l', 'b'):
        {
            metadataKey = kKeyAlbum;
            break;
        }
        case FOURCC(0xa9, 'A', 'R', 'T'):
        {
            metadataKey = kKeyArtist;
            break;
        }
        case FOURCC('a', 'A', 'R', 'T'):
        {
            metadataKey = kKeyAlbumArtist;
            break;
        }
        case FOURCC(0xa9, 'd', 'a', 'y'):
        {
            metadataKey = kKeyYear;
            break;
        }
        case FOURCC(0xa9, 'n', 'a', 'm'):
        {
            metadataKey = kKeyTitle;
            break;
        }
        case FOURCC(0xa9, 'w', 'r', 't'):
        {
            metadataKey = kKeyWriter;
            break;
        }
        case FOURCC('c', 'o', 'v', 'r'):
        {
            metadataKey = kKeyAlbumArt;
            break;
        }
        case FOURCC('g', 'n', 'r', 'e'):
        {
            metadataKey = kKeyGenre;
            break;
        }
        case FOURCC(0xa9, 'g', 'e', 'n'):
        {
            metadataKey = kKeyGenre;
            break;
        }
        case FOURCC('c', 'p', 'i', 'l'):
        {
            if (size == 9 && flags == 21) {
                char tmp[16];
                sprintf(tmp, "%d",
                        (int)buffer[size - 1]);

                mFileMetaData->setCString(kKeyCompilation, tmp);
            }
            break;
        }
        case FOURCC('t', 'r', 'k', 'n'):
        {
            if (size == 16 && flags == 0) {
                char tmp[16];
                sprintf(tmp, "%d/%d",
                        (int)buffer[size - 5], (int)buffer[size - 3]);

                mFileMetaData->setCString(kKeyCDTrackNumber, tmp);
            }
            break;
        }
        case FOURCC('d', 'i', 's', 'k'):
        {
            if (size == 14 && flags == 0) {
                char tmp[16];
                sprintf(tmp, "%d/%d",
                        (int)buffer[size - 3], (int)buffer[size - 1]);

                mFileMetaData->setCString(kKeyDiscNumber, tmp);
            }
            break;
        }
        case FOURCC('-', '-', '-', '-'):
        {
            buffer[size] = '\0';
            switch (mPath[5]) {
                case FOURCC('m', 'e', 'a', 'n'):
                    mLastCommentMean.setTo((const char *)buffer + 4);
                    break;
                case FOURCC('n', 'a', 'm', 'e'):
                    mLastCommentName.setTo((const char *)buffer + 4);
                    break;
                case FOURCC('d', 'a', 't', 'a'):
                    mLastCommentData.setTo((const char *)buffer + 8);
                    break;
            }
            if (mLastCommentMean == "com.apple.iTunes"
                    && mLastCommentName == "iTunSMPB"
                    && mLastCommentData.length() != 0) {
                int32_t delay, padding;
                if (sscanf(mLastCommentData,
                           " %*x %x %x %*x", &delay, &padding) == 2) {
                    mLastTrack->meta->setInt32(kKeyEncoderDelay, delay);
                    mLastTrack->meta->setInt32(kKeyEncoderPadding, padding);
                }
                mLastCommentMean.clear();
                mLastCommentName.clear();
                mLastCommentData.clear();
            }
            break;
        }

        default:
            break;
    }

    if (size >= 8 && metadataKey) {
        if (metadataKey == kKeyAlbumArt) {
            mFileMetaData->setData(
                    kKeyAlbumArt, MetaData::TYPE_NONE,
                    buffer + 8, size - 8);
        } else if (metadataKey == kKeyGenre) {
            if (flags == 0) {
                // uint8_t genre code, iTunes genre codes are
                // the standard id3 codes, except they start
                // at 1 instead of 0 (e.g. Pop is 14, not 13)
                // We use standard id3 numbering, so subtract 1.
                int genrecode = (int)buffer[size - 1];
                genrecode--;
                if (genrecode < 0) {
                    genrecode = 255; // reserved for 'unknown genre'
                }
                char genre[10];
                sprintf(genre, "%d", genrecode);

                mFileMetaData->setCString(metadataKey, genre);
            } else if (flags == 1) {
                // custom genre string
                buffer[size] = '\0';

                mFileMetaData->setCString(
                        metadataKey, (const char *)buffer + 8);
            }
        } else {
            buffer[size] = '\0';

            mFileMetaData->setCString(
                    metadataKey, (const char *)buffer + 8);
        }
    }

    delete[] buffer;
    buffer = NULL;

    return OK;
}

sp<MediaSource> MoofExtractor::getTrack(size_t index) {
    status_t err;
    if ((err = readMetaData()) != OK) {
        return NULL;
    }

    ALOGV("getTrack %d", index);
    Track *track = mFirstTrack;
    int i = index;
    while (i > 0) {
        if (track == NULL) {
            ALOGD("null track %d", i);
            return NULL;
        }

        track = track->next;
        --i;
    }

    if (track == NULL) {
        ALOGD("no matching track %d", index);
        return NULL;
    }

    ALOGV("new MoofSource track %d %p timescale track->timescale %d", index, track, track->timescale);
    MoofSource *moofSource = new MoofSource(index, track->meta, mDataSource, track->timescale, mMoofSampleTable);
    ALOGV("moofSource %p", moofSource);
    return moofSource;
}

// static
status_t MoofExtractor::verifyTrack(Track *track) {
    const char *mime;
    CHECK(track->meta->findCString(kKeyMIMEType, &mime));

    uint32_t type;
    const void *data;
    size_t size;
    if (!strcasecmp(mime, MEDIA_MIMETYPE_VIDEO_AVC)) {
        if (!track->meta->findData(kKeyAVCC, &type, &data, &size)
                || type != kTypeAVCC) {
            return ERROR_MALFORMED;
        }
    } else if (!strcasecmp(mime, MEDIA_MIMETYPE_VIDEO_MPEG4)
            || !strcasecmp(mime, MEDIA_MIMETYPE_AUDIO_AAC)) {
        if (!track->meta->findData(kKeyESDS, &type, &data, &size)
                || type != kTypeESDS) {
            return ERROR_MALFORMED;
        }
    }

    return OK;
}

status_t MoofExtractor::updateAudioTrackInfoFromESDS_MPEG4Audio(
        const void *esds_data, size_t esds_size) {
    ESDS esds(esds_data, esds_size);

    uint8_t objectTypeIndication;
    if (esds.getObjectTypeIndication(&objectTypeIndication) != OK) {
        return ERROR_MALFORMED;
    }

    if (objectTypeIndication == 0xe1) {
        // This isn't Moof audio at all, it's QCELP 14k...
        mLastTrack->meta->setCString(kKeyMIMEType, MEDIA_MIMETYPE_AUDIO_QCELP);
        return OK;
    }

    if (objectTypeIndication  == 0x6b) {
        // The media subtype is MP3 audio
        // Our software MP3 audio decoder may not be able to handle
        // packetized MP3 audio; for now, lets just return ERROR_UNSUPPORTED
        ALOGE("MP3 track in MP4/3GPP file is not supported");
        return ERROR_UNSUPPORTED;
    }

    const uint8_t *csd;
    size_t csd_size;
    if (esds.getCodecSpecificInfo(
                (const void **)&csd, &csd_size) != OK) {
        return ERROR_MALFORMED;
    }

#if 0
    printf("ESD of size %d\n", csd_size);
    hexdump(csd, csd_size);
#endif

    if (csd_size == 0) {
        // There's no further information, i.e. no codec specific data
        // Let's assume that the information provided in the mpeg4 headers
        // is accurate and hope for the best.

        return OK;
    }

    if (csd_size < 2) {
        return ERROR_MALFORMED;
    }

    ABitReader br(csd, csd_size);
    uint32_t objectType = br.getBits(5);

    if (objectType == 31) {  // AAC-ELD => additional 6 bits
        objectType = 32 + br.getBits(6);
    }

    uint32_t freqIndex = br.getBits(4);

    int32_t sampleRate = 0;
    int32_t numChannels = 0;
    if (freqIndex == 15) {
        if (csd_size < 5) {
            return ERROR_MALFORMED;
        }
        sampleRate = br.getBits(24);
        numChannels = br.getBits(4);
    } else {
        static uint32_t kSamplingRate[] = {
            96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050,
            16000, 12000, 11025, 8000, 7350
        };

        if (freqIndex == 13 || freqIndex == 14) {
            return ERROR_MALFORMED;
        }

        sampleRate = kSamplingRate[freqIndex];
        numChannels = br.getBits(4);
    }

    if (numChannels == 0) {
        return ERROR_UNSUPPORTED;
    }

    int32_t prevSampleRate;
    CHECK(mLastTrack->meta->findInt32(kKeySampleRate, &prevSampleRate));

    if (prevSampleRate != sampleRate) {
        ALOGV("mpeg4 audio sample rate different from previous setting. "
             "was: %d, now: %d", prevSampleRate, sampleRate);
    }

    mLastTrack->meta->setInt32(kKeySampleRate, sampleRate);

    int32_t prevChannelCount;
    CHECK(mLastTrack->meta->findInt32(kKeyChannelCount, &prevChannelCount));

    if (prevChannelCount != numChannels) {
        ALOGV("mpeg4 audio channel count different from previous setting. "
             "was: %d, now: %d", prevChannelCount, numChannels);
    }

    mLastTrack->meta->setInt32(kKeyChannelCount, numChannels);

    return OK;
}

////////////////////////////////////////////////////////////////////////////////

MoofSource::MoofSource(
    size_t trackIndex,
    const sp<MetaData> &format,
    const sp<DataSource> &dataSource,
    uint32_t timeScale,
    sp<MoofSampleTable> &moofSampleTable)
    : mFormat(format),
      mDataSource(dataSource),
      mTimescale(timeScale),
      mCurrentSampleIndex(1),
      mLastTimeReadUs(0),
      mTrackIndex(trackIndex),
      mSampleTable(moofSampleTable),
      mIsAVC(false),
      mNALLengthSize(0),
      mStarted(false),
      mGroup(NULL),
      mBuffer(NULL),
      mWantsNALFragments(false),
      mSrcBuffer(NULL) {

    const char *mime;
    bool success = mFormat->findCString(kKeyMIMEType, &mime);
    CHECK(success);

    mIsAVC = !strcasecmp(mime, MEDIA_MIMETYPE_VIDEO_AVC);
    ALOGV("MoofSource isAVC %d %s", mIsAVC, mime);
    mFormat->dumpToLog();

    if (mIsAVC) {
        uint32_t type;
        const void *data;
        size_t size;
        CHECK(format->findData(kKeyAVCC, &type, &data, &size));

        const uint8_t *ptr = (const uint8_t *)data;

        CHECK(size >= 7);
        CHECK_EQ((unsigned)ptr[0], 1u);  // configurationVersion == 1

        // The number of bytes used to encode the length of a NAL unit.
        mNALLengthSize = 1 + (ptr[4] & 3);
    }
}

MoofSource::~MoofSource() {
    if (mStarted) {
        stop();
    }
}

status_t MoofSource::start(MetaData *params) {
    Mutex::Autolock autoLock(mLock);

    CHECK(!mStarted);

    int32_t val;
    if (params && params->findInt32(kKeyWantsNALFragments, &val)
        && val != 0) {
        mWantsNALFragments = true;
    } else {
        mWantsNALFragments = false;
    }
    CHECK(!mWantsNALFragments);

    mGroup = new MediaBufferGroup;

    int32_t max_size;
    CHECK(mFormat->findInt32(kKeyMaxInputSize, &max_size));

    mGroup->add_buffer(new MediaBuffer(max_size));

    ALOGV("MoofSource::start max_size %d wantsNALFragments=%d", max_size, mWantsNALFragments);

    mSrcBuffer = new uint8_t[max_size];
    mSrcBufferSize = max_size;

    mStarted = true;

    return OK;
}

status_t MoofSource::stop() {
    Mutex::Autolock autoLock(mLock);

    CHECK(mStarted);

    ALOGV("MoofSource::stop");

    if (mBuffer != NULL) {
        mBuffer->release();
        mBuffer = NULL;
    }

    delete[] mSrcBuffer;
    mSrcBuffer = NULL;

    delete mGroup;
    mGroup = NULL;

    mStarted = false;
    mCurrentSampleIndex = 1;

    return OK;
}

sp<MetaData> MoofSource::getFormat() {
    Mutex::Autolock autoLock(mLock);

    ALOGV("MoofSource::getFormat");
    return mFormat;
}

status_t MoofSource::read(
        MediaBuffer **out, const ReadOptions *options)
{
    Mutex::Autolock autoLock(mLock);

    CHECK(mStarted);

    uint64_t cur = mTimeSource.getRealTimeUs();
    ALOGD("MoofSource::read track %d sample %d cur %lld (us) delta %lld (us)",
          mTrackIndex, mCurrentSampleIndex, cur, cur - mLastTimeReadUs);
    mLastTimeReadUs = cur;

    *out = NULL;

    int64_t targetSampleTimeUs = -1;

    int64_t seekTimeUs = 0;
    ReadOptions::SeekMode mode;
    if (options && options->getSeekTo(&seekTimeUs, &mode)) {
        uint32_t findFlags = 0;
        ALOGV("MoofSource::read seek mode %d seekTimeUs %lld", mode, seekTimeUs);
        // fall through
    }

    off64_t offset = 0;
    size_t size;
    uint32_t cts = 0;
    bool isSyncSample = false;
    int err;
    CHECK(mBuffer == NULL);


    bool found = false;
    err = -1;
    MoofSample sample;
    err = mSampleTable->getSample(mTrackIndex, mCurrentSampleIndex, &sample);
    offset = sample.sampleOffset;
    size = sample.sampleSize;
    cts = sample.sampleTime;
    targetSampleTimeUs = sample.sampleTime;
    isSyncSample = ((sample.sampleFlags & 0x00010000) == 0);
    found = true;

    mBuffer = sample.takeBuffer();

    if (0) {
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);
        ALOGV("track %d sample %d delay %.3f offset %lld size %d time %lld err %d",
              mTrackIndex, mCurrentSampleIndex, (tv.tv_sec + tv.tv_usec / 1.0e6) - sample.mfts,
              offset, size, targetSampleTimeUs, err);
    }
    if (0) {
        ALOGD("buffer size %d range_offset %d range_length %d",
              mBuffer->size(), mBuffer->range_offset(), mBuffer->range_length());
        if (mBuffer->range_offset() > 0
            || (mBuffer->range_offset() + mBuffer->range_length()) > mBuffer->size())
            ALOGE("bad buffer?");
        char prefix[128];
        snprintf(prefix, sizeof(prefix), "this %p track %d sample %d: ",
                 this, mTrackIndex, mCurrentSampleIndex);
        hexdump(prefix, mBuffer->data(), mBuffer->range_length());
    }

    if (err != OK) {
        return err;
    }

    ++mCurrentSampleIndex;
    *out = mBuffer;
    mBuffer = NULL;

    return OK;
}

MoofExtractor::Track *MoofExtractor::findTrackByMimePrefix(
        const char *mimePrefix) {
    for (Track *track = mFirstTrack; track != NULL; track = track->next) {
        const char *mime;
        if (track->meta != NULL
                && track->meta->findCString(kKeyMIMEType, &mime)
                && !strncasecmp(mime, mimePrefix, strlen(mimePrefix))) {
            return track;
        }
    }

    return NULL;
}

static bool LegacySniffMoof(
        const sp<DataSource> &source, String8 *mimeType, float *confidence) {
    uint8_t header[8];

    ssize_t n = source->readAt(4, header, sizeof(header));
    if (n < (ssize_t)sizeof(header)) {
        return false;
    }

    if (!memcmp(header, "ftypdash", 8)) {
        *mimeType = MEDIA_MIMETYPE_CONTAINER_DASH;
        *confidence = 0.4;

        return true;
    }

    return false;
}

static bool isCompatibleBrand(uint32_t fourcc) {
    static const uint32_t kCompatibleBrands[] = {
        FOURCC('d', 'a', 's', 'h'),
    };

    for (size_t i = 0;
         i < sizeof(kCompatibleBrands) / sizeof(kCompatibleBrands[0]);
         ++i) {
        if (kCompatibleBrands[i] == fourcc) {
            return true;
        }
    }

    return false;
}

// Attempt to actually parse the 'ftyp' atom and determine if a suitable
// compatible brand is present.
// Also try to identify where this file's metadata ends
// (end of the 'moov' atom) and report it to the caller as part of
// the metadata.
static bool BetterSniffMoof(
        const sp<DataSource> &source, String8 *mimeType, float *confidence,
        sp<AMessage> *meta) {
    // We scan up to 128 bytes to identify this file as an MP4.
    static const off64_t kMaxScanOffset = 128ll;

    off64_t offset = 0ll;
    bool foundGoodFileType = false;
    off64_t moovAtomEndOffset = -1ll;
    bool done = false;

    while (!done && offset < kMaxScanOffset) {
        uint32_t hdr[2];
        if (source->readAt(offset, hdr, 8) < 8) {
            return false;
        }

        uint64_t chunkSize = ntohl(hdr[0]);
        uint32_t chunkType = ntohl(hdr[1]);
        off64_t chunkDataOffset = offset + 8;

        if (chunkSize == 1) {
            if (source->readAt(offset + 8, &chunkSize, 8) < 8) {
                return false;
            }

            chunkSize = ntoh64(chunkSize);
            chunkDataOffset += 8;

            if (chunkSize < 16) {
                // The smallest valid chunk is 16 bytes long in this case.
                return false;
            }
        } else if (chunkSize < 8) {
            // The smallest valid chunk is 8 bytes long.
            return false;
        }

        off64_t chunkDataSize = offset + chunkSize - chunkDataOffset;

        switch (chunkType) {
            case FOURCC('f', 't', 'y', 'p'):
            {
                if (chunkDataSize < 8) {
                    return false;
                }

                uint32_t numCompatibleBrands = (chunkDataSize - 8) / 4;
                for (size_t i = 0; i < numCompatibleBrands + 2; ++i) {
                    if (i == 1) {
                        // Skip this index, it refers to the minorVersion,
                        // not a brand.
                        continue;
                    }

                    uint32_t brand;
                    if (source->readAt(
                                chunkDataOffset + 4 * i, &brand, 4) < 4) {
                        return false;
                    }

                    brand = ntohl(brand);

                    if (isCompatibleBrand(brand)) {
                        foundGoodFileType = true;
                        break;
                    }
                }

                if (!foundGoodFileType) {
                    return false;
                }

                break;
            }

            case FOURCC('m', 'o', 'o', 'v'):
            {
                moovAtomEndOffset = offset + chunkSize;

                done = true;
                break;
            }

            default:
                break;
        }

        offset += chunkSize;
    }

    if (!foundGoodFileType) {
        return false;
    }

    *mimeType = MEDIA_MIMETYPE_CONTAINER_DASH;
    *confidence = 0.8f;

    if (moovAtomEndOffset >= 0) {
        *meta = new AMessage;
        (*meta)->setInt64("meta-data-size", moovAtomEndOffset);

        ALOGV("found metadata size: %lld mimetype dash",
              moovAtomEndOffset);
    }

    return true;
}

bool SniffMoof(
        const sp<DataSource> &source, String8 *mimeType, float *confidence,
        sp<AMessage> *meta) {
    if (BetterSniffMoof(source, mimeType, confidence, meta)) {
        return true;
    }

    if (LegacySniffMoof(source, mimeType, confidence)) {
        ALOGW("Identified supported mpeg4 through LegacySniffMoof.");
        return true;
    }

    return false;
}

}  // namespace android
