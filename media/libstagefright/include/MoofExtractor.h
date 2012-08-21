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

#ifndef MOOF_EXTRACTOR_H_

#define MOOF_EXTRACTOR_H_

#include <media/stagefright/MediaExtractor.h>
#include <utils/List.h>
#include <utils/Vector.h>
#include <utils/String8.h>

namespace android {

struct AMessage;
class DataSource;
 class MoofDataSource;
class MoofSampleTable;
class String8;
struct Track;

class MoofExtractor : public MediaExtractor {
public:
    // Extractor assumes ownership of "source".
    MoofExtractor(const sp<DataSource> &source);

    virtual size_t countTracks();
    virtual sp<MediaSource> getTrack(size_t index);
    virtual sp<MetaData> getTrackMetaData(size_t index, uint32_t flags);

    virtual sp<MetaData> getMetaData();
    virtual uint32_t flags() const;

    // for DRM
    virtual char* getDrmTrackInfo(size_t trackID, int *len);

protected:
    virtual ~MoofExtractor();

private:
    struct Track {
        Track *next;
        sp<MetaData> meta;
        uint32_t timescale;
        bool includes_expensive_metadata;
        bool skipTrack;
    };

    sp<MoofDataSource> mDataSource;
    sp<MoofSampleTable> mMoofSampleTable;
    status_t mInitCheck;
    bool mHasVideo;

    int mNumTracks;
    Track *mFirstTrack, *mLastTrack;
    uint64_t mOffset;

    sp<MetaData> mFileMetaData;

    Vector<uint32_t> mPath;
    String8 mLastCommentMean;
    String8 mLastCommentName;
    String8 mLastCommentData;

    status_t readMetaData();
    status_t parseChunk(off64_t *offset, int depth);
    status_t parseMetaData(off64_t offset, size_t size);

    status_t updateAudioTrackInfoFromESDS_MPEG4Audio(
            const void *esds_data, size_t esds_size);

    static status_t verifyTrack(Track *track);

    struct SINF {
        SINF *next;
        uint16_t trackID;
        uint8_t IPMPDescriptorID;
        ssize_t len;
        char *IPMPData;
    };

    SINF *mFirstSINF;

    bool mIsDrm;
    static int32_t readSize(off64_t offset,
                            const sp<DataSource> DataSource, uint8_t *numOfBytes);

    status_t parseDrmSINF(off64_t *offset, off64_t data_offset);

    status_t parseTrackHeader(off64_t data_offset, off64_t data_size);

    Track *findTrackByMimePrefix(const char *mimePrefix);

    MoofExtractor(const MoofExtractor &);
    MoofExtractor &operator=(const MoofExtractor &);
};

bool SniffMoof(
        const sp<DataSource> &source, String8 *mimeType, float *confidence,
        sp<AMessage> *);

}  // namespace android

#endif  // MOOF_EXTRACTOR_H_
