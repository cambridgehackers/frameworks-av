/*
 * Copyright (C) 2010 The Android Open Source Project
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

#include "SineSource.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/socket.h>
#include <netdb.h>
#include <errno.h>

#include <binder/ProcessState.h>
#include <camera/Camera.h>
#include <gui/Surface.h>
#include <media/stagefright/foundation/ADebug.h>
#include <media/stagefright/AudioPlayer.h>
#include <media/stagefright/CameraSource.h>
#include <media/stagefright/MediaBufferGroup.h>
#include <media/stagefright/MediaDefs.h>
#include <media/stagefright/MetaData.h>
#include <media/stagefright/MPEG4Writer.h>
#include <media/stagefright/MoofWriter.h>
#include <media/stagefright/OMXClient.h>
#include <media/stagefright/OMXCodec.h>
#include <media/MediaPlayerInterface.h>
#include <media/mediarecorder.h>
#include <gui/SurfaceComposerClient.h>

using namespace android;

// Print usage showing how to use this utility to record videos
static void usage(const char *me) {
    fprintf(stderr, "usage: %s\n", me);
    fprintf(stderr, "       -h(elp)\n");
    fprintf(stderr, "       -o outputfile  (default: /sdcard/output.mp4)\n");
    fprintf(stderr, "       -a hostname to which to post video (default: unspecified)\n");
    fprintf(stderr, "       -p port to which to post video (default: 80)\n");
    fprintf(stderr, "       -b bit rate in bits per second (default: 300000)\n");
    fprintf(stderr, "       -f frame rate in frames per second (default: 30)\n");
    fprintf(stderr, "       -i I frame interval in seconds (default: 1)\n");
    fprintf(stderr, "       -n number of frames to be recorded (default: 300)\n");
    fprintf(stderr, "       -w width in pixels (default: 720)\n");
    fprintf(stderr, "       -t height in pixels (default: 480)\n");
    fprintf(stderr, "       -v video codec: [0] H264 [1] MPEG_4_SP [2] H263 (default: 0)\n");
    fprintf(stderr, "       -c container format: [0] MPEG4 [1] DASH (default: 0)\n");
    fprintf(stderr, "The output file is /sdcard/output.mp4\n");
    exit(1);
}

volatile bool sDone = false;

class MyMediaRecorderListener : public MediaRecorderListener
{
public:
    virtual void notify(int msg, int ext1, int ext2) {
        fprintf(stderr, "MyMediaRecorderListener msg=%x ext1=%x ext2=%x", msg, ext1, ext2);
        if (msg == MEDIA_RECORDER_TRACK_EVENT_INFO) {
            if (ext1 & MEDIA_RECORDER_TRACK_INFO_COMPLETION_STATUS)
                sDone = true;
        }
    }
};

int connectToHost(const char *hostname, int port)
{
    struct addrinfo *res;

    char servname[64];
    snprintf(servname, sizeof(servname), "%d", port);

    fprintf(stderr, "connectToHost hostname=%s servname=%s\n", hostname, servname);

    int error = getaddrinfo(hostname, servname, NULL, &res);
    if (error)
	fprintf(stderr, "getaddrinfo error=%d %s\n", error, gai_strerror(error));

    //int sock = socket(res[0].ai_family, res[0].ai_socktype, res[0].ai_protocol);
    int sock = socket(PF_INET, SOCK_STREAM, 0);
    if (sock < 0)
	fprintf(stderr, "sock=%d errno=%d\n", sock, errno);

    struct sockaddr address;
    error = connect(sock, res[0].ai_addr, res[0].ai_addrlen);
    if (error < 0) {
	fprintf(stderr, "failed to connect: errno=%d\n", errno);
	return -1;
    }
    fprintf(stderr, "connected with socket %d\n", sock);
    return sock;
}

void httpPost(int sock, const char *hostname, const char *path)
{
    char line[1024];
    int nchars = snprintf(line, sizeof(line), "POST %s HTTP/1.1\r\n", path);
    nchars += snprintf(line+nchars, sizeof(line)-nchars, "Host: %s\r\n", hostname);
    nchars += snprintf(line+nchars, sizeof(line)-nchars, "User-Agent: nrcc-webcam\r\n");
    nchars += snprintf(line+nchars, sizeof(line)-nchars, "Content-Type: video/avc\r\n");
    nchars += snprintf(line+nchars, sizeof(line)-nchars, "Content-Length: 2000000000\r\n");
    nchars += snprintf(line+nchars, sizeof(line)-nchars, "\r\n");
    int nWritten = write(sock, line, nchars);
    if (nWritten < nchars) {
	fprintf(stderr, "short write %d < %d errno=%d", nWritten, nchars, errno);
    } else {
        fprintf(stderr, "posted %d characters\n", nWritten);
    }
}

class CameraRecorder : public virtual RefBase {
public:
    CameraRecorder(int cameraNumber)
        : mCameraNumber(cameraNumber) {
    }
    
    void connect() {
        status_t err = OK;
        mCamera = Camera::connect(0);
        mICamera = mCamera->remote();
    }

    struct Params {
        int width;
        int height;
        int codec;
        int output_format;
        int frameRateFps;
        int nFrames;        
    };

    void startPreview(const Params &params) {
        sp<SurfaceComposerClient> surfaceComposerClient = new SurfaceComposerClient();
        sp<SurfaceControl> surfaceControl =
            surfaceComposerClient->createSurface(0, params.width, params.height, PIXEL_FORMAT_RGB_888);

        SurfaceComposerClient::openGlobalTransaction();
        surfaceControl->setLayer(0x40000000);
        SurfaceComposerClient::closeGlobalTransaction();

        mSurface = surfaceControl->getSurface();

        mCamera->setPreviewDisplay(mSurface);
        mCamera->startPreview();
    }

    void startRecording(const Params &params, int fd) {
        Size videoSize(params.width, params.height);
        sp<ICameraRecordingProxy> cameraRecordingProxy = mCamera->getRecordingProxy();
        mCamera->unlock();

        int encoder = VIDEO_ENCODER_H264;
        switch (params.codec) {
        case 2:
            encoder = VIDEO_ENCODER_H263;
            break;
        case 1:
            encoder = VIDEO_ENCODER_MPEG_4_SP;
            break;
        case 0:
        default:
            encoder = VIDEO_ENCODER_H264;
            break;
        }

        sp<MediaRecorder> recorder = new MediaRecorder();
        mRecorder = recorder;
        recorder->setCamera(mICamera, cameraRecordingProxy);
        recorder->setVideoSource(VIDEO_SOURCE_CAMERA);
        recorder->setOutputFormat(params.output_format);

        recorder->setVideoSize(params.width, params.height);
        recorder->setVideoFrameRate(params.frameRateFps);
        recorder->setVideoEncoder(encoder);
        char rparams[64];
        snprintf(rparams, sizeof(rparams), "max-duration=%d", params.nFrames/params.frameRateFps*1000);
        recorder->setParameters(String8(rparams));

        recorder->setOutputFile(fd, 0, 0);

        recorder->setPreviewSurface(mSurface);

        recorder->setListener(new MyMediaRecorderListener());

        recorder->prepare();

        recorder->start();

    }
    
    void stopRecording() {
        close(mFd);
        mRecorder->close();
        mRecorder->release();
        mCamera->disconnect();
    }

private:
    int mCameraNumber;
    int mFd;
    sp<Camera> mCamera;
    sp<ICamera> mICamera;
    sp<Surface> mSurface;
    sp<MediaRecorder> mRecorder;
};

int main(int argc, char **argv) {

    // Default values for the program if not overwritten
    int frameRateFps = 30;
    int width = 720;
    int height = 480;
    int bitRateBps = 300000;
    int iFramesIntervalSeconds = 1;
    int nFrames = 300;
    int codec = 0;
    enum {
        CONTAINER_MPEG4 = 0,
        CONTAINER_DASH = 1
    };
    output_format output_format = OUTPUT_FORMAT_THREE_GPP;
    const char *fileName = "/sdcard/output.mp4";
    const char *hostname = 0;
    const char *streamNames[4];
    sp<CameraRecorder> recorders[4];
    int numStreams = 0;
    int port = 80;

    android::ProcessState::self()->startThreadPool();
    int res;
    while ((res = getopt(argc, argv, "a:b:c:f:i:n:w:t:l:p:s:v:h")) >= 0) {
        switch (res) {
            case 'b':
            {
                bitRateBps = atoi(optarg);
                break;
            }

            case 'f':
            {
                frameRateFps = atoi(optarg);
                break;
            }

            case 'a':
            {
                hostname = optarg;
                break;
            }

            case 'i':
            {
                iFramesIntervalSeconds = atoi(optarg);
                break;
            }

            case 'n':
            {
                nFrames = atoi(optarg);
                break;
            }

            case 'o':
            {
                fileName = optarg;
                break;
            }

            case 'p':
            {
                port = atoi(optarg);
                break;
            }

            case 's':
            {
                streamNames[numStreams++] = optarg;
                break;
            }

            case 'w':
            {
                width = atoi(optarg);
                break;
            }

            case 't':
            {
                height = atoi(optarg);
                break;
            }

            case 'v':
            {
                codec = atoi(optarg);
                if (codec < 0 || codec > 2) {
                    usage(argv[0]);
                }
                break;
            }

            case 'c':
            {
                int param = atoi(optarg);
                if (param < 0 || param > 2) {
                    usage(argv[0]);
                }
                if (CONTAINER_DASH == param)
                    output_format = OUTPUT_FORMAT_MPEG_4;
                break;
            }

            case 'h':
            default:
            {
                usage(argv[0]);
                break;
            }
        }
    }

    for (int cam = 0; cam < numStreams; cam++) {
        sp<CameraRecorder> recorder = new CameraRecorder(cam);
        recorders[cam] = recorder;

        recorder->connect();
        CameraRecorder::Params params;
        params.width = width;
        params.height = height;
        params.codec = codec;
        params.output_format = output_format;
        params.frameRateFps = frameRateFps;
        params.nFrames = nFrames;

        recorder->startPreview(params);
        
        int fd = -1;
        if (hostname == 0) {
            fd = open(fileName, O_WRONLY|O_CREAT, 0660);
        } else {
            fd = connectToHost(hostname, port);
            char path[256];
            snprintf(path, sizeof(path), "/%s", streamNames[cam]);
            httpPost(fd, hostname, path);
        }
        recorder->startRecording(params, fd);
    }

    int64_t start = systemTime();
    while (!sDone) {
        // wait for video capture to complete
    }
    int64_t end = systemTime();

    for (int cam = 0; cam < numStreams; cam++) {
        recorders[cam]->stopRecording();
    }

    fprintf(stderr, "$\n");

    fprintf(stderr, "encoding %d frames in %lld us\n", nFrames, (end-start)/1000);
    fprintf(stderr, "encoding speed is: %.2f fps\n", (nFrames * 1E9) / (end-start));
    return 0;
}


