using VideoIO

function readImages(videoFile::String, eventTimings)

    """

        readImages(videoFile, eventTimings)

    Function to extract frames from particular times given in an array as seconds. 

    Inputs are the video file location and array containing time instances of event in seconds and as a float.
    Output is an array of all the images of type ::Array{RGB{Normed{UInt8,8}},nx,ny,nImages}    

    """


    nImages = length(eventTimings)
    f = VideoIO.openvideo(videoFile)
    nx, ny = raw_frame_size(f)
    # seekstart(f)
    # img = read(f)
    # nx,ny=size(img)
    @show nx, ny, nImages
    
    images = zeros(out_frame_eltype(f),nx,ny,nImages) # #images = zeros(RGB{Normed{UInt8,8}},nx,ny,nImages)
    # images = Array{typeof(img),1}(undef, nImages)
    for i = 1:nImages
        seek(f, eventTimings[i])
        img = read(f)
        images[:,:,i].=img[:,:]
        # images[i]=read(f)
    end
    close(f)
    images
end

function readImage(vidFname::String, event_time)
    f = VideoIO.openvideo(vidFname)
    seek(f, event_time)
    img = read(f)
    close(f)
    return img
end

function readImage(vid::VideoIO.VideoReader, event_time)
    # f = VideoIO.openvideo(vidFname)
    try
        seek(vid, event_time)
        img = read(vid)
        return img
    catch err
        @error "Failed at reading image frame from VIDEO"
        if eof(vid)
            @error "Seeked to eof, unable to retrieve image frame!!"
        end
        @error (err)
        @error event_time
        return []
    end
end
