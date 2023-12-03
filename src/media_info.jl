using FFMPEG
using Glob

function get_fps(file::AbstractString, streamno::Integer = 0)
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
    fps_strs = FFMPEG.exe(
        `-v 0 -of compact=p=0 -select_streams v:0 -show_entries stream=r_frame_rate $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug fps_strs
	try
		fps = split(fps_strs[1], '=')[2]
		if occursin("No such file or directory", fps)
			error("Could not find file $file")
		elseif occursin("N/A", fps)
			return nothing
		end
		return reduce(//, parse.(Int, split(fps,'/')) )
		# return round(reduce(/, parse.(Float64, split(fps,'/')) ), digits=3)
	catch err
		@debug(err)
		return NaN
	end
end

function get_framerate(file::AbstractString, streamno::Integer = 0, video_or_audio="v")
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
	if video_or_audio == "v"
		entries = "r_frame_rate"
	elseif video_or_audio == "a"
		entries = "sample_rate"
	end
	
    fps_strs = FFMPEG.exe(
        `-v 0 -of compact=p=0 -select_streams $video_or_audio:$streamno -show_entries stream="$entries" $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug fps_strs
	try
		fps = split(fps_strs[1], '=')[2]
		if occursin("No such file or directory", fps)
			error("Could not find file $file")
		elseif occursin("N/A", fps)
			return nothing
		end
		return reduce(//, parse.(Int, split(fps,'/')) )
		# return round(reduce(/, parse.(Float64, split(fps,'/')) ), digits=3)
	catch err
		@debug(err)
		return NaN
	end
end

function get_whatever(file::AbstractString, streamno::Integer = 0, video_or_audio="v"; entries_custom=nothing)
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
	if video_or_audio == "v"
		entries = "r_frame_rate"
	elseif video_or_audio == "a"
		entries = "sample_rate"
	end
	if !isnothing(entries)
		entries = entries_custom
	end
	
    fps_strs = FFMPEG.exe(
        `-v 0 -of compact=p=0 -select_streams $video_or_audio:$streamno -show_entries stream="$entries" $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug fps_strs
	try
		fps = split(fps_strs[1], '=')[2]
		if occursin("No such file or directory", fps)
			error("Could not find file $file")
		elseif occursin("N/A", fps)
			return nothing
		end
		if occursin('/', fps)
			return reduce(//, parse.(Int, split(fps,'/')) )
		else
			return parse.(Float64, fps) 
		end
		# return round(reduce(/, parse.(Float64, split(fps,'/')) ), digits=3)
	catch err
		@debug(err)
		return NaN
	end
end

"""
    get_number_frames(file [, streamno])
Query the the container `file` for the number of frames in video stream
`streamno` if applicable, instead returning `nothing` if the container does not
report the number of frames. Will not decode the video to count the number of
frames in a video.
"""
function get_number_frames(file::AbstractString, streamno::Integer = 0)
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
    frame_strs = FFMPEG.exe(
		`-v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets $file`, #-hide_banner
        # `-v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets -of csv=p=0 $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug frame_strs
	frame_str = frame_strs[1]
	# num_frames = parse(Int, split(frame_str,'=')[end])
    if occursin("No such file or directory", frame_str)
        error("Could not find file $file")
    elseif occursin("N/A", frame_str)
        return NaN
    end
	
	try
		frame_str = frame_strs[2]
	    return parse(Int, split(frame_str,'=')[end])
	catch err
		@debug (err)
		return NaN
	end
end

function get_duration(file::AbstractString, streamno::Integer = 0)
	try
		streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))

		
		frame_strs = FFMPEG.exe(
			`-show_entries format=duration -v quiet -of csv="p=0" $file`,
			# `-v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets -of csv=p=0 $file`,
			command = FFMPEG.ffprobe,
			collect = true,
		)
		@debug frame_strs
		frame_str = frame_strs[1]
		# num_frames = parse(Int, split(frame_str,'=')[end])
		if occursin("No such file or directory", frame_str)
			error("Could not find file $file")
		elseif occursin("N/A", frame_str)
			@debug "manually calculate duration"
			return get_number_frames(file) / get_fps(file) |> Float64
		end
	
	# try
		frame_str = frame_str
	    return parse(Float64, split(frame_str,'=')[end])
	catch err
		@debug (err)
		return NaN
	end
end

function get_videos_audiodata(vidfname)
	strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -f s16le -acodec pcm_s16le -loglevel error -`)# .|> Int16;
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# vid_audioFS =  get_framerate(vidfname, 0, "a")
	# length(strs)/vid_audioFS/2 #vid_auDur
	# get_duration(vidfname) #vid_Dur

	# # method 1, fast and creates array
	# vid_audiodata = Array{Int16}(undef, Int(length(strs)/2)) # 8bits to 16bits per frame
	# map!( x -> strs[2x[1]-1] + Int16(256)*strs[2x[1]] , vid_audiodata, 1:length(vid_audiodata)) # convert to 8bits little endian to Int16 merging each 2 bytes to 1 frame
	# vid_audiodata, get_framerate(vidfname, 0, "a")
	
	# # method 2, slow but easy to read, create array
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# method 3, fastest, doesnt create array, does the job
	vid_audiodata = reinterpret(Int16, strs)#, get_framerate(vidfname, 0, "a")
	if get_whatever(vidfname, 0, "a"; entries_custom="channels")>1
		vid_audiodata = reshape(vid_audiodata, get_whatever(vidfname, 0, "a"; entries_custom="channels")|>Int, :)'
	end
	vid_audiodata, get_framerate(vidfname, 0, "a")
end

function get_videos_audiodata_direct(vidfname)
	strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -f s16le -acodec pcm_s16le -loglevel error -`)# .|> Int16;
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# vid_audioFS =  get_framerate(vidfname, 0, "a")
	# length(strs)/vid_audioFS/2 #vid_auDur
	# get_duration(vidfname) #vid_Dur

	# # method 1, fast and creates array
	# vid_audiodata = Array{Int16}(undef, Int(length(strs)/2)) # 8bits to 16bits per frame
	# map!( x -> strs[2x[1]-1] + Int16(256)*strs[2x[1]] , vid_audiodata, 1:length(vid_audiodata)) # convert to 8bits little endian to Int16 merging each 2 bytes to 1 frame
	# vid_audiodata, get_framerate(vidfname, 0, "a")
	
	# # method 2, slow but easy to read, create array
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# method 3, fastest, doesnt create array, does the job
	reinterpret(Int16, strs), get_framerate(vidfname, 0, "a")
end

## FIXME: output codec to automatically convert to raw form
function get_videos_audiodata_all(vidfname, streamno=:all)
	if streamno == :all
		info = get_media_info(vidfname)
		streamno = 0:length(info["streams"])-1
		d = get_videos_audiodata_all.(Ref(vidfname), streamno)
		return hcat(map(x->x[1],d)...), d[1][2]
	end
	strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -map 0:$streamno -f s16le -acodec pcm_s16le -loglevel error -`)# .|> Int16;
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# vid_audioFS =  get_framerate(vidfname, 0, "a")
	# length(strs)/vid_audioFS/2 #vid_auDur
	# get_duration(vidfname) #vid_Dur

	# # method 1, fast and creates array
	# vid_audiodata = Array{Int16}(undef, Int(length(strs)/2)) # 8bits to 16bits per frame
	# map!( x -> strs[2x[1]-1] + Int16(256)*strs[2x[1]] , vid_audiodata, 1:length(vid_audiodata)) # convert to 8bits little endian to Int16 merging each 2 bytes to 1 frame
	# vid_audiodata, get_framerate(vidfname, 0, "a")
	
	# # method 2, slow but easy to read, create array
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# method 3, fastest, doesnt create array, does the job
	vid_audiodata = reinterpret(Int16, strs)#, get_framerate(vidfname, 0, "a")
	if get_whatever(vidfname, streamno, "a"; entries_custom="channels")>1
		vid_audiodata = reshape(vid_audiodata, get_whatever(vidfname, 0, "a"; entries_custom="channels")|>Int, :)'
	end
	vid_audiodata, get_framerate(vidfname, 0, "a")
end

function get_media_info(fname)
	@ffmpeg_env read(`$ffprobe $fname -loglevel error -v quiet -print_format json -show_format -show_streams`, String) |> JSON.parse
end

function get_duration_smart(fname; vidtype=r".mkv|.MP4|.avi|.mp4", autype=r".wav|.mat|.flac|.mp3|.aac")
	occursin(vidtype, fname) && return get_number_frames(fname) / get_fps(fname) |> Float64
	occursin(autype, fname)  && return get_duration(fname)
end

@info "LOADED!\tmedia_info"