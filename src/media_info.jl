using FFMPEG
using JSON
# using Glob

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
		return missing
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
		return missing
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
		return missing
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
        return missing
    end
	
	try
		frame_str = frame_strs[2]
	    return parse(Int, split(frame_str,'=')[end])
	catch err
		@debug (err)
		return missing
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
		return missing
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
	try 
		@ffmpeg_env read(`$ffprobe $fname -loglevel error -v quiet -print_format json -show_format -show_streams`, String) |> JSON.parse
	catch e
		@warn "Error getting media info for $fname, with Julia's FFMPEG package"
		@info "try with system default ffprobe instead"
		read(`$ffprobe $fname -loglevel error -v quiet -print_format json -show_format -show_streams`, String) |> JSON.parse
	end
end

function get_duration_smart(fname; vidtype=r".mkv|.avi|.mp4", autype=r".wav|.mat|.flac|.mp3|.aac")
	fname_lowered = lowercase(fname)
	occursin(vidtype, fname_lowered) && return get_number_frames(fname) / get_fps(fname) |> Float64
	occursin(autype, fname_lowered)  && return get_duration(fname)
end

function mediatype(filename; vidtypes=vidtypes, autypes=autypes)
	filetype = splitext(filename)[2]
	if ((filetype .== vidtypes) |> sum ) == 1
		return "video"
	elseif ((filetype .== autypes) |> sum ) == 1
		return "audio"
	else
		return "other"
	end
end

"""
  get_ffmpeg_metadata(filename::AbstractString)

Obtain metadata from a media file using ffprobe (part of ffmpeg). Returns a parsed JSON object with all available metadata.
"""
function get_ffmpeg_metadata(fname::AbstractString)
  output = @ffmpeg_env read(`ffprobe -v quiet -print_format json -show_format -show_streams $fname`, String)
  return JSON.parse(output)
end

"""
concat_media(filelist::Vector{String}, outputfile_dir::String; flag_output_auto=true, flag_overwrite=false)

Concatenate multiple media files into a single output file using ffmpeg's concat demuxer,
while attempting to preserve metadata from the first input file.

Arguments
- filelist::Vector{String}: Non-empty vector of input file paths to concatenate. All files should
	share a compatible container/codec for stream-copy concatenation (the concat demuxer requirement).
- outputfile_dir::String: If `flag_output_auto` is true this is treated as the directory where the
	automatically named output file will be written. If `flag_output_auto` is false this is treated
	as the explicit path for the output file.

Keyword arguments
- flag_output_auto::Bool = true: When true, the output filename is auto-generated as
	joinpath(outputfile_dir, "combined__\$(join(basename.(filelist), '_'))<ext>") where <ext> is the
	extension of the first input file. When false, `outputfile_dir` is used as the output file path.
- flag_overwrite::Bool = false: When true, ffmpeg is invoked with overwrite enabled (passes -y)
	so existing output files are replaced.

Behavior / Implementation notes
- A temporary file named "temp_filelist.txt" is written beside the intended output location. Each line
	is formatted as:  file 'path/to/input'
	This file is used as input to ffmpeg's concat demuxer (-f concat -safe 0 -i <tempfile>).
- The function also adds the first input file as a second ffmpeg input and uses `-map_metadata 1`
	to copy global metadata from that first file into the concatenated output while `-map 0` copies
	the concatenated streams.
- The primary ffmpeg invocation uses `@ffmpeg_env run(...)` (FFMPEG.jl environment). If that run
	throws an error, the function retries with a system call and forces overwrite (-y).
- The temporary file removal is present in the code but commented out; the temporary list file may
	remain on disk after execution.

Returns
- outputfile::String: The full path to the concatenated output file (the generated or provided path).

Errors and side effects
- Throws an exception if ffmpeg fails on both the primary and fallback attempts.
- Requires a working ffmpeg installation in PATH and optionally FFMPEG.jl for `@ffmpeg_env`.
- Inputs must be compatible for stream-copy concatenation. If not, ffmpeg may fail or produce invalid output.
- Temporary file "temp_filelist.txt" may persist if cleanup is not enabled or if the process is interrupted.

Example
- concat_media(["a.mp4","b.mp4"], "/out/dir"; flag_output_auto=true, flag_overwrite=false)
	-> creates "/out/dir/combined__a.mp4_b.mp4.mp4" (extension based on first input) and returns that path.
"""
function concat_media(filelist::Vector{String}, outputfile_dir::String; flag_output_auto=true, flag_overwrite=false)
	# create a temporary text file listing the input files
	ext = splitext(filelist[1])[2]
	if flag_output_auto
		outputfile = joinpath(outputfile_dir, "combined__$(join(basename.(filelist), '_'))$ext")
	else
		outputfile = outputfile_dir
	end
	tmpfile = joinpath(dirname(outputfile), "temp_filelist.txt")
	open(tmpfile, "w") do io
		for f in filelist
			write(io, "file '$f'\n")
		end
	end

	overwrite = flag_overwrite ? " -y" : ""

	# use the concat demuxer as input 0 and add the first file as input 1,
	# then copy streams from input 0 and copy metadata from input 1
	firstfile = filelist[1]
	cmd = `ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i $tmpfile -i $firstfile -c copy -map 0 -map_metadata 1 $outputfile`
	# _cmd = `ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i $tmpfile -i $firstfile -c copy -map 0 -map_metadata 1 -map_chapters 1 -movflags use_metadata_tags $outputfile`
	
	# run the command
	try
		@ffmpeg_env run(`$cmd$overwrite`)
	catch err
		@warn "Concatenation failed with Julia's FFMPEG, trying again with system call: $err"
		run(`$cmd -y`)
	finally
		# delete the temporary file
		# rm(tmpfile; force=true)
	end
	return outputfile
end

"""
	combine_vidau(newvidname, aufname_list; vidau_syncdiff=0,
				  MERGE_VID_AU_DYNAMIC_NORM=false, rx_vect=rx_vect,
				  flag_rm_oldfile=false, kwargs...)

Combine a video file with one or more audio files, optionally normalizing audio
and applying an audio/video offset.

Behavior
- If `aufname_list` is a Vector of audio paths and contains more than one file,
  the audio files will first be concatenated (via `concat_media`) into a single
  temporary audio file which is then used for muxing.
- If `MERGE_VID_AU_DYNAMIC_NORM` is true, the function applies FFmpeg's
  loudness normalization filter (`loudnorm`) to the audio while remuxing the
  video and produces an output file appended with `_DYnormalized-audio.mp4`.
- If `MERGE_VID_AU_DYNAMIC_NORM` is false, the function:
  1. Runs FFmpeg's `astats=metadata=1` to extract per-channel peak levels.
  2. Selects the relevant channels based on `rx_vect` via
	 `get_relevant_channels(rx_vect)`.
  3. Computes a gain to normalize the loudest relevant channel to 0 dBFS and
	 applies this gain using the `volume` audio filter.
  4. Produces an output file appended with `_normalized-audio.mp4`.
- The muxing command includes `-itsoffset vidau_syncdiff` to shift the audio
  start relative to the video by `vidau_syncdiff` seconds.
- The function attempts to run FFmpeg via the `@ffmpeg_env` wrapper and falls
  back to a system call if that invocation fails.

Arguments
- newvidname::AbstractString
	Path to the input video file (and base for the output filename).
- aufname_list::Union{AbstractString, Vector{<:AbstractString}}
	A single audio path or a vector of audio paths to be concatenated/muxed.

Keyword arguments
- vidau_syncdiff::Real=0
	Time offset (in seconds) to apply to the audio stream relative to the
	video when muxing (`-itsoffset`).
- MERGE_VID_AU_DYNAMIC_NORM::Bool=false
	If true, use FFmpeg's loudnorm filter (two-pass style loudness
	normalization) instead of a simple peak-based volume adjustment.
- rx_vect
	Receiver/channel selection vector used by `get_relevant_channels` to pick
	which channels' peaks are considered when computing normalization gain.
- flag_rm_oldfile::Bool=false
	If true, remove the original `newvidname` file after successful creation of
	the normalized/muxed output.
- kwargs...
	Additional keyword arguments are accepted but not consumed by the current
	implementation (preserved for forward compatibility).

Return
- String
	Path to the created muxed/normalized MP4 file:
	either "<newvidname>_DYnormalized-audio.mp4" or
	"<newvidname>_normalized-audio.mp4" depending on `MERGE_VID_AU_DYNAMIC_NORM`.

Side effects and cleanup
- Creates a temporary file when analyzing audio with FFmpeg's `astats`, and
  removes it after parsing.
- May create a concatenated temporary audio file when multiple audio inputs are
  provided; that temporary file is removed after muxing.
- May delete the original video file if `flag_rm_oldfile` is true and the new
  output file exists.
- Prints the FFmpeg command to stdout before running it.
- Relies on external functions/variables: `ffmpeg` (binary or command wrapper),
  `@ffmpeg_env`, `concat_media`, and `get_relevant_channels`. These must be
  available in the calling scope.

Errors
- Propagates errors from FFmpeg if both the primary (`@ffmpeg_env`) and the
  fallback system invocation fail.
- If `aufname_list` is not a string or array of strings, behavior is undefined.

Example
	# Single audio file, no dynamic normalization, 0.2s audio delay
	out = combine_vidau("video.mp4", "audio.wav"; vidau_syncdiff=0.2)

	# Multiple audio files, dynamic loudness normalization, remove old video
	out = combine_vidau("video.mp4", ["a1.wav","a2.wav"];
					   MERGE_VID_AU_DYNAMIC_NORM=true,
					   flag_rm_oldfile=true)
"""
#~ combine video and audio
function combine_vidau(newvidname, aufname_list; vidau_syncdiff=0, MERGE_VID_AU_DYNAMIC_NORM=false, rx_vect=rx_vect, flag_rm_oldfile=false, kwargs...) #TODO: swap audio channel according to rx_vect location to correspond Left, Right, Center
    if aufname_list isa Array
        if length(aufname_list)==1
            aufname = aufname_list[1]
        else
            aufname = concat_media(aufname_list, dirname(newvidname))
        end
    end
    
    if MERGE_VID_AU_DYNAMIC_NORM
        cmd = `$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af loudnorm=I=-16:LRA=11:TP=-1.5 "$newvidname""_DYnormalized-audio.mp4"`
        # println(`$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af loudnorm=I=-16:LRA=11:TP=-1.5 -f matroska "$newvidname""_DYnormalized-audio.mkv"`)
        # output = @ffmpeg_env run(`$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af loudnorm=I=-16:LRA=11:TP=-1.5 "$newvidname""_DYnormalized-audio.mp4"`)
        # run(`ffmpeg -i "$newvidname" -i "$aufname" -map 0:v -map 1:a -vcodec copy -af loudnorm=I=-16:LRA=11:TP=-1.5 -f matroska "$newvidname""_normalized-audio.mkv"`)
    else
        # Get max volume and normalize
        tempfile = tempname()
        # read(pipeline(`ffmpeg -i $aufname -filter:a volumedetect -f null /dev/null`; stderr = tempfile))
        @ffmpeg_env read(pipeline(`ffmpeg -i $aufname -af astats=metadata=1 -f null /dev/null`; stderr = tempfile))

        ss = read(tempfile, String)
        rm(tempfile)
        m = eachmatch(r"Channel: (\d+)", ss)
        channel_list = [parse(Int, m.captures[1]) for m in m]
        m = eachmatch(r"Peak level dB: (.*)", ss)
        peak_db_list = [parse(Float64, m.captures[1]) for m in m]
        norm_gain = -maximum(peak_db_list[get_relevant_channels(rx_vect)])

        # m = match(r"max_volume: (.*) dB", ss)
        # max_volume = m !== nothing ? parse(Float64, m.captures[1]) : nothing
        # norm_gain = -max_volume
        cmd = `$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af "volume=$(norm_gain)dB" "$newvidname""_normalized-audio.mp4"`
        
    end
    println(cmd)
    try
        output = @ffmpeg_env run(cmd)
    catch err
        @warn "Combining video and audio failed with Julia's FFMPEG, trying again with system call"
        output = run(`$cmd -y`)
    end

    isfile("$newvidname"*"_normalized-audio.mp4") && flag_rm_oldfile && rm(newvidname) # delete video without audio
    
    aufname_list isa Array && length(aufname_list) >1 && rm(aufname) # delete concatenated audio file

    return "$newvidname"*"_normalized-audio.mp4"
    # if aufname_old isa String; rm(aufname); end
end


@info "LOADED!\tmedia_info"