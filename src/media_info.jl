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

# """
#     choose_audio_output_codec(fname::AbstractString, streamno::Integer=0)

# Inspect `fname` (ffprobe JSON via `get_media_info`) and pick an output PCM codec that
# preserves the input encoding or rounds up to the next common PCM depth.
# Returns a tuple `(codec::String, bits::Int)` suitable for passing as `-acodec $codec`.

# Behavior:
# - If the stream reports `bits_per_raw_sample` or `bits_per_sample`, that is used.
# - Float input (`sample_fmt` containing "flt" or "dbl") maps to `pcm_f32le`/`pcm_f64le`.
# - Integer bit-depths round up to 16, 24, 32 as needed.
# - Falls back to `pcm_s16le` if no info is available.
# """
function choose_audio_output_codec(fname::AbstractString, streamno::Integer=0)
    info = get_media_info(fname)
    streams = get(info, "streams", nothing)
    streams !== nothing || throw(ArgumentError("No streams found in media info for $fname"))

    # ffprobe stream indices are zero-based in our callers; JSON array is 1-based
    idx = streamno + 1
    if idx < 1 || idx > length(streams)
        # try to pick the first audio stream as fallback
        found = findfirst(s->get(s, "codec_type", "") == "audio", streams)
        found === nothing && throw(ArgumentError("streamno out of range and no audio stream found"))
        st = streams[found]
    else
        st = streams[idx]
        if get(st, "codec_type", "") != "audio"
            found = findfirst(s->get(s, "codec_type", "") == "audio", streams)
            found !== nothing && (st = streams[found])
        end
    end

    # Prefer explicit bits info if available
    bits = nothing
    if haskey(st, "bits_per_raw_sample")
        bits = try parse(Int, st["bits_per_raw_sample"]) catch; nothing end
    elseif haskey(st, "bits_per_sample")
        bits = try parse(Int, st["bits_per_sample"]) catch; nothing end
    end

    sample_fmt = get(st, "sample_fmt", "")

    # decide codec and output bit-depth (rounded up to common PCM sizes)
    if sample_fmt != "" && occursin("dbl", sample_fmt)
        return ("pcm_f64le", 64)
    elseif sample_fmt != "" && occursin("flt", sample_fmt)
        return ("pcm_f32le", 32)
    elseif bits !== nothing
        if bits <= 16
            return ("pcm_s16le", 16)
        elseif bits <= 24
            return ("pcm_s24le", 24)
        elseif bits <= 32
            return ("pcm_s32le", 32)
        else
            # uncommon >32-bit integer: cap to 32-bit PCM signed
            return ("pcm_s32le", 32)
        end
    else
        # fallback mapping from sample_fmt when bits not present
        if occursin("s24", sample_fmt)
            return ("pcm_s24le", 24)
        elseif occursin("s32", sample_fmt) || occursin("sint32", sample_fmt)
            return ("pcm_s32le", 32)
        elseif occursin("s16", sample_fmt) || sample_fmt != ""
            return ("pcm_s16le", 16)
        else
            # ultimate fallback
            return ("pcm_s16le", 16)
        end
    end
end


"""
    get_audio_data_auto(fname::AbstractString, streamno::Union{Integer, Symbol}=0)

Read raw audio samples from `fname` automatically selecting an output PCM codec
that preserves or rounds-up input encoding (uses `choose_audio_output_codec`).

Returns a tuple `(data, samplerate)` where `data` is:
- a Vector{Int16}/Vector{Int32}/Vector{Float32}/Vector{Float64} for single-channel streams, or
- a Matrix with shape (frames, channels) for multichannel streams (same element type as above).

If `streamno == :all` the function returns a combined matrix built by horizontally
concatenating all audio streams' data (and samplerate from the first stream).
"""
function get_audio_data_auto(fname::AbstractString, streamno::Union{Integer, Symbol}=0) #FIXME: haven't check, might need work on 24 bits, or make it more useful
    if streamno == :all
        info = get_media_info(fname)
        # pick only audio stream indices (0-based)
        aud_idxs = [i-1 for (i,s) in enumerate(info["streams"]) if get(s, "codec_type", "") == "audio"]
        d = get_audio_data_auto.(Ref(fname), aud_idxs)
        return hcat(map(x->x[1], d)...), d[1][2]
    end

    # pick codec + bit depth to request
    codec, bits = choose_audio_output_codec(fname, streamno)
    # derive raw format token ffmpeg expects (e.g. s16le, s24le, s32le, f32le, f64le)
    fmt = occursin("pcm_", codec) ? replace(codec, "pcm_" => "") : codec

    # read raw bytes from ffmpeg
    raw = @ffmpeg_env read(`$ffmpeg -i $fname -map 0:$streamno -f $fmt -acodec $codec -loglevel error -`)

    # convert bytes to typed samples
    data = nothing
    if bits == 16
        data = reinterpret(Int16, raw)
    elseif bits == 32 && occursin("pcm_f", codec)
        data = reinterpret(Float32, raw)
    elseif bits == 64 && occursin("pcm_f", codec)
        data = reinterpret(Float64, raw)
    elseif bits == 32
        # 32-bit integer PCM
        data = reinterpret(Int32, raw)
    elseif bits == 24
        # 3 bytes per sample, little-endian signed 24-bit -> sign-extend to Int32
        bytes = Vector{UInt8}(raw)
        nsamples = length(bytes) ÷ 3
        out = Array{Int32}(undef, nsamples)
        @inbounds for i in 1:nsamples
            b1 = Int32(bytes[3i-2])
            b2 = Int32(bytes[3i-1]) << 8
            b3 = Int32(bytes[3i])   << 16
            v = b1 | b2 | b3
            # sign-extend 24->32
            if (v & 0x800000) != 0
                v |= Int32(0xFF000000)
            end
            out[i] = v
        end
        data = out
    else
        # fallback: try to reinterpret as Int16
        try
            data = reinterpret(Int16, raw)
        catch
            data = reinterpret(UInt8, raw)
        end
    end

    # reshape to frames x channels if multichannel
    ch = get_whatever(fname, streamno, "a"; entries_custom="channels")
    if !(ch === nothing) && ch |> Int > 1
        channels = Int(ch)
        # ensure length is divisible by channels
        nsamps = length(data) ÷ channels
        data = reshape(data, channels, nsamps)'  # frames x channels
    end

    return data, get_framerate(fname, streamno, "a")
end


## FIXME: output codec to automatically convert to raw form
function get_videos_audiodata_all2(vidfname, streamno=:all; data_format=Int16)
	if streamno == :all
		info = get_media_info(vidfname)
		streamno = 0:length(info["streams"])-1
		d = get_videos_audiodata_all2.(Ref(vidfname), streamno)
		return hcat(map(x->x[1],d)...), d[1][2]
	end
	codec = choose_audio_output_codec(vidfname, streamno)
	@info codec
	# choose a compatible Julia element type (data_format) based on the chosen codec,
	# and adjust the codec request if necessary (e.g. promote 24-bit -> 32-bit for easy reinterpret)
	cname, cbits = codec
	cname_l = lowercase(String(cname))
	# try to extract alphabetic token and numeric bits from codec name with regex
	m = match(r"^([a-z_]+).*?(\d+)", cname_l)
	alpha = m === nothing ? "" : m.captures[1]
	# prefer explicit cbits when present, otherwise take from codec name
	bits = cbits !== nothing ? cbits : (m !== nothing ? try parse(Int, m.captures[2]) catch; nothing end : nothing)

	# special-case 24-bit: promote to 32-bit PCM for byte alignment
	if bits == 24 || occursin(r"s24", cname_l)
		@info "24-bit PCM requested — promoting to 32-bit PCM output for byte alignment"
		cname = "pcm_s32le"
		codec = (cname, 32)
		data_format = Int32
	else
		# detect float codecs by alpha token or name containing 'f' / 'flt' / 'dbl'
		isfloat = occursin(r"\bf|flt|dbl|float\b", cname_l) || occursin("pcm_f", cname_l)
		if isfloat
			# map float bit-depths
			if bits == 64 || occursin("f64", cname_l)
				data_format = Float64
			else
				# default to 32-bit float for float-like codecs
				data_format = Float32
			end
		else
			# integer PCM: map bit-depth to appropriate Int type (common cases)
			if bits === nothing
				@warn "Could not detect bit depth for codec $codec; defaulting to pcm_s16le / Int16"
				cname = "pcm_s16le"
				codec = (cname, 16)
				data_format = Int16
			elseif bits <= 16
				data_format = Int16
			else
				# anything >16 (and not 24 which was handled) -> use 32-bit integer container
				data_format = Int32
			end
		end
	end
	@info "Using data format $data_format for codec $codec"
	
	strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -map 0:$streamno -f $(split(codec[1],'_')[2]) -acodec $(codec[1]) -loglevel error -`)


	# strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -map 0:$streamno -f s16le -acodec pcm_s16le -loglevel error -`)# .|> Int16;
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
	vid_audiodata = reinterpret(data_format, strs)#, get_framerate(vidfname, 0, "a")
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
function concat_media(filelist::Vector{String}, outputfile_dir::String;
	flag_output_auto=true, flag_overwrite=false, flag_return_cmd=false,
	force_extension_type=nothing, kwargs...)

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

	# use the concat demuxer as input 0 and add the first file as input 1,
	# then copy streams from input 0 and copy metadata from input 1
	firstfile = filelist[1]
	cmd = `ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i $tmpfile -i $firstfile -c copy -map 0 -map_metadata 1 $outputfile`
	# _cmd = `ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i $tmpfile -i $firstfile -c copy -map 0 -map_metadata 1 -map_chapters 1 -movflags use_metadata_tags $outputfile`
	# remux if requested extension differs from input
	if !isnothing(force_extension_type)
		# normalize extension to start with a dot
		# newext = String(force_extension_type)
		# startswith(newext, ".") || (newext = "." * newext)
		if force_extension_type != ext
			outputfile = splitext(outputfile)[1] * force_extension_type
			cmd = `ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i $tmpfile -i $firstfile -map 0 -map_metadata 1 $outputfile`
			@info "Forced Ext: $cmd"
			# `ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i $tmpfile -i $firstfile -c copy -map 0 -map_metadata 1 $new_outputfile`
			# outputfile = new_outputfile
		end
	end

	out = ""
	# run the command
	try
		if flag_overwrite
			@info "Running ffmpeg with overwrite enabled: $cmd -y"
			out = @ffmpeg_env read(`$cmd -y`)
			# buf = IOBuffer()
			# run(pipeline(`$cmd -y 1>&2`, stdout=buf, stderr=buf))
			# out = String(take!(buf))
		else
			out = @ffmpeg_env read(cmd)
		end
	catch err
		@warn "Concatenation failed with Julia's FFMPEG, trying again with system call: $err"
		out = read(`$cmd -y`)
	finally
		# delete the temporary file
		rm(tmpfile; force=true)
	end
	if flag_return_cmd
		return outputfile, out
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
function combine_vidau(newvidname, aufname_list; vidau_syncdiff=0, MERGE_VID_AU_DYNAMIC_NORM=false, rx_vect=rx_vect, 
	flag_rm_oldfile=false, flag_rm_concataudio=false, kwargs...) #TODO: swap audio channel according to rx_vect location to correspond Left, Right, Center
    if aufname_list isa Array
        if length(aufname_list)==1
            aufname = aufname_list[1]
        else
            aufname = concat_media(aufname_list, dirname(newvidname); kwargs...)
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
    
    flag_rm_concataudio && aufname_list isa Array && length(aufname_list) >1 && rm(aufname) # delete concatenated audio file

    return "$newvidname"*"_normalized-audio.mp4"
    # if aufname_old isa String; rm(aufname); end
end


@info "LOADED!\tmedia_info"