using WAV
using FileIO
using FLAC
import FLAC
using JSON
using SignalAnalysis, SignalAnalysis.Units
using FFMPEG
# require package MAT for reading .mat file
using MAT
include("utils.jl")
include("media_info.jl")
include("config.jl")
using Dates
include("PAM.jl")


function FLAC.save(f::File{format"FLAC"}, data::Array{T,2}, samplerate; bits_per_sample = 24, compression_level = 3) where T<:Real
    # @info "NEW SAVE FLAC v2023-12-03T15:04:45.469 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    
    encoder = StreamEncoderPtr()

    # Set encoder parameters
    FLAC.set_compression_level(encoder, compression_level)
    num_samples = prod(size(data))
    FLAC.set_total_samples_estimate(encoder, num_samples)
    FLAC.set_channels(encoder, size(data)[2])
    FLAC.set_sample_rate(encoder, samplerate)
    FLAC.set_bits_per_sample(encoder, bits_per_sample)

    # Open file, make sure encoder was properly initialized
    initfile!(encoder, f.filename)
    if FLAC.get_state(encoder) != FLAC.EncoderOK
        throw(InvalidStateException("Encoder state after init_file is $(get_state(en))"))
    end

    # Shove interleaved samples into the encoder by transposing, and convert to Int32
    data_t = data'
    if eltype(data) <: AbstractFloat
        data_t = round.(Int32, data'*2^(bits_per_sample - 1))
    elseif eltype(data_t) != Int32
        data_t = Int32.(data_t)
    end

    blocksize = FLAC.get_blocksize(encoder)
    for idx in 1:div(num_samples, blocksize)
        block_idxs = ((idx - 1) * blocksize + 1):(idx * blocksize)
        FLAC.process_interleaved(encoder, data_t[block_idxs])
    end
    FLAC.process_interleaved(encoder, data_t[end - rem(num_samples, blocksize) + 1:end])
    FLAC.finish(encoder)
    return nothing
end

function readAudio(aufname::Array{String,1}; kwargs...)
    d = readAudio.(aufname; kwargs...)
    return vcat(map( x -> x[1], d)...), d[1][2], d[1][3], d[1][4], d[1][5]
end

function readAudio(aufname; fname2timestamp_func=DEFAULT_fname2timestamp_func, 
    channels=DEFAULT_bin_channels, datatype=Float64, fs=500_000)
    opt = nothing; nbits = nothing; timestamp=nothing;
    filetype = splitext(aufname)[2] |> lowercase
    if ".wav" == filetype 
        data, fs, nbits, opt = wavread(aufname, format="native")
        try
            if !isnothing(opt)
                timestamp = wav_info_read(opt)[:INAM]
            end
        catch
            @info "No timestamp in metadata"
        end
    elseif ".mat" == filetype
        d = nothing
        try
            d = load(aufname)
        catch err
            @warn "Normal .mat loading failed:"
            @warn(err)
            d = matread(aufname)
        end
        data = d["data"]
        haskey(d, "fs") ? fs=d["fs"] : fs=500_000
        if length(fs)==1; fs=fs[1]; end
        
        if aufname[end-5:end-4] == "_1"
            if isfile( splitext(aufname)[1][1:end-1]*"2.mat" )
                try
                    @debug("Loading 2nd File: "*splitext(aufname)[1][1:end-1]*"2.mat")
                    d = load(splitext(aufname)[1][1:end-1]*"2.mat")
                    data = vcat(data, d["data"])
                catch err
                    @warn "Failed to open 2nd File.............."
                    @warn err
                end
            end
        end
    elseif filetype in [".flac", ".ogg"]
        @info "Reading flac/ogg file via ffmpeg..."
        data, fs = get_videos_audiodata_all(aufname)
    elseif ".bin" == filetype
        data = loadDataBin2(aufname; channels=channels, datatype=datatype, fs=fs)
    else
        data, fs = load(aufname)
    end
    if !isnothing(fname2timestamp_func)
        try
            timestamp = fname2timestamp_func(aufname)
        catch err
            @warn "Cant convert filename to timestamp"
            @warn err
        end
    end

    return data, Int(fs), nbits, opt, timestamp
end

simple_fname2dt(aufname) = DateTime(basename(aufname)[1:17], dateformat"yyyymmdd_H.M.S")
fname2dt_dashdot(aufname) = DateTime(basename(aufname)[1:19], dateformat"yyyy-mm-dd_H.M.S")
fname2dt_ls1(aufname) = DateTime(basename(aufname)[1:15], dateformat"yyyymmdd_HHMMSS")

function writeWAV(data, fpath; Fs=1)
    nbits = parse(Int64, string(eltype(data))[end-1:end])
    @debug nbits
    wavwrite(data, fpath; Fs=Fs, nbits=nbits)
end  

#~ matlab files
mat2wav(filepath, Fs, outfilepath=filepath; kwargs...) = mat2wav(filepath; Fs=Fs, outfilepath, kwargs)

function mat2wav(filepath; Fs=500_000, outfilepath=filepath, normalization_factor=nothing, skipdone=false)
    if isdir(filepath)
        @info "Directory! Recursively converting entire directory"
        return mat2wav.(readdir(filepath; join=true) |> skiphiddenfiles; Fs=Fs, outfilepath=outfilepath,  normalization_factor= normalization_factor)
        # broadcast(mat2wav, readdir(filepath; join=true) |> skiphiddenfiles, Fs,  joinpath.(Ref(outfilepath),(readdir(filepath)|> skiphiddenfiles) .*".wav") )
        # mat2wav.(filepath.*readdir(filepath); Fs=Fs, outfilepath=outfilepath)
    end
    
    filepath[end-3:end] != ".mat" && return
    filepath[end-4] == '2' && return
    @info filepath|>basename

    
    mkpath(outfilepath)
    if isnothing(outfilepath)
        outfilepath,_ = splitext(filepath)
        outfilepath = outfilepath*".wav"
    elseif isdir(outfilepath)
        @debug basename(filepath)[1:end-3] *"wav"
        outfilepath = joinpath(outfilepath, basename(filepath)[1:end-3] *"wav")
    end
    @debug outfilepath
    if skipdone && isfile(outfilepath) 
        @info("SKIPPED: "*outfilepath);
        return
    end

    data=nothing;fs=nothing;opt=nothing;timestamp=nothing
    try
        data,fs,_,opt,timestamp = readAudio(filepath);
    catch err
        @error(err)
        return
    end
    # vars = matread(filepath)
    # if haskey(vars, "fs")
    #     Fs = vars["fs"]
    # end

    @info extrema(data, dims=1)
    if !isnothing(normalization_factor)
        data = data ./ normalization_factor
    end
    @debug outfilepath
    wavwrite(data, outfilepath; Fs=fs)
end

# data type
using DataStructures
function find_DataDynamicRange(aufname::String; kwargs...)
    data, fs, nbits, opt, timestamp = readAudio(aufname; kwargs...)
    find_DataDynamicRange(data; kwargs...)
end


function find_DataDynamicRange(data; kwargs...)
    # data, fs, nbits, opt, timestamp = readAudio(aufname; kwargs...)
    dd=OrderedDict() #Dict()
    # Threads.@threads 
    for ele in data
        if isnothing( get(dd, ele, nothing) )
            dd[ele]=1
        else
            dd[ele]+=1
        end
    end

    return dd

    # arr = [x for (x,y) in dd]
    # sort!(arr)
    # dif = vcat(@view(arr[2:end]), 0) .- arr
    # sort(dif)

    # plot(dif)
    # histogram(dif; nbins=1000)|>display

    # # median(dif)
    # # histogram( sort(dd.keys)|>diff; nbins=2^16)
    # return arr
end

function find_DataDynamicRange_multich(aufname::String; kwargs...)
    data, fs, nbits, opt, timestamp = readAudio(aufname; kwargs...)
    find_DataDynamicRange_multich(data; kwargs...)
end

function find_DataDynamicRange_multich(data; kwargs...)
    # data, fs, nbits, opt, timestamp = readAudio(aufname; kwargs...)
    # dd=OrderedDict() #Dict()
    chs = size(data,2)
    dd = [OrderedDict() for i in 1:chs]

    Threads.@threads for ch in 1:chs
        for ele in data[:,ch]
            if isnothing( get(dd[ch], ele, nothing) )
                dd[ch][ele]=1
            else
                dd[ch][ele]+=1
            end
        end
    end

    return dd

    for ch in 1:length(dd)
        histogram( sort(dd[ch].keys)|>diff; nbins=2^16, title=string(ch)) 
        xlims!(0,0.001) |> display
        # title!(string(ch))
    end
end

# 1/2^15 * 10 .* (1:10)' |> collect

using StatsBase
# h = fit(Histogram, sort(dd[1].keys)|>diff, nbins=2^16); 
# h.edges[1][sortperm(h.weights) |> reverse]

using Pipe:@pipe

min_vals(vallist) = [(filter( <(0), vallist) |> sort)[end]; (filter( >(0), vallist) |> sort)[1]]

function voltage2binary_find(vallist; int_type=Int16)#, normal_zero_centered=true)
    try # FIXME: this is a hack to deal with the case where the data is biased from 0, can be improved to do better stats analyzing the best output with the entire data value list
        length(vallist) == 1 && return [1; vallist[1]]
        minval = min_vals(vallist) 
        small_arg = minval .|> abs |> argmin
        bias = minval[small_arg] #median(vallist) #vallist[ (vallist .- (extrema(vallist) |> mean)) .|> abs |> argmin] ]
        # return [diff(minval); (minval[1]-diff(minval)[1]*3)]
        # return [diff(minval); small_arg==1 ? minval[small_arg] : minval[small_arg] ]
        dif = diff(minval)[1]
        min_dif = (vallist |> sort |> diff |> sort)[1]
        @debug dif
        @debug min_dif
        # return [abs(min_dif-dif)<1e-10 ? dif : min_dif;
        #     minval[small_arg] ]
        default_dif = 1/2^15 * 10;
        correction = [ abs(dif-default_dif)/default_dif < 0.11 ? dif : min_dif;
            bias]

        if voltage2binary( extrema(vallist) .|> abs |>maximum, correction)  > 2^(sizeof(int_type)*8 - 1)
            @warn "voltage2binary_find: dynamic range too large, trying to find better correction..."
            @info voltage2binary( extrema(vallist) .|> abs |>maximum, correction) 
            correction[2] = vallist[ (vallist .- (extrema(vallist) |> mean)) .|> abs |> argmin] # = voltage2binary_find( map(x->x.keys,dd[ch_exceed_dynamicrange]); normal_zero_centered=false)
        end
        return correction
    catch err
        # if length(vallist) == 1
        #     return [1; vallist[1]]
        # end
        @warn "bias in data from 0"
        dif = (vallist|> sort |> diff |> sort)[1]#[1:3]
        if sum(dif |> x -> x .- x[1]) == 0
            dif = dif[1]
        end

        return [ dif; vallist[length(vallist)÷2] ]
    end
end

function voltage2binary(vallist, correction)
    (vallist .- correction[2]) ./ correction[1]
end

voltage2binary_round(vallist, correction; int_type=Int16) = round.(int_type,voltage2binary(vallist, correction)) 

binary2voltage(vallist, correction) = vallist .* correction[1] .+ correction[2]

function flac2signal(aufname::String; kwargs...)
    info = @ffmpeg_env read(`ffprobe "$aufname" -loglevel error  -v quiet -print_format json -show_format -show_streams`, String) |> JSON.parse
    
    if haskey(info["format"], "tags")
        correction = JSON.parse(info["format"]["tags"]["comment"])["correction"]
    else
        correction = JSON.parse(info["streams"][1]["tags"]["comment"])["correction"]
    end
    # correction = map(x->(x["tags"]["comment"]|>JSON.parse)["correction"], info["streams"])

    data, fs = get_videos_audiodata_all(aufname);
    hcat(binary2voltage.(eachcol(data), correction)...), fs
end

# a = voltage2binary.(map(x->x.keys, dd), map(x->x.keys, dd) .|> voltage2binary_find)
# ind = 4;
# sort(a[ind]) .- round.(sort(a[ind])) .|> abs |> plot; vline!([argmin(abs.(sort(a[ind])))]); title!(string(ind))
# [filter( <(0), a[ind]) |> sort |> x->x[end-9:end] filter( >=(0), a[ind]) |> sort |> x->x[1:10]]

# binary2voltage(voltage2binary_round(@view(data[:,ind]), correct[ind]), correct[ind]) .- data[:,ind] |> maximum


volt2binary(val, x, nbits_m1) = (val - x[2]) / x[1] * 2^nbits_m1
f_convertQuantization(x,dd, nbits_m1=15) = @pipe (dd.keys .- x[2]) ./ x[1] .* 2^nbits_m1
f_dynamic_cost(x,dd, nbits_m1=15) = @pipe (dd.keys .- x[2]) ./ x[1] .* 2^nbits_m1 .|> _ .- round(_) .|> abs |> sum

# f_dynamic(dd; trial_vrange_offset=[10.0 0.1]) = optimize(x -> f_dynamic_cost(x, dd), trial_vrange_offset, ParticleSwarm()) #, GradientDescent(); autodiff = :forward)

# using Distributed
# @everywhere function f_dynamic(dd; trial_vrange_offset=[10.0 0.01], optimizer=NelderMead())
#     return optimize(x -> f_dynamic_cost(x, dd), trial_vrange_offset, optimizer)
# end


# @time dd = find_DataDynamicRange_multich(data)
# @time correction = voltage2binary_find.(map(x->x.keys,dd))

# # data_new = Array{Int}(undef, size(data))#similar(data)
# # @time Threads.@threads for ind = 1:size(data,2)
# #     data_new[:,ind] = voltage2binary_round(@view(data[:,ind]), correction[ind])
# # end

# @time data_new = hcat( voltage2binary_round.(eachcol(data), correction)...)
# save(joinpath(dirname(aufname), splitext(aufname)[1] * ".flac"), data_new, fs; bits_per_sample=16, raw_Int_data=true)


# @time redo_data = hcat(binary2voltage.(eachcol(data_new), correction)...)
# extrema(redo_data - data; dims=1)

bin2flac(args...; kwargs...) = mat2flac(args...; filetype=".bin", kwargs...)
bin2flac_check(args...; kwargs...) = mat2flac_check(args...; filetype=".bin", kwargs...)

function mat2flac_check(filepath; kwargs...)
    accum_res = [];
    mat2flac(filepath; accum_res=accum_res, kwargs...)
    showall(accum_res)
    return accum_res
end

# using JSON
mat2flac(filepath, Fs, outfilepath=filepath; kwargs...) = mat2flac(filepath; Fs=Fs, outfilepath, kwargs...)

# FIXME: implement reduction of bit-depth if dynamic range is found to be small
# using Base64
function mat2flac(filepath; Fs=500_000, outfilepath=filepath, normalization_factor=nothing, skipdone=false, binary_channel_list=nothing, remove_original=false, remove_original_errortolerance=1.4e-4, accum_res=nothing,
    filetype=".mat", int_type=Int16, metadata=(device_ID=FILE_device_ID, location_ID=FILE_location_ID, gain_setting=FILE_gain_setting, device_location=FILE_device_location), kwargs...)
    @debug "version 2024-01-23T12:00"
    GC.gc()
    if isdir(filepath)
        @info "Directory! Recursively converting entire directory"
        @info "-- " * filepath
        return mat2flac.(readdir(filepath; join=true) |> skiphiddenfiles; Fs=Fs, outfilepath=outfilepath,  normalization_factor= normalization_factor, skipdone=skipdone, binary_channel_list=binary_channel_list, remove_original=remove_original, remove_original_errortolerance=remove_original_errortolerance, accum_res=accum_res, filetype=filetype, kwargs...)
        # broadcast(mat2wav, readdir(filepath; join=true) |> skiphiddenfiles, Fs,  joinpath.(Ref(outfilepath),(readdir(filepath)|> skiphiddenfiles) .*".wav") )
        # mat2wav.(filepath.*readdir(filepath); Fs=Fs, outfilepath=outfilepath)
    end
    
    filepath[end-3:end] != filetype && return
    filepath[end-5:end-4] == "_2" && return
    @info filepath|>basename

    outfilepath == :inplace && (outfilepath=dirname(filepath))
    mkpath(outfilepath)
    if isnothing(outfilepath)
        outfilepath,_ = splitext(filepath)
        outfilepath = outfilepath*".flac"
    elseif isdir(outfilepath)
        @debug basename(filepath)[1:end-3] *"flac"
        outfilepath = joinpath(outfilepath, splitext(basename(filepath))[1] *".flac")
    end
    #~ for ophk when file is split into 2 parts. #FIXME
    if splitext(outfilepath)[1][end-2:end] == "._1"
        splitted = splitext(outfilepath)
        outfilepath = splitted[1][1:end-3] * splitted[2]
        @info "___removed '._1' from output filename"
    end

    @debug outfilepath
    if skipdone && (isfile(outfilepath) || isfile(splitext(outfilepath)[1] * ".ogg") )
        @info("SKIPPED: "*filepath);
        return
    end

    data=nothing;fs=nothing;opt=nothing;timestamp=nothing
    try
        data,fs,_,opt,timestamp = readAudio(filepath; kwargs...);
        @debug size(data), fs
        # plot(signal(data,fs)); title!("1") |> display
    catch err
        @error(err)
        return
    end
    # vars = matread(filepath)
    # if haskey(vars, "fs")
    #     Fs = vars["fs"]
    # end

    data_old=nothing; data_binchan=nothing;
    if !isnothing(binary_channel_list)
        data_old = data;
        data_binchan = @view(data_old[:,binary_channel_list])
        maxi = maximum(data_binchan; dims=1)
        for i in 1:size(data_binchan,2)
            maxi = maximum(data_binchan[:,i])
            data_binchan[:,i] = maxi==0 ? Int16.(data_binchan[:,i]) : (data_binchan[:,i] ./ maxi .|> Int16)
        end
        # data_binchan = maxi==0 ? Int16.(data_binchan) : (data_binchan ./ maxi .|> Int16)
        data = @view(data_old[:, setdiff(1:end, binary_channel_list)])
    end
    
    # @info extrema(data, dims=1)
    @info "Finding Dynamic Range............"
    @time dd = find_DataDynamicRange_multich(data)
    @debug size(dd)
    correction = voltage2binary_find.(map(x->x.keys,dd); int_type=int_type)
    @debug correction
    #~ check correction of dd extrema is still within Nbit range
    # ch_exceed_dynamicrange = findall( >(sizeof(int_type)*8), voltage2binary.( map( x-> extrema(x.keys).|>abs|>maximum, dd), correction))
    # correction[ch_exceed_dynamicrange] = voltage2binary_find.( map(x->x.keys,dd[ch_exceed_dynamicrange]); normal_zero_centered=false)

    # data_new = Array{Int}(undef, size(data))#similar(data)
    # @time Threads.@threads for ind = 1:size(data,2)
    #     data_new[:,ind] = voltage2binary_round(@view(data[:,ind]), correction[ind])
    # end
    try
        data_new = hcat( voltage2binary_round.(eachcol(data), correction)...)
    catch err
        @error "FAILED: Convertion to newdata!! Skipping......"
        @error filepath
        @error exception=(err, catch_backtrace())
        @info extrema(hcat( voltage2binary.(eachcol(data), correction)...); dims=1)
        @info correction
        return
    end
    if !isnothing(binary_channel_list)
        data_new = hcat(data_new, data_binchan)
        data = data_old
        push!(correction, [1;0])
    end
    
    comments = JSON.json(Dict("correction"=>correction)) |> string
    if size(data_new,2) < 9
        if isfile(outfilepath)
            println("File($outfilepath) already exists. Overwrite? (y/n): ")
            s = readline()
            if lowercase(s) != "y"
                @warn "not overwriting. Skipping file......"
                return
            end
            @warn "Overwriting file: $outfilepath"
        end

        outfilepath_temp = joinpath(tempdir(), basename(outfilepath * "_temp.flac"))
        @debug size(data_new), fs
        @debug typeof(data_new), typeof(fs)
        # plot(signal(data_new,fs)); title!("2") |> display
        if reduce(|, size(data_new,2) .== [5,6,7])
            wavfile_temp = joinpath(tempdir(), basename(splitext(outfilepath_temp)[1] * "_temp.wav"))
            writeWAV(data_new, wavfile_temp; Fs=fs)
            @ffmpeg_env run(`ffmpeg -i $wavfile_temp -acodec flac $outfilepath_temp -loglevel error -y`)
            rm(wavfile_temp)
            # @ffmpeg_env run(`ffmpeg -i $temp_wavfile -metadata comment="$comments" -acodec flac $outfilepath_temp -loglevel error -y`)
        else
            save(outfilepath_temp, data_new, fs; bits_per_sample=16)#, raw_Int_data=true)
        end
        # FIXME: combine the two steps above and below into one
        if isnothing(timestamp)
            if isnothing(metadata)
                @ffmpeg_env run(`ffmpeg -i "$outfilepath_temp" -metadata comment="$comments" -acodec copy "$outfilepath" -loglevel error -y`)
            else
                @ffmpeg_env run(`ffmpeg -i "$outfilepath_temp" -metadata comment="$comments" -metadata artist="$(metadata.device_ID)" -metadata genre="$(metadata.location_ID)" -metadata track="$(metadata.gain_setting)" -metadata album="$(metadata.device_location)" -acodec copy "$outfilepath" -loglevel error -y`)
            end
        else
            if isnothing(metadata)
                @ffmpeg_env run(`ffmpeg -i "$outfilepath_temp" -metadata comment="$comments" -metadata title="$timestamp" -acodec copy "$outfilepath" -loglevel error -y`)
            else
                @ffmpeg_env run(`ffmpeg -i "$outfilepath_temp" -metadata comment="$comments" -metadata title="$timestamp" -metadata artist="$(metadata.device_ID)" -metadata genre="$(metadata.location_ID)" -metadata track="$(metadata.gain_setting)" -metadata album="$(metadata.device_location)" -acodec copy "$outfilepath" -loglevel error -y`)
            end
        end
        # mv(outfilepath*"_meta.flac", outfilepath; force=true)
        try 
            rm(outfilepath_temp)
        catch err
            @warn "Failed to remove temp file"
            @warn err
        end
    else # merge into .ogg file from multiple .flac due to max 8 channels max limit for flac
        fnames = []
        for ind = 1:div(size(data_new,2), 8, RoundUp)
            max_channel = ind*8
            max_channel > size(data_new,2) && (max_channel = size(data_new,2))  

            fname = joinpath(pwd(), splitext(basename(outfilepath))[1]*"_n$ind.flac")  #splitext(outfilepath)[1]*"_n$ind.flac"
            save( fname, data_new[:,((ind-1)*8)+1:max_channel], fs; bits_per_sample=16)#, raw_Int_data=true)
            push!(fnames, fname)
        end
        # join_fnames = join(fnames, "' -i '")
        # join_fnames = join(["$(fname)'" for fname in fnames], " -i ")
        cmds=[]; cmds2=[];
        for fname in fnames
            push!(cmds, "-i")
            push!(cmds, "$fname")

            val = length(cmds2)÷2;
            push!(cmds2, "-map")
            push!(cmds2, val)
        end
        outfilepath = splitext(outfilepath)[1] * ".ogg"
        if !isnothing(timestamp)
            if !isnothing(metadata)
                @ffmpeg_env run(`ffmpeg $cmds $cmds2 -c:a flac -metadata comment="$comments" -metadata title="$timestamp" -metadata artist="$(metadata.device_ID)" -metadata genre="$(metadata.location_ID)" -metadata track="$(metadata.gain_setting)" -metadata album="$(metadata.device_location)" $outfilepath -loglevel error -y`)
            else
                @ffmpeg_env run(`ffmpeg $cmds $cmds2 -c:a flac -metadata comment="$comments" -metadata title="$timestamp" $outfilepath -loglevel error -y`)
            end
        else
            if isnothing(metadata)
                @ffmpeg_env run(`ffmpeg $cmds $cmds2 -c:a flac -metadata comment="$comments" $outfilepath -loglevel error -y`)
            else
                @ffmpeg_env run(`ffmpeg $cmds $cmds2 -c:a flac -metadata comment="$comments" -metadata artist="$(metadata.device_ID)" -metadata genre="$(metadata.location_ID)" -metadata track="$(metadata.gain_setting)" -metadata album="$(metadata.device_location)" $outfilepath -loglevel error -y`)
            end
        end
        rm.(fnames)
    end


    # Convert data_new to a binary string
    # binary_data = base64encode(reinterpret(UInt16, data_new));
    # # Create a pipeline of commands
    # cmds = pipeline(`echo -n $binary_data`, `base64 --decode`, `ffmpeg -f s16le -ar 500000 -ac 1 -i pipe:0 output.flac`);
    # run(cmds);


    redo_data = hcat(binary2voltage.(eachcol(data_new), correction)...)
    conversion_error = extrema(redo_data - data; dims=1)
    @info "Conversion Error:"
    @info conversion_error

    maxi = map( x-> abs.(x) |> maximum, conversion_error) |>  maximum
    print_color = :red; error_word = "FAILED!"
    if maxi < remove_original_errortolerance
        print_color = :green
        error_word = "Passed!"
        printstyled( "$error_word _________________________________________  $remove_original_errortolerance > error: $maxi\n"; color=print_color)
    else
        @warn "Conversion Error too large!!!!"
        printstyled( "$error_word _________________________________________  $remove_original_errortolerance < error: $maxi\n"; color=print_color)
        try
            open(outfilepath * "__error.json", "a") do io
                write(io, "{\"$outfilepath\": $conversion_error }")
            end
            mv(outfilepath, splitext(outfilepath)[1] * "__error" * splitext(outfilepath)[2]; force=true)
        catch err
            @warn "Failed to write error file"
            @warn err
        end
    end

    @debug outfilepath
    @debug "remove_original:    .........    .........."
    @debug remove_original
    if remove_original
        @info "check & Removing original file: $filepath....................."
        
        @info maxi
        if maxi < remove_original_errortolerance
            rm(filepath)
            @info "-- deleted: " * filepath
        else
            @error "Conversion Error > $remove_original_errortolerance, not removing original file: $filepath....................."
        end
    end
    isnothing(accum_res) || push!(accum_res, (filepath, maxi))
    return conversion_error #data_new, correction
end


# dd = find_DataDynamicRange_multich(data)
# extrema(data, dims=1)
# for ch in 1:length(dd)
#     histogram( sort(dd[ch].keys)|>diff; nbins=2^16, title=string(ch)) 
#     xlims!(0,0.001) |> display
#     # title!(string(ch))
# end




# r = f_dynamic.(dd; trial_vrange_offset=[10.0 0.01])
# vrange_offset_sa = Optim.minimizer.(r)

# ddo = find_DataDynamicRange(data)
# histogram( sort(ddo.keys)|>diff; nbins=2^16); xlims!(0,0.001)

# hcat((map(x -> (x.keys |> sort |> diff |> sort)[1:100]', dd))'...)


# quantum = 1/2^15 * 10
# a=map(x -> (x.keys |> sort |> diff |> sort)[1:1000], dd) .|> mean

# delta_quant = a ./ quantum
# data2 = data ./ delta_quant'
# ddo = find_DataDynamicRange(data2)
# histogram( sort(ddo.keys)|>diff; nbins=2^16); xlims!(0,0.001)



# ch=3
# h = fit(Histogram, sort(dd[ch].keys)|>diff, nbins=2^16)
# h.edges[1][sortperm(h.weights; rev=true)]
