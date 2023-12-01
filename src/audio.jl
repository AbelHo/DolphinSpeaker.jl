using WAV
using FileIO
using FLAC
import FLAC
using JSON
using SignalAnalysis, SignalAnalysis.Units
using FFMPEG
# require package MAT for reading .mat file
include("utils.jl")
include("media_info.jl")
include("config.jl")
using Dates


function FLAC.save(f::File{format"FLAC"}, data::Array{T,2}, samplerate; bits_per_sample = 24, compression_level = 3) where T<:Real
    @info "NEW SAVE FLAC 1003 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    
    encoder = StreamEncoderPtr()

    # Set encoder parameters
    set_compression_level(encoder, compression_level)
    num_samples = prod(size(data))
    set_total_samples_estimate(encoder, num_samples)
    set_channels(encoder, size(data)[2])
    set_sample_rate(encoder, samplerate)
    set_bits_per_sample(encoder, bits_per_sample)

    # Open file, make sure encoder was properly initialized
    initfile!(encoder, f.filename)
    if get_state(encoder) != EncoderOK
        throw(InvalidStateException("Encoder state after init_file is $(get_state(en))"))
    end

    # Shove interleaved samples into the encoder by transposing, and convert to Int32
    data_t = data'
    if eltype(data) <: AbstractFloat
        data_t = round.(Int32, data'*2^(bits_per_sample - 1))
    elseif eltype(data_t) != Int32
        data_t = Int32.(data_t)
    end

    blocksize = get_blocksize(encoder)
    for idx in 1:div(num_samples, blocksize)
        block_idxs = ((idx - 1) * blocksize + 1):(idx * blocksize)
        process_interleaved(encoder, data_t[block_idxs])
    end
    process_interleaved(encoder, data_t[end - rem(num_samples, blocksize) + 1:end])
    finish(encoder)
    return nothing
end

function readAudio(aufname::Array{String,1}; kwargs...)
    d = readAudio.(aufname; kwargs...)
    return vcat(map( x -> x[1], d)...), d[1][2], d[1][3], d[1][4], d[1][5]
end

function readAudio(aufname; fname2timestamp_func=nothing)
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
    elseif filetype == ".flac"
        data, fs = get_videos_audiodata_all(aufname)
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

function writeWAV(data, fpath; Fs=1)
    nbits = parse(Int64, string(eltype(data))[end-1:end])
    @debug nbits
    wavwrite(data, fpath; Fs=Fs, nbits=64)
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

function voltage2binary_find(vallist)
    minval = min_vals(vallist) 
    small_arg = minval .|> abs |> argmin
    # return [diff(minval); (minval[1]-diff(minval)[1]*3)]
    return [diff(minval); small_arg==1 ? minval[small_arg] : minval[small_arg] ]
end

function voltage2binary(vallist, correction)
    (vallist .- correction[2]) ./ correction[1]
end

voltage2binary_round(vallist, correction; int_type=Int16) = round.(int_type,voltage2binary(vallist, correction)) 

binary2voltage(vallist, correction) = vallist .* correction[1] .+ correction[2]

function flac2signal(aufname::String; kwargs...)
    info = @ffmpeg_env read(`$ffprobe "$aufname" -loglevel error  -v quiet -print_format json -show_format -show_streams`, String) |> JSON.parse
    # correction = JSON.parse(info["format"]["tags"]["comment"])["correction"]
    correction = JSON.parse(info["streams"][1]["tags"]["comment"])["correction"]
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

# using JSON
mat2flac(filepath, Fs, outfilepath=filepath; kwargs...) = mat2flac(filepath; Fs=Fs, outfilepath, kwargs)

# using Base64
function mat2flac(filepath; Fs=500_000, outfilepath=filepath, normalization_factor=nothing, skipdone=false)
    if isdir(filepath)
        @info "Directory! Recursively converting entire directory"
        return mat2flac.(readdir(filepath; join=true) |> skiphiddenfiles; Fs=Fs, outfilepath=outfilepath,  normalization_factor= normalization_factor, skipdone=skipdone)
        # broadcast(mat2wav, readdir(filepath; join=true) |> skiphiddenfiles, Fs,  joinpath.(Ref(outfilepath),(readdir(filepath)|> skiphiddenfiles) .*".wav") )
        # mat2wav.(filepath.*readdir(filepath); Fs=Fs, outfilepath=outfilepath)
    end
    
    filepath[end-3:end] != ".mat" && return
    filepath[end-4] == '2' && return
    @info filepath|>basename

    
    mkpath(outfilepath)
    if isnothing(outfilepath)
        outfilepath,_ = splitext(filepath)
        outfilepath = outfilepath*".flac"
    elseif isdir(outfilepath)
        @debug basename(filepath)[1:end-3] *"flac"
        outfilepath = joinpath(outfilepath, basename(filepath)[1:end-3] *"flac")
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

    # @info extrema(data, dims=1)
    @info "Finding Dynamic Range............"
    @time dd = find_DataDynamicRange_multich(data)
    @debug size(dd)
    @time correction = voltage2binary_find.(map(x->x.keys,dd))
    @debug correction
    # data_new = Array{Int}(undef, size(data))#similar(data)
    # @time Threads.@threads for ind = 1:size(data,2)
    #     data_new[:,ind] = voltage2binary_round(@view(data[:,ind]), correction[ind])
    # end

    @time data_new = hcat( voltage2binary_round.(eachcol(data), correction)...)
    comments = JSON.json(Dict("correction"=>correction)) |> string
    if size(data_new,2) < 9
        save(outfilepath, data_new, fs; bits_per_sample=16)#, raw_Int_data=true)
    
        @ffmpeg_env run(`ffmpeg -i "$outfilepath" -metadata comment="$comments" -acodec copy "$outfilepath"_meta.flac -loglevel error`)
        mv(outfilepath*"_meta.flac", outfilepath; force=true)
    else
        fnames = []
        for ind = 1:div(size(data,2), 8, RoundUp)
            max_channel = ind*8
            max_channel > size(data,2) && (max_channel = size(data,2))  

            fname = splitext(outfilepath)[1]*"_n$ind.flac"
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
        outfilepath_new = splitext(outfilepath)[1] * ".ogg"
        @ffmpeg_env run(`ffmpeg $cmds $cmds2 -c:a flac -metadata comment="$comments" $outfilepath_new -loglevel error -y`)
        rm.(fnames)
    end


    # Convert data_new to a binary string
    # binary_data = base64encode(reinterpret(UInt16, data_new));
    # # Create a pipeline of commands
    # cmds = pipeline(`echo -n $binary_data`, `base64 --decode`, `ffmpeg -f s16le -ar 500000 -ac 1 -i pipe:0 output.flac`);
    # run(cmds);


    @time redo_data = hcat(binary2voltage.(eachcol(data_new), correction)...)
    conversion_error = extrema(redo_data - data; dims=1)
    @info "Conversion Error:"
    @info conversion_error

    @debug outfilepath
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