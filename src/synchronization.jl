using SignalAnalysis, SignalAnalysis.Units
using WAV
# using FFMPEG
include("utils.jl")
include("media_info.jl")
include("audio.jl")
include("dsp.jl")
# using Plots

showall(x) = show(stdout, "text/plain", x)

function skiphiddenfiles(list)
    filter(!startswith('.') ∘  basename, list)
end

function findTrigger(data, fs; threshold_percentMAX=0.75, plot_window_inS=nothing, ref_channel=size(data,2), band_pass=nothing, argmax_len=0.5) # plot_window_inS=[-.1 .1]
    # @debug ref_channel
    # @debug abs.(data[:,ref_channel])
    !isnothing(band_pass) && (data = filter_simple(data, band_pass; fs=fs))
    trigger = findfirst(abs.(data[:,ref_channel]) .> maximum(data[:,ref_channel])*threshold_percentMAX )
    @debug trigger
    isnothing(argmax_len) || iszero(argmax_len) || (trigger = argmax(@view(data[trigger:trigger+round(Int, fs*argmax_len),ref_channel])) -1 +trigger)
    @debug trigger
    isnothing(plot_window_inS) && return trigger
    # @debug fs
    # @debug "before plot"
    plot( (1:size(data,1))./fs, data[:, ref_channel]; label="sync_channel", xlabel="Time (s)", ylabel="Amplitude")
    # @debug "after plot"
    scatter!([trigger/fs], [data[trigger,ref_channel]]; label="Sync Point")
    # @debug "after scatter"
    xlims!(trigger/fs+plot_window_inS[1], trigger/fs+plot_window_inS[2]) |> display
    # xlims!(trigger+fs*plot_window_inS[1], trigger+fs*plot_window_inS[2]) |> display
    # @debug "after xlims"
    return trigger
end

function findVidAudioBlip(fname; threshold_percentMAX=0.60, plot_window_inS=[-.1 .1], flag_savefig=false, kwargs...)
    try
        data, fs = get_videos_audiodata(fname)
        @debug (basename(fname), extrema(data))

        trigger = findTrigger(data, fs; threshold_percentMAX=threshold_percentMAX, plot_window_inS=plot_window_inS, kwargs...)
        flag_savefig!=false && savefig(joinpath(flag_savefig,basename(fname)*".png"))
        return (trigger-1)/fs
        # isnothing(plot_window_inS) && return (trigger-1)/fs
        # @debug trigger
        # plot!(trigger; title=basename(fname))
    catch err
        missing
    end
end

function findAudioBlip(fname; threshold_percentMAX=0.81, plot_window_inS=[-.1 .1], flag_savefig=false, kwargs...)
    try
        # data, fs = wavread(fname; subrange=0)
        data, fs, _, _, _ = readAudio(fname)
        try
            data = data[1:fs*4,:]
        catch
            @debug "file too short"
        end
        # try 
        #     data, fs = wavread(fname; subrange=3*fs)
        # catch err
        #     @warn "file too short"
        #     data, fs = wavread(fname)
        # end
        # @debug "in"
        trigger = findTrigger(data, fs; threshold_percentMAX=threshold_percentMAX, plot_window_inS=plot_window_inS, kwargs...)
        flag_savefig!=false && savefig(joinpath(flag_savefig,basename(fname)*".png"))
        return (trigger-1)/fs
        # isnothing(plot_window_inS) && return (trigger-1)/fs
        # plot!(; title=basename(fname))

        # plot(signal(data,fs)[:,end][1.5s:2.2s]; title=basename(fname))
        # plot!(signal(data,fs)[:,end][1.5s:2.2s] .|> abs; title=basename(fname))
        # ylims!(-0.001, .101)
        # vline!([500])
    catch err
        @warn err
        missing
    end
end

function findBlip_bothVidAudio(folname; filename_only=false, vidtype=r".mkv|.MP4|.avi|.mp4", autype=r".wav|.flac|.mp3", kwargs...)
    # folname = "/Volumes/My Passport/Bahamas_2022/2022.06.25/0002"
    flist = readdir(folname; join=true) |> skiphiddenfiles
    
    vids = filter( x -> occursin(vidtype, x), flist)
    aus = filter( x -> occursin(autype, x), flist)

    if filename_only
        trigger_times = [basename.(vids) findVidAudioBlip.(vids; kwargs...) get_duration_smart.(vids) basename.(aus) findAudioBlip.(aus; kwargs...) get_duration_smart.(aus)]
    else
        trigger_times = [findVidAudioBlip.(vids; kwargs...) findAudioBlip.(aus; kwargs...)]
    end
    return trigger_times
end

function split_vid_au(folname; vidtype=r".mkv|.MP4|.avi|.mp4", autype=r".wav|.mat|.flac|.mp3")
    flist = readdir(folname; join=true) |> skiphiddenfiles
    filter( x -> occursin(vidtype, x), flist), filter( x -> occursin(autype, x), flist)
end

find_vidVSau_sync(vidfname, aufname; band_pass=[2900 3100]) = findVidAudioBlip(vidfname; plot_window_inS=nothing, band_pass=band_pass) - findAudioBlip(aufname; plot_window_inS=nothing, band_pass=band_pass)


@info "end"

# findBlip_bothVidAudio(folname; filename_only=true) |> showall



#~ #FUTURE
# skipsomething(itr) = SkipSomething(itr)

# struct SkipSomething{T}
#     x::T
# end
# IteratorSize(::Type{<:SkipSomething}) = SizeUnknown()
# IteratorEltype(::Type{SkipSomething{T}}) where {T} = IteratorEltype(T)
# eltype(::Type{SkipSomething{T}}) where {T} = nonsomethingtype(eltype(T))

# function iterate(itr::SkipSomething, state...)
#     y = iterate(itr.x, state...)
#     y === nothing && return nothing
#     item, state = y
#     while f(item)
#         y = iterate(itr.x, state)
#         y === nothing && return nothing
#         item, state = y
#     end
#     item, state
# end

# IndexStyle(::Type{<:SkipSomething{T}}) where {T} = IndexStyle(T)
# eachindex(itr::SkipSomething) =
#     Iterators.filter(i -> @inbounds(itr.x[i]) |> f, eachindex(itr.x))
# keys(itr::SkipSomething) =
#     Iterators.filter(i -> @inbounds(itr.x[i]) |> f, keys(itr.x))
# @propagate_inbounds function getindex(itr::SkipSomething, I...)
#     v = itr.x[I...]
#     f(v) && throw(SomethingException("the value at index $I is something"))
#     v
# end