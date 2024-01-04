# using SignalAnalysis, SignalAnalysis.Units
# using WAV
# using FFMPEG
# using Dates
# import TimeZones
# include("audio.jl")
import FileIO.load

showall(x) = show(stdout, "text/plain", x)

function skiphiddenfiles(list)
    filter(!startswith('.') ∘  basename, list)
end

function savejld(savefname; kwargs...)
    # savefname = splitext(basename(aufname))[1] *"_t"*string(impulsive_autothreshold_median_ratio)*"_d"*string(dist_impulsive)*".jld2"
    if !Sys.islinux()
        jldsave(savefname; kwargs...)
        # jldsave(joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".jld2"); pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect,  vidfname, aufname, res_dir, detector_set, pixel_related)
    #@error(err)
    #@error("Cant save jld2 file, saving locally and copying instead")
    #rm(joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".jld2"))
    else
        jldsave(joinpath("", savefname|>basename); kwargs...)
        mv(joinpath("", savefname|>basename),
            savefname)
    end
    @info "Saved to: " * savefname
end

function FileIO.load(filenames::Array{String,1}; kwargs...)
    merge( load.(filenames; kwargs...)...)
end

# function findTrigger(data, fs; threshold_percentMAX=0.75, plot_window_inS=nothing, ref_channel=size(data,2)) # plot_window_inS=[-.1 .1]
#     # @debug ref_channel
#     # @debug abs.(data[:,ref_channel])
#     trigger = findfirst(abs.(data[:,ref_channel]) .> maximum(data[:,ref_channel])*threshold_percentMAX )
#     # @debug trigger
#     isnothing(plot_window_inS) && return trigger
#     # @debug fs
#     plot(data[:, ref_channel])
#     @debug "after plot"
#     scatter!([trigger], [data[trigger,ref_channel]])
#     xlims!(trigger+fs*plot_window_inS[1], trigger+fs*plot_window_inS[2])
#     # return p
# end

# function findVidAudioBlip(fname; threshold_percentMAX=0.60, plot_window_inS=[-.1 .1])
#     try
#         data, fs = get_videos_audiodata(fname)
#         @debug (basename(fname), extrema(data))

#         trigger = findTrigger(data, fs; threshold_percentMAX=threshold_percentMAX, plot_window_inS=plot_window_inS)
        
#         isnothing(plot_window_inS) && return (trigger-1)/fs
#         # @debug trigger
#         plot!(trigger; title=basename(fname))
#     catch err
#         missing
#     end
# end

# function findAudioBlip(fname; threshold_percentMAX=0.81, plot_window_inS=[-.1 .1])
#     try
#         if occursin(".wav", fname)
#             data, fs = wavread(fname; subrange=0)
#             try 
#                 data, fs = wavread(fname; subrange=3*fs)
#             catch err
#                 @warn "file too short"
#                 data, fs = wavread(fname)
#             end
#         else
#             data, fs, nbits, opt, timestamp = readAudio(fname);
#             data = data[1:fs*3,:]
#         end

#         trigger = findTrigger(data, fs; threshold_percentMAX=threshold_percentMAX, plot_window_inS=plot_window_inS)

#         isnothing(plot_window_inS) && return (trigger-1)/fs
#         plot!(; title=basename(fname))

#         # plot(signal(data,fs)[:,end][1.5s:2.2s]; title=basename(fname))
#         # plot!(signal(data,fs)[:,end][1.5s:2.2s] .|> abs; title=basename(fname))
#         # ylims!(-0.001, .101)
#         # vline!([500])
#     catch err
#         missing
#     end
# end

# function findBlip_bothVidAudio(folname; filename_only=false, vidtype=r".mkv|.MP4|.avi|.mp4", autype=r".wav|.mat|.flac.mp3")
#     # folname = "/Volumes/My Passport/Bahamas_2022/2022.06.25/0002"
#     flist = readdir(folname; join=true) |> skiphiddenfiles
    
#     vids = filter( x -> occursin(vidtype, x), flist)
#     aus = filter( x -> occursin(autype, x), flist)

#     if filename_only
#         trigger_times = [basename.(vids) findVidAudioBlip.(vids; plot_window_inS=nothing) get_duration_smart.(vids) basename.(aus) findAudioBlip.(aus; plot_window_inS=nothing) get_duration_smart.(aus)]
#     else
#         trigger_times = [findVidAudioBlip.(vids; plot_window_inS=nothing) findAudioBlip.(aus; plot_window_inS=nothing)]
#     end
#     return trigger_times
# end


# function TimeZones.ZonedDateTime(dt::AbstractString)
#     plus_index = findlast('+', dt)
#     tz_ind = isnothing(plus_index) ? findlast('-',dt) : plus_index
#     return ZonedDateTime( DateTime(dt[1:tz_ind-1]), TimeZone(dt[tz_ind:end]))
# end

# function TimeZones.ZonedDateTime(dt::String)
# # function zdt(dt::String)
#     if '.' in dt
#         if length(dt) > 27
#             plus_index = findlast('+', dt)
#             tz_ind = isnothing(plus_index) ? findlast('-',dt) : plus_index

#             ts =  split(dt, ['-','T',':','.','+'])
#             return ZonedDateTime( parse.(Int, ts[1:7])..., TimeZone(dt[tz_ind:end]))
#         end
#     end

#     # ZonedDateTime( (parse.(Int,split(ts[1:23], ['-','T',':','.'])))..., TimeZone("-0400"))
# end

@info "end utils.jl"

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