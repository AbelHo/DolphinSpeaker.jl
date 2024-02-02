using Pipe:@pipe
# using FFMPEG
# include("detector_impulsive.jl")
# include("localization.jl")
# threshold_impulsive = nothing
# include("detector_tonal.jl")

if false
    aufname = "/Volumes/One Touch/calf/clicker_calib/20231005/20231005_16.27.34_log._1.mat"
    vidfname = "/Volumes/One Touch/calf/clicker_calib/20231005/20231005_16.27.34_log.mkv"
    res_dir = "/Volumes/One Touch/res/calf/clickertest/20231005"

    aufname = "/Volumes/One Touch/calf/clicker_calib/20231005/20231005_16.23.55_log.flac"


    analyse_dir = "/Volumes/One Touch/calf/clicker_calib/20211214"
    res_dir = "/Volumes/One Touch/res/calf/clickertest/20211214"
    analyse_dir = "/Volumes/One Touch/calf/coop/20231120"
    res_dir = "/Volumes/One Touch/res/calf/coop/20231120"

    auftype = r".wav|.flac";
    mkpath(res_dir)

    aufnames = @pipe readdir(analyse_dir; join=true) |> skiphiddenfiles |> filter( x-> occursin(auftype, x), _)
    res_impulses = detect_impulse.(aufnames, Ref(res_dir); return_datafilt=true)
    angle_list = map(res_impulse -> detection2angle(res_impulse.data_filt, res_impulse.pind_good, rx_vect; fs=res_impulse.fs, return_residual=true, window=window_impulsive), res_impulses)

    plotlyjs()
    flag_savefig = res_dir #nothing
    for ind = 1:length(angle_list)
        scatter(res_impulses[ind].pind_good_inS, angle_list[ind][1][1].|>rad2deg; markershape=:xcross, alpha=res_impulses[ind].ppeak.^2/maximum(res_impulses[ind].ppeak.^2), labels=["azimuth" "inclination"]); 
        title!(aufnames[ind]|>basename; xlabel="Time(s)", ylabel="Angle(°)") |> display
        isnothing(flag_savefig) && continue;
        savefig(joinpath(res_dir, splitext(aufnames[ind]|>basename)[1]*"_angles.png"))
        savefig(joinpath(res_dir, splitext(aufnames[ind]|>basename)[1]*"_angles.html"))

        scatter(angle_list[ind][1][1].|>rad2deg; markershape=:xcross, alpha=res_impulses[ind].ppeak.^2/maximum(res_impulses[ind].ppeak.^2), labels=["azimuth" "inclination"]); 
        title!(aufnames[ind]|>basename; xlabel="Detection Index", ylabel="Angle(°)") |> display
        savefig(joinpath(res_dir, splitext(aufnames[ind]|>basename)[1]*"_angles-index.png"))
        savefig(joinpath(res_dir, splitext(aufnames[ind]|>basename)[1]*"_angles-index.html"))
    end

    ind = 1
    res = res_impulses[ind]
    aufname = aufnames[ind]
    vidfname = splitext(aufname)[1]*".mkv"
    include("plotting.jl"); include("readImages.jl")
    using VideoIO
    plot_all_clicks(joinpath(res_dir,splitext(aufname|>basename)[1]), vidfname, res.pind_good_inS, res.data_filt, res.pind_good, angle_list[ind][2], window_impulsive, angle_list[ind][1][1]; plotsize=(1280,720))



    scatter(angs[1].|>rad2deg; markershape=:xcross, alpha=res_impulse.ppeak.^2/maximum(res_impulse.ppeak.^2))

end

################### ^  ###################  ^ ################### ^  ###################   ###################   
# # aufname = "/Users/abel/Documents/data/calf/coop/20231120/20231120_14.20.14_log.flac"
# vidfname = splitext(aufname)[1]*".mkv"
# # res_dir = "/Users/abel/Documents/data_res/calf/coop"
# # res_impulse = res_impulses[2]
# mkpath(res_dir)

# res_impulse = detect_impulse(aufname, res_dir; return_datafilt=true); #band_pass=[500 Inf], dist=80, threshold=.01)
# res_impulsetrain = detect_impulsetrain(res_impulse, res_dir);
# res_tonal = detect_tonal(aufname, res_dir)
# res = res_impulse;

# include("localization.jl")
# angs, tdoas = detection2angle(res_impulse.data_filt, res_impulse.pind_good, rx_vect; fs=res_impulse.fs, return_residual=true, window=window_impulsive)#, ref_channel=ref_channel, channels_relevant=1:size(rx_vect,2),
# # getTDOA_func=get_tdoa_raw, solver_func=default_tdoa2dir_solver, cost_tdoa2ang=cost_tdoa2ang, return_residual=false)

# scatter(res_impulse.pind_good_inS, angs[1].|>rad2deg; markershape=:xcross, alpha=res_impulse.ppeak.^2/maximum(res_impulse.ppeak.^2), labels=["azimuth" "inclination"]); 
# title!(aufname|>basename; xlabel="Time(s)", ylabel="Angle(°)") |> display


# ## using Plots
# ## p1 = plot(angs[1] .|> rad2deg);
# ## p2 = plot(angs[2]);
# ## p3 = plot(res_impulse.ppeak);
# ## plot(p1, p2, p3, layout=(3,1)) 
# #
# ## plot(angs .|> rad2deg)#, tdoas, seriestype=:scatter, xlabel="angle", ylabel="tdoa", legend=false)
# ## plot(tdoas[:,1:2])#, seriestype=:scatter, xlabel="angle", ylabel="tdoa", legend=false)

# include("plotting.jl"); include("readImages.jl")
# using VideoIO
# outfol = joinpath(res_dir,splitext(aufname|>basename)[1])
# plot_all_clicks(outfol, vidfname, res.pind_good_inS, res.data_filt, res.pind_good, tdoas, window_impulsive, angs[1]; plotsize=(1280,720))

# outfol_img = joinpath(outfol, "%d.png")
# outfol_clicks = outfol * "_clicks.mp4"

# @ffmpeg_env run(`ffmpeg -framerate 30 -i $outfol_img -c:v libx264 -r 30 -pix_fmt yuv420p $outfol_clicks`)

# plotlyjs()
# scatter(res_impulse.pind_good_inS, angs[1].|>rad2deg; markershape=:xcross, alpha=res_impulse.ppeak.^2/maximum(res_impulse.ppeak.^2), labels=["azimuth" "inclination"]); 
# p=title!(aufname|>basename; xlabel="Time(s)", ylabel="Angle(°)");# |> display
# savefig(outfol * "_clicks.png")
# savefig(outfol * "_clicks.html")
# display(p)


if false

include("run_example.jl");include("detector_impulsive.jl")
aufname = "/Volumes/One Touch/calf/coop/20231120/20231120_14.20.14_log._1.mat"
vidfname = "/Volumes/One Touch/calf/coop/20231120/20231120_14.20.14_log.mkv"
res_dir = "/Volumes/One Touch/res/calf/coop"
band_pass = tonal_band_pass; threshold_tonal = -60; threshold_impulsive = 0.050800082760299394

process_one_set(vidfname, aufname, res_dir)


### showcase
findAudioBlip("/Volumes/One Touch/calf/clicker_calib/20231006/20231006_11.37.13_log.flac"; plot_window_inS=[-1 1], band_pass=[2900 3100])
title!("Acoustic Synchronization Channel")
savefig("Acoustic_sync.html")
savefig("Acoustic_sync.png")

findVidAudioBlip("/Volumes/One Touch/calf/clicker_calib/20231006/20231006_11.37.13_log.mkv"; band_pass=[2900 3100], plot_window_inS=[-1 1])
savefig("Video_sync.html")
savefig("Video_sync.png")

mat2flac("/Volumes/One Touch/calf/clicker_calib"; skipdone=true)
mat2flac("/Volumes/One Touch/calf/coop/20231120"; skipdone=true)

tt = findBlip_bothVidAudio.(readdir("/Volumes/One Touch/calf/clicker_calib";join=true)[4:end] |> skiphiddenfiles; plot_window_inS=[-1 1], filename_only=true, band_pass=[2900 3100], flag_savefig="/Volumes/One Touch/res/calf/clickertest/sync2")

end

include("localization.jl")
include("media_info.jl")
include("detector.jl")
include("audio.jl")
include("plotting.jl")

# aufname = "/Users/abel/Documents/data/calf/coop/20231128/20231128_15.16.48_log.flac"
# vidfname = splitext(aufname)[1]*".mkv"
# res_dir = "/Users/abel/Documents/data_res/calf/coop/20231128_test-20231215_tt2"
function process_detections(aufname, vidfname; res_dir=nothing)
    processed_skip_flag = false
    @info aufname
    data, fs, nbits, opt, timestamp = readAudio(aufname)
    res = detect_impulseNtonal((aufname, data, fs, timestamp), res_dir; threshold_tonal=nothing, processed_skip_flag=processed_skip_flag, return_datafilt=true);
    # isnothing(res) && (@info("___no results found!..."); return nothing)

    @info "Calculating angles..."
    #~ impulsive to angle to pixel
    res_new = res.res_impulsetrain
    window = window_impulsive
    # detection2angle(res_new.data_filt, res_impulse.pind_good, rx_vect; fs=res_new.fs, return_residual=true, window=window)
    # angs = detection2angle(res_new.data_filt, res_new.pind_good, rx_vect; fs=fs, window=window, return_residual=true)
    angs = detection2angle(res.res_impulse.data_filt, res_new.pind_good, rx_vect; fs=fs, window=window, return_residual=true)
    # res = (;  Base.structdiff(res, NamedTuple{(:res_impulse,)})..., res_impulse=Base.structdiff(res.res_impulse, NamedTuple{(:data_filt,)}))
    ang = angs[1][1]
    # plot_ang(res_new, ang; label=["azimuth" "inclination"], type="Power")
    # for list
    # angle_list = map(res_impulse -> detection2angle(res_impulse.data_filt, res_impulse.pind_good, rx_vect; fs=res_new.fs, return_residual=true, window=window), res_new)

    # convert to pixel
    p_pixels = angle2px(ang, fov_angle)
    pind_vidframes = round.(Int, res_new.pind_good_inS * get_fps(vidfname)) .+ 1

    # writedlm( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".csv"), ["p_pixel" "px" "py"], ',')
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_Impulse_t"*string(res.res_impulse.threshold)*"_d"*string(res.res_impulse.dist)*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes p_pixels], ',')
    end

    [pind_vidframes, p_pixels]



    #~ impulse xcorr
    res_new = res.res_impulsetrain
    window = -130:130#window_impulsive
    angs2 = detection2angle(res.res_impulse.data_filt, res_new.pind_good, rx_vect;
        fs=fs, window=window, return_residual=true,
        getTDOA_func=get_tdoa_raw_MaxEnergyRefChannel)
    ang2 = angs2[1][1]

    # convert to pixel
    p_pixels2 = angle2px(ang2, fov_angle)
    pind_vidframes2 = round.(Int, res_new.pind_good_inS * get_fps(vidfname)) .+ 1

    pixel_related_impulsive2 = [pind_vidframes2, p_pixels2]
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_Impulse2_t"*string(res.res_impulse.threshold)*"_d"*string(res.res_impulse.dist)*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, pixel_related_impulsive2, ',')
    end

    #~ impulse lowpassed hilbert
    res_new = res.res_impulsetrain
    window = window_impulsive#-130:130#
    angs3 = detection2angle(res.res_impulse.data_filt, res_new.pind_good, rx_vect;
        fs=fs, window=window, return_residual=true,
        getTDOA_func=get_tdoa_envelope_filtered)
    ang3 = angs3[1][1]

    # convert to pixel
    p_pixels3 = angle2px(ang3, fov_angle)
    pind_vidframes3 = round.(Int, res_new.pind_good_inS * get_fps(vidfname)) .+ 1

    pixel_related_impulsive3 = [pind_vidframes3, p_pixels3]
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_Impulse3_t"*string(res.res_impulse.threshold)*"_d"*string(res.res_impulse.dist)*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, pixel_related_impulsive3, ',')
    end


    res = (;  Base.structdiff(res, NamedTuple{(:res_impulse,)})..., res_impulse=Base.structdiff(res.res_impulse, NamedTuple{(:data_filt,)}))



    #~ tonal to angle to pixel
    res_new = res.res_tonalsegment
    windows_tonal = map( i -> res_new.train_start[i]:res_new.train_end[i], eachindex(res_new.train_start))
    data_filt = filter_simple(data, tonal_band_pass; fs=fs)
    ang_tonal, tdoas_tonal = detection2angle(data_filt, windows_tonal, rx_vect;
        getTDOA_func = get_tdoa_raw_flexi,        
        fs=fs, window=window_impulsive)#, getTDOA_func=get_tdoa_max)

    # plot_ang(res_new, ang_tonal; label=["azimuth" "inclination"], type="Power")

    # convert to pixel
    p_pixels_tonal = angle2px(ang_tonal, fov_angle)
    pind_vidframes_tonal = round.(Int, res_new.pind_good_inS * get_fps(vidfname)) .+ 1

    res_new = res.res_tonal
    # writedlm( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".csv"), ["p_pixel" "px" "py"], ',')
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_Tonal_t"*string(res_new.threshold)*"tonal_bp"*string(res_new.band_pass[1])*"_"*string(res_new.band_pass[2])*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes_tonal p_pixels_tonal], ',')
    end

    #~ tonal short window
    ang_tonal_short, tdoas_tonal_short = detection2angle(data_filt, windows_tonal, rx_vect;
        getTDOA_func = get_tdoa_raw_flexi,        
        fs=fs, window=window_impulsive)#, getTDOA_func=get_tdoa_max)
    ang_tonal_short = detection2angle(data_filt, res_new.pind_good, rx_vect; fs=fs, window=window, return_residual=true)
    ang_tonal_short = ang_tonal_short[1][1]

    # plot_ang(res_new, ang_tonal; label=["azimuth" "inclination"], type="Power")

    # convert to pixel
    p_pixels_tonal_short = angle2px(ang_tonal_short, fov_angle)
    pind_vidframes_tonal_short = round.(Int, res_new.pind_good_inS * get_fps(vidfname)) .+ 1
    
    pixel_related_tonal_short = [pind_vidframes_tonal_short, p_pixels_tonal_short]

    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_Tonal-short_t"*string(res_new.threshold)*"tonal_bp"*string(res_new.band_pass[1])*"_"*string(res_new.band_pass[2])*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes_tonal_short p_pixels_tonal_short], ',')
    end


    pixel_related_impulsive = [pind_vidframes, p_pixels]
    pixel_related_tonal = [pind_vidframes_tonal, p_pixels_tonal]

    @info res_dir
    @info res.fname
    if !isnothing(res_dir)
        savejld(joinpath(res_dir, splitext(res.fname)[1]*"_angles.jld2"); ang_impulsive=angs, ang_tonal=ang_tonal, ang_tonal_short=ang_tonal_short,
            pixel_related_impulsive=pixel_related_impulsive, pixel_related_tonal=pixel_related_tonal, pixel_related_tonal_short=pixel_related_tonal_short)
    end
    
    # return [pixel_related_tonal, pixel_related_impulsive, pixel_related_tonal_short]
    return [pixel_related_tonal, pixel_related_impulsive, pixel_related_impulsive2, pixel_related_impulsive3]

end

function detectNangles(aufname, vidfname; res_dir=nothing)
    processed_skip_flag = false
    @info aufname
    data, fs, nbits, opt, timestamp = readAudio(aufname)
    res = detect_impulseNtonal((aufname, data, fs, timestamp), res_dir; threshold_tonal=nothing, processed_skip_flag=processed_skip_flag, return_datafilt=true);
    # isnothing(res) && (@info("___no results found!..."); return nothing)

    @info "Calculating angles..."
    #~ impulsive to angle to pixel
    res_new = res.res_impulsetrain
    window = window_impulsive
    # detection2angle(res_new.data_filt, res_impulse.pind_good, rx_vect; fs=res_new.fs, return_residual=true, window=window)
    # angs = detection2angle(res_new.data_filt, res_new.pind_good, rx_vect; fs=fs, window=window, return_residual=true)
    angs = detection2angle(res.res_impulse.data_filt, res_new.pind_good, rx_vect; fs=fs, window=window, return_residual=true)
    res = (;  Base.structdiff(res, NamedTuple{(:res_impulse,)})..., res_impulse=Base.structdiff(res.res_impulse, NamedTuple{(:data_filt,)}))
    ang = angs[1][1]
    # plot_ang(res_new, ang; label=["azimuth" "inclination"], type="Power")
    # for list
    # angle_list = map(res_impulse -> detection2angle(res_impulse.data_filt, res_impulse.pind_good, rx_vect; fs=res_new.fs, return_residual=true, window=window), res_new)
    ang, res_new.pind_good_inS, splitext(basename(aufname))[1] *"_Impulse_t"*string(res.res_impulse.threshold)*"_d"*string(res.res_impulse.dist)
    


    res_new = res.res_tonalsegment
    windows_tonal = map( i -> res_new.train_start[i]:res_new.train_end[i], eachindex(res_new.train_start))
    data_filt = filter_simple(data, tonal_band_pass; fs=fs)
    ang_tonal, tdoas_tonal = detection2angle(data_filt, windows_tonal, rx_vect;
        getTDOA_func = get_tdoa_raw_flexi,        
        fs=fs, window=window_impulsive)#, getTDOA_func=get_tdoa_max)

    # plot_ang(res_new, ang_tonal; label=["azimuth" "inclination"], type="Power")

    # convert to pixel
    p_pixels_tonal = angle2px(ang_tonal, fov_angle)
    pind_vidframes_tonal = round.(Int, res_new.pind_good_inS * get_fps(vidfname)) .+ 1

    res_new = res.res_tonal
    # writedlm( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".csv"), ["p_pixel" "px" "py"], ',')
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_Tonal_t"*string(res_new.threshold)*"tonal_bp"*string(res_new.band_pass[1])*"_"*string(res_new.band_pass[2])*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes p_pixels], ',')
    end

    pixel_related_impulsive = [pind_vidframes, p_pixels]
    pixel_related_tonal = [pind_vidframes_tonal, p_pixels_tonal]

    return [ang, res.res_impulsetrain.pind_good_inS, splitext(basename(aufname))[1] *"_Impulse_t"*string(res.res_impulse.threshold)*"_d"*string(res.res_impulse.dist),
            ang_tonal, res.res_tonal.pind_good_inS, splitext(basename(aufname))[1] *"_Tonal_t"*string(res_new.threshold)*"tonal_bp"*string(res_new.band_pass[1])*"_"*string(res_new.band_pass[2])]
    [pixel_related_tonal, pixel_related_impulsive]

end

function ang2pixel(ang, pind_good_inS, vidfname, outfname)
    # convert to pixel
    p_pixels = angle2px(ang, fov_angle)
    pind_vidframes = round.(Int, pind_good_inS * get_fps(vidfname)) .+ 1

    # writedlm( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".csv"), ["p_pixel" "px" "py"], ',')
    open( joinpath(res_dir, outfname*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes p_pixels], ',')
    end

    pind_vidframes, p_pixels
end

# angs, tdoas = detection2angle(res_impulse.data_filt, res_impulse.pind_good, rx_vect; fs=res_impulse.fs, return_residual=true, window=window_impulsive)#, ref_channel=ref_channel, channels_relevant=1:size(rx_vect,2),
# # getTDOA_func=get_tdoa_raw, solver_func=default_tdoa2dir_solver, cost_tdoa2ang=cost_tdoa2ang, return_residual=false)

# scatter(res_impulse.pind_good_inS, angs[1].|>rad2deg; markershape=:xcross, alpha=res_impulse.ppeak.^2/maximum(res_impulse.ppeak.^2), labels=["azimuth" "inclination"]); 
# title!(aufname|>basename; xlabel="Time(s)", ylabel="Angle(°)") |> display

