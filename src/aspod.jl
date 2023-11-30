using Pipe:@pipe
using FFMPEG
include("detector_impulsive.jl")
include("localization.jl")
# threshold_impulsive = nothing
include("detector_tonal.jl")

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
# aufname = "/Users/abel/Documents/data/calf/coop/20231120/20231120_14.20.14_log.flac"
vidfname = splitext(aufname)[1]*".mkv"
# res_dir = "/Users/abel/Documents/data_res/calf/coop"
# res_impulse = res_impulses[2]
mkpath(res_dir)

res_impulse = detect_impulse(aufname, res_dir; return_datafilt=true); #band_pass=[500 Inf], dist=80, threshold=.01)
res_impulsetrain = detect_impulsetrain(res_impulse, res_dir);
res_tonal = detect_tonal(aufname, res_dir)
res = res_impulse;

include("localization.jl")
angs, tdoas = detection2angle(res_impulse.data_filt, res_impulse.pind_good, rx_vect; fs=res_impulse.fs, return_residual=true, window=window_impulsive)#, ref_channel=ref_channel, channels_relevant=1:size(rx_vect,2),
# getTDOA_func=get_tdoa_raw, solver_func=default_tdoa2dir_solver, cost_tdoa2ang=cost_tdoa2ang, return_residual=false)

scatter(res_impulse.pind_good_inS, angs[1].|>rad2deg; markershape=:xcross, alpha=res_impulse.ppeak.^2/maximum(res_impulse.ppeak.^2), labels=["azimuth" "inclination"]); 
title!(aufname|>basename; xlabel="Time(s)", ylabel="Angle(°)") |> display


# using Plots
# p1 = plot(angs[1] .|> rad2deg);
# p2 = plot(angs[2]);
# p3 = plot(res_impulse.ppeak);
# plot(p1, p2, p3, layout=(3,1)) 

# plot(angs .|> rad2deg)#, tdoas, seriestype=:scatter, xlabel="angle", ylabel="tdoa", legend=false)
# plot(tdoas[:,1:2])#, seriestype=:scatter, xlabel="angle", ylabel="tdoa", legend=false)

include("plotting.jl"); include("readImages.jl")
using VideoIO
outfol = joinpath(res_dir,splitext(aufname|>basename)[1])
plot_all_clicks(outfol, vidfname, res.pind_good_inS, res.data_filt, res.pind_good, tdoas, window_impulsive, angs[1]; plotsize=(1280,720))

outfol_img = joinpath(outfol, "%d.png")
outfol_clicks = outfol * "_clicks.mp4"

@ffmpeg_env run(`ffmpeg -framerate 30 -i $outfol_img -c:v libx264 -r 30 -pix_fmt yuv420p $outfol_clicks`)


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