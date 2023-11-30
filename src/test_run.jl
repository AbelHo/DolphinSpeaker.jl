include("calibration_aspod4_2.jl")
include("config.jl")
using Pipe
using DelimitedFiles
using WAV
using ProgressMeter

# rx = 0.14722/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) )
# rx = 0.4/sqrt(3) .* exp.(im.* deg2rad.([30 -90 150]) )
# rx_vect = [real(rx); imag(rx); zeros(1,3)]

# aufname = "/Volumes/dd/Bahamas_2022/2022.06.27/0000/Aud_2022-06-27_11.16.17.wav"
# vidfname = "/Volumes/dd/Bahamas_2022/2022.06.27/0000/Vid_2022-06-27_11.16.16.mkv"
# res_dir = "/Users/abel/Documents/data_res/aspod/real/bahamas_2022"
#= USAGE:
pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect =
process_audioVideo(aufname, vidfname, "/Users/abel/Documents/data_res/aspod/real/bahamas_2022")
=#

# aufname = "/Users/abel/Documents/data/aspod/field/bahamas_2022/Aud_2022-06-25_10.33.28.wav"
# vidfname = "/Users/abel/Documents/data/aspod/field/bahamas_2022/Vid_2022-06-25_10.33.28.mkv"
function readAudio(aufname)
    filetype = splitext(aufname)[2] |> lowercase
    if ".wav" == filetype 
        data, fs = wavread(aufname, format="native")
    elseif ".mat" == filetype
        d = load(aufname)
        data = d["data"]
        haskey(d, "fs") ? fs=d["fs"] : fs=500_000

        if aufname[end-5:end-4] == "_1"
            if isfile( splitext(aufname)[1][1:end-1]*"2.mat" )
                @debug("Loading 2nd File: "*splitext(aufname)[1][1:end-1]*"2.mat")
                d = load(splitext(aufname)[1][1:end-1]*"2.mat")
                data = vcat(data, d["data"])
            end
        end
    else
        data, fs = load(aufname)
    end
    return data, Int(fs)
end


# aufname = "/Users/abel/Documents/aspod/data/2022-05-13_nus-pool/0015_wav/Aud_20131219_101319_20131219_101320.wav"
# vidfname = "/Users/abel/Documents/aspod/data/2022-05-13_nus-pool/0015/Vid_20131219_101319.mkv"
#"/Volumes/PortableSSD/aspod/aspod4_2/2022-09-11_BoatTest_wav/Aud_B001C0044_20220912040444_0001_20131219_021518.wav"
#"/Users/abel/Downloads/temp/data/aspod/maui_2022/2022.03.04/Aud_A003C0154_20220304130251_0001_20131219_100827.wav" #"/Users/abel/Downloads/temp/data/aspod/test/a4_2/2022-09-11_BoatTest_flac/Aud_B001C0044_20220912040444_0001_20131219_021518.flac"

function process_audioVideo(aufname::String, vidfname, res_dir; ref_channel=ref_channel, c=c, dist=dist_impulseive, window=window_impulsive, thresh=threshold_impulsive)
    data, fs = readAudio(aufname)
    process_audioVideo((aufname,data,fs), vidfname, res_dir; ref_channel=ref_channel, c=c, dist=dist_impulseive, window=window_impulsive, thresh=threshold_impulsive)
end

function process_audioVideo(aufname_data_fs::Tuple, vidfname, res_dir; ref_channel=ref_channel, c=c, dist=dist_impulsive, window=window_impulsive, thresh=threshold_impulsive)
    aufname, data, fs = aufname_data_fs
    # rx = 0.14722/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) )
    # rx_vect = [real(rx); imag(rx); zeros(1,3)]
    #=
    0.07361    5.2046e-18  -0.07361
    -0.0424988  0.0849975   -0.0424988
    0.0        0.0          0.0
    =#

    # ref_channel = 3
    # c = 1540 # m/s speed of sound
    # dist = 800# 15000#pinger
    # scatter(rx_vect[1,:],rx_vect[2,:],rx_vect[3,:])
    # xlabel!("X");ylabel!("Y")

    # vid = VideoIO.openvideo(vidfname)
    # imsize = raw_frame_size(vid)
    # data, fs = wavread(aufname, format="native")
    @info ("Duration: " * string(size(data,1)/fs) *"seconds")
    res_impulse = detect_impulse((aufname, data, fs), res_dir; return_datafilt=true)#; band_pass=[500 Inf], dist=80, threshold=.01)
    pind=res_impulse.pind_good; ppeak=res_impulse.ppeak; threshold_indices=res_impulse.threshold_indices; pind_good=res_impulse.pind_good; pind_good_inS=res_impulse.pind_good_inS;# data_filt=res_impulse.data_filt

    # pind, ppeak = findPings(data|>hilbert.|>abs; ref_channel=ref_channel, dist=dist)
    # # p = findPings(data|>hilbert.|>abs; ref_channel=3, dist=15000)

    # # window = -50:200; #-5000:10000 #-50:300
    # # thresh = .003#calf_hk .005*32767#aspod2 #0.01 # 0.1#pinger/clickler #0.2
    # threshold_indices = findall(>(thresh), ppeak)
    # pind_good = pind[threshold_indices]
    # pind_good_inS = (pind_good.-1) ./fs
    @info "Num of Clicks Detected: " * string(length(pind_good_inS))
    # overthresh=filter(x -> x>thresh, ppeak)

    # res_dir = "/Users/abel/Documents/data_res/aspod/real/bahamas_2022"
    isdir(res_dir) || mkpath(res_dir)
    audacity_label(pind_good./fs, joinpath(res_dir, splitext(aufname)[1]*"_t"*string(thresh)*"_d"*string(dist)*".txt" |> basename))

    tdoa = get_tdoa_envelope(data, pind_good; window=window, ref_channel=ref_channel)
    tdoa_raw = get_tdoa_raw(data, pind_good ; window=window, ref_channel=ref_channel)
    ang_env = tdoa2dir(tdoa, rx_vect,fs)
    ang = tdoa2dir(tdoa_raw, rx_vect,fs)

    p_pixels = angle2px(ang, fov_angle)
    pind_vidframes = round.(Int, pind_good_inS * get_fps(vidfname)) .+ 1

    # writedlm( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".csv"), ["p_pixel" "px" "py"], ',')
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes p_pixels], ',')
    end

    # run(`ffplay -f lavfi -i "sine=frequency=1000:duration=1" -autoexit -nodisp`)
    return pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect, fs#, data, fs
end

# laskdjfasdfasdf dsskdjfasdfasdf dsskdjfasdfasdf dsskdjfasdfasdf dsskdjfasdfasdf dsskdjfasdfasdf ds

# event_plots_dir = "$res_dir/bahamas_2022"

# event_plots_dir = joinpath(res_dir, basename(vidfname)*"_clicks"*"_"*"_t"*string(thresh)*"_d"*string(dist))
# plot_all_clicks(event_plots_dir*"_tdoa-raw", vidfname, pind_good_inS, data, pind[threshold_indices], tdoa_raw, window, ang)
# plot_all_clicks(joinpath(res_dir, splitext(basename(aufname))[1] ), vidfname, detected_tonal_inS, data_filt, (detected_tonal_inS.*fs).|>round.|>Int, hcat(tdoa...)', 1:1000, vcat(pxs...) .|> deg2rad)
function plot_all_clicks(event_plots_dir, vidfname, pind_good_inS, data, pind_threshold_indices, tdoa, window, ang; func2=plotTDOA_raw)
    vid = VideoIO.openvideo(vidfname)
    imsize = raw_frame_size(vid)
    mkpath(event_plots_dir)
    @time @showprogress "Writing each detection to image..." for i in eachindex(pind_good_inS)
        # @show angle2px(ang[i,:]') .|> round .|> Int |> string
        p = plotImg( readImage(vid, pind_good_inS[i]), ang, i)
        p = Plots.plot!(title=
            string(i)*"_"*string(pind_good_inS[i])*"s tdoa:"*string(tdoa[i,:]) *"\n"*
            string(@pipe ang[i,:] .|> rad2deg .|> round(_; digits=3))*"\n"*
            string(angle2px(ang[i,:]') .|> round .|> Int),
            legend=false)
        p = Plots.xlims!( 1, imsize[1])
        p = Plots.ylims!( 1, imsize[2])

        # savefig(p, joinpath(event_plots_dir, string(i)*".png"))

        p2 = func2(data, i, pind_threshold_indices, tdoa; window=window, func=gr)
        # p2 = plot!(title=
        #     string(i)*"_"*string(pind_good_inS[i])*"s tdoa:"*string(tdoa[i,:]),# *"\n"*
        #     # string(ang[i,:] .|> rad2deg),
        #     legend=false)

        pnew = Plots.plot(p, p2, layout=@layout [a b])
        savefig(pnew, joinpath(event_plots_dir, string(i)*".png"))
    end
end
# i=1; plotImg( readImage(vid, pind_good_inS[i]), ang, i)
# plot( pind_good/fs, (mod.(directions,2*pi) .|> rad2deg) .* sign.(directions) )


# event_signals_dir = "/Users/abel/Documents/data/aspod/field/bahamas_2022"
# event_signals_dir = joinpath(event_signals_dir, basename(vidfname)*"_signals")
# mkpath(event_signals_dir)
# for i in eachindex(pind_good_inS)
#     p = plotTDOA_raw(data, i, pind[threshold_indices], tdoa; window=window, func=gr)
#     p = plot!(title=
#         string(i)*"_"*string(pind_good_inS[i])*"s tdoa:"*string(tdoa[i,:]),# *"\n"*
#         # string(ang[i,:] .|> rad2deg),
#         legend=false)
#     savefig(p, joinpath(event_signals_dir, string(i)*".png"))
# end

# run(`ffplay -f lavfi -i "sine=frequency=1000:duration=1" -autoexit -nodisp`)

# return pind_vidframes, p_pixels