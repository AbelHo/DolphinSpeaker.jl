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