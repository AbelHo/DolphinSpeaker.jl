using FFTW
using DSP, SignalAnalysis
using JLD2
include("config.jl")

# aufname = "/Volumes/dd/Bahamas_2022/2022.06.27/0000/Aud_2022-06-27_11.16.17.wav"
# vidfname = "/Volumes/dd/Bahamas_2022/2022.06.27/0000/Vid_2022-06-27_11.16.16.mkv"
#res_dir
#fs
#data, assumise last channel is synch channel and not used for processing

function process_audioVideo_tonal1(aufname::String, vidfname, res_dir)
    filetype = splitext(aufname)[2] |> lowercase
    if ".wav" == filetype 
        data, fs = wavread(aufname, format="native")
    elseif ".mat" == filetype
        d = load(aufname)
        data = d["data"]
        haskey(d, "fs") ? fs=d["fs"] : fs=500_000
    else
        data, fs = load(aufname)
    end

    process_audioVideo_tonal1( (aufname,data,fs), vidfname, res_dir)
end

#     process_audioVideo_tonal1(, vidfname, res_dir)

function process_audioVideo_tonal1(aufname_data_fs::Tuple, vidfname, res_dir;
    ref_channel = ref_channel, thresh_tonal = threshold_tonal, 
    nfft_inS = nfft_inS, band_pass = band_pass, 
    butterworth_size = butterworth_size,c = c, winlen=winlen)

    aufname, data, fs = aufname_data_fs
    # ref_channel = 3
    # thresh_tonal = -20
    # nfft_inS = 0.01
    # band_pass = [5000, 14000] # [5000, 25000] # 
    # butterworth_size = 4
    # c = 1540

    # winlen=0.01s # tdoa estimation window length

    # data, fs = wavread(aufname, format="native")
    @info ("Duration: " * string(size(data,1)/fs) *"seconds")

    freq_filt = Int(round(band_pass[1]/100)):Int(round(band_pass[2]/100)) #50:140
    # FIXME automatically find according to actual frequency
    filter_weight = digitalfilter(Bandpass(band_pass[1],band_pass[2], fs=fs), Butterworth(butterworth_size)) #Bandpass(5000,14000, fs=fs)

    nfft = nextfastfft(nfft_inS*fs)
    specs = mapslices( x -> spectrogram(x, nfft, 0; fs=fs), data[:,1:3], dims=1)
    mag_ft = broadcast( x -> pow2db.(x.power) ,specs)

    time_index = first(specs).time#[window]
    fft_max = maximum(mag_ft[ref_channel][freq_filt,:]; dims=1)'[:]

    tonal_indices = findall( >(thresh_tonal), fft_max)
    tonal_indices = tonal_indices[1:end-1] #FIXME #tempsolution

    @info "Num of Whistles Detected: " * string(length(tonal_indices))
    detected_tonal_inS = time_index[tonal_indices]
    mkpath(res_dir)
    audacity_label(detected_tonal_inS, joinpath(res_dir, splitext(aufname)[1]*"_tonal_t"*string(thresh_tonal)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".txt" |> basename) )

    data_filt = mapslices( x -> filtfilt( filter_weight, x), data[:,1:end-1], dims=1) #data, assumise last channel is synch channel and not used for processing
    d_filt=signal(data_filt, fs)

    tdoa = twindow2tdoa.( (detected_tonal_inS)s; fs=fs, d_filt=d_filt, winlen=winlen)
    pxs = twindow2dird.( (detected_tonal_inS)s; fs=fs, d_filt=d_filt, winlen=winlen)

    ## create return variables:
    ang = vcat(pxs...) #.|> deg2rad
    p_pixels = angle2px(ang)

    pind_vidframes = round.(Int, detected_tonal_inS * get_fps(vidfname)) .+ 1    

    # writedlm( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh_tonal)*"tonal_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".csv"), ["p_pixel" "px" "py"], ',')
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh_tonal)*"tonal_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes p_pixels], ',')
    end
    
    return pind_vidframes, p_pixels, thresh_tonal, band_pass, vcat(pxs...) .|> deg2rad, tdoa, nothing, winlen, tonal_indices, nothing, detected_tonal_inS, nothing, fft_max[tonal_indices], ref_channel, c, rx_vect, fs
end



# tdoa = twindow2tdoa.( (detected_tonal_inS)s; d_filt=d_filt, winlen=winlen)
# plot_all_clicks(joinpath(res_dir, splitext(basename(aufname))[1] ), vidfname, detected_tonal_inS, data_filt, (detected_tonal_inS.*fs).|>round.|>Int, hcat(tdoa...)', 1:500, vcat(pxs...) .|> deg2rad)

function twindow2tdoa(win1; fs=fs, d_filt=d_filt, winlen=winlen, ref_channel=3, rx_vect=rx_vect)
    window = win1:(win1+winlen);
    tdoa = [finddelay(d_filt[window,1],d_filt[window,ref_channel]); finddelay(d_filt[window,2],d_filt[window,ref_channel]); finddelay(d_filt[window,3],d_filt[window,ref_channel])];
end

function twindow2dird(win1; fs=fs, d_filt=d_filt, winlen=winlen, ref_channel=3, rx_vect=rx_vect)
    window = win1:(win1+winlen);
    tdoa = [finddelay(d_filt[window,1],d_filt[window,ref_channel]); finddelay(d_filt[window,2],d_filt[window,ref_channel]); finddelay(d_filt[window,3],d_filt[window,ref_channel])];
    # tdoa_ori = [finddelay(d[window,1],d[window,3]); finddelay(d[window,2],d[window,3]); finddelay(d[window,3],d[window,3])];
    # tdoas = [tdoa tdoa_ori]
    # tdoa2dir(tdoas', rx_vect) .|> rad2deg
    tdoa2dir(tdoa', rx_vect, fs) #.|> rad2deg
    # tdoa2dir(tdoa_ori', rx_vect) .|> rad2deg
end