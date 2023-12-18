using DSP
include("config.jl")
include("audacity.jl") # if required reading from audiofile directly instead of just calling detector


function detect_boat(aufname::String, res_dir=nothing; kwargs...)
    data, fs = readAudio(aufname)
    detect_boat((aufname,data,fs), res_dir; kwargs...)
end

function detect_boat(aufname_data_fs::Tuple, res_dir=nothing;
    ref_channel = ref_channel, thresh = threshold_boat,#7e-6,#2.5e-6,#threshold_boat, 
    nfft_inS = nfft_inS, band_pass = band_pass_boat, window=window_tonal)

    isinf(thresh) && return (;ppeak=nothing, time_index=nothing, tonal_indices=nothing, num_detection=-1, detected_tonal_inS=nothing, thresh=thresh, band_pass, nfft_inS, ref_channel)
    nfft_multiplier = 10

    aufname, data, fs = aufname_data_fs 
    filter_weight = digitalfilter(Bandpass(band_pass[1], band_pass[2], fs=fs), Butterworth(butterworth_size))
    data_filt = mapslices( x -> filtfilt( filter_weight, x), data, dims=1)

    lp = arraysplit(data_filt[:,ref_channel], Int(nfft_inS*fs*nfft_multiplier), 0)

    isempty(lp) && return (;num_detection=nothing, detected_tonal_inS=nothing, thresh, band_pass, nfft_inS, ref_channel)

    ppeak = energy.(lp)./length(lp[1])
    time_index = 0:nfft_inS*10:(nfft_inS*10*(length(lp)-1)) 

    tonal_indices = findall( >(thresh), ppeak)
    tonal_indices = tonal_indices[1:end-1] #FIXME #tempsolution

    num_detection = length(tonal_indices)
    @info "Num of Noise Detected: " * string(num_detection)
    detected_tonal_inS = time_index[tonal_indices]
    # mkpath(res_dir)
    # @debug "audacity label: "*joinpath(res_dir, splitext(aufname)[1]*"_tonal_t"*string(thresh_tonal)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".txt" |> basename)
    isnothing(res_dir) || audacity_label(detected_tonal_inS, joinpath(res_dir, splitext(aufname)[1]*"___noise_t"*string(thresh)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2]) *"_nfftS"*string(nfft_inS) *"_refCh"*string(ref_channel)* ".txt" |> basename) )
    

    return (;ppeak, time_index, tonal_indices, num_detection, detected_tonal_inS, thresh, band_pass, nfft_inS, ref_channel)
end