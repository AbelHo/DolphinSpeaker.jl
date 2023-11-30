using DSP

function filter_simple(data, band_pass; fs=1, butterworth_size=butterworth_size)
    data_filt = data;
    if !iszero(band_pass[1]) || !isinf(band_pass[2])
        filter_type = nothing
        if iszero(band_pass[1]) 
            filter_type = Lowpass(band_pass[2]; fs=fs)
        elseif isinf(band_pass[2]) 
            filter_type = Highpass(band_pass[1]; fs=fs)
        else
            filter_type = Bandpass(band_pass[1], band_pass[2]; fs=fs)
        end

        filter_weight = digitalfilter(filter_type, Butterworth(butterworth_size))
        data_filt = mapslices( x -> filtfilt( filter_weight, x), data, dims=1)
    end
    return data_filt
end