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
        # data_filt = mapslices( x -> filtfilt( filter_weight, x), data, dims=1)
        if size(data,2) > 1
            data_filt = mapslices2( x -> filtfilt( filter_weight, x), data)
        else
            data_filt = filtfilt( filter_weight, data)
        end
    end
    return data_filt
end

function extrema_index(arr)
    min_val, min_idx = arr[1], 1
    max_val, max_idx = arr[1], 1

    for i in 2:length(arr)
        if arr[i] < min_val
            min_val, min_idx = arr[i], i
        elseif arr[i] > max_val
            max_val, max_idx = arr[i], i
        end
    end

    return (min_val, max_val, min_idx, max_idx)
end


function mapslices2(func, arr::Array{T, 2} where T)
    # extrema_indices = Array{Tuple{T, Int, T, Int}, 1}(undef, size(arr, 2))
    extrema_indices = Array{Any}(undef,size(arr, 2))
    # lock = ReentrantLock()

    Threads.@threads for j in 1:size(arr, 2)
        # result = func(@view(arr[:,j]))
        # lock(extrema_indices) do
        #     extrema_indices[j] = result
        # end
        # min_val, min_idx = arr[1, j], 1
        # max_val, max_idx = arr[1, j], 1

        # for i in 2:size(arr, 1)
        #     if arr[i, j] < min_val
        #         min_val, min_idx = arr[i, j], i
        #     elseif arr[i, j] > max_val
        #         max_val, max_idx = arr[i, j], i
        #     end
        # end

        extrema_indices[j] = func(@view(arr[:,j])) #(min_val, min_idx, max_val, max_idx)
    end

    return hcat(extrema_indices...)
end


function extrema_and_indices(arr::Array{T, 2} where T)
    # extrema_indices = Array{Tuple{T, Int, T, Int}, 1}(undef, size(arr, 2))
    extrema_indices = Array{Tuple}(undef,size(arr, 2))

    Threads.@threads for j in 1:size(arr, 2)
        min_val, min_idx = arr[1, j], 1
        max_val, max_idx = arr[1, j], 1

        for i in 2:size(arr, 1)
            if arr[i, j] < min_val
                min_val, min_idx = arr[i, j], i
            elseif arr[i, j] > max_val
                max_val, max_idx = arr[i, j], i
            end
        end

        extrema_indices[j] = (min_val, max_val, min_idx, max_idx)
    end

    return extrema_indices
end