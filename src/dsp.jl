using DSP
import SignalAnalysis.findsignal

function filter_simple(data, band_pass, band_stop=nothing; fs=1, butterworth_size=butterworth_size)
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
        if !isnothing(band_stop)
            if eltype(band_stop) <: Number
                @debug band_stop
                filter_type = Bandstop(band_stop[1], band_stop[2]; fs=fs)
                filter_weight = filter_weight * digitalfilter(filter_type, Butterworth(butterworth_size))
            else
                for bs in band_stop
                    @debug bs
                    filter_type = Bandstop(bs[1], bs[2]; fs=fs)
                    filter_weight = filter_weight * digitalfilter(filter_type, Butterworth(butterworth_size))
                end
            end

        end
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

function extrema_and_indices(arr::Array{T, 1} where T)
    # extrema_indices = Array{Tuple{T, Int, T, Int}, 1}(undef, size(arr, 2))
    extrema_indices = Array{Tuple}(undef,size(arr, 1))

    # Threads.@threads for j in 1:size(arr, 2)
        min_val, min_idx = arr[1], 1
        max_val, max_idx = arr[1], 1

        for i in 2:size(arr, 1)
            if arr[i] < min_val
                min_val, min_idx = arr[i], i
            elseif arr[i] > max_val
                max_val, max_idx = arr[i], i
            end
        end

        extrema_indices = (min_val, max_val, min_idx, max_idx)
    # end

    return extrema_indices
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

# """
# 	funcOnWindows(data, windows = default_window ; func=x->x, kwargs...)

# Apply a function `func` on windows of data.

# # Arguments
# - `data`: The input data matrix.
# - `windows`: An array of window indices.
# - `func`: The function to apply on each window. Default is the identity function.
# - `kwargs`: Additional keyword arguments to be passed to `func`.

# # Returns
# - `out`: An array of the same size as `windows` and `data` containing the results of applying `func` on each window.

# # Examples:
# ```julia-repl
# julia> out = funcOnWindows(res.res_impulse.data_filt, map(x-> x.+ (-300:300), res.res_impulse.pind_good); func=maximum, dims=1)
# ```

# """
function funcOnWindows(data, windows = default_window, args... ; func=x->x, kwargs...)
	out = Array{Any}(undef,length(windows),size(data,2))
	Threads.@threads for i in eachindex(windows)
		window = windows[i]
		for ch in 1:size(data,2)
			if ch==ref_channel
				out[i,ch]=0
				continue
			end
			out[i,:] = func(data[window, :], args...; kwargs...)
		end
	end
	return out
end


import Optim: optimize, minimizer, BFGS

"""
Finds up to `n` copies of reference signal `r` in signal `s`. The reference
signal `r` should have a delta-like autocorrelation for this function to work
well. If the keyword parameter `coarse` is set to `true`, approximate arrival
times are computed based on a matched filter. If it is set to `false`, an
iterative optimization is performed to find more accruate arrival times.

Returns named tuple `(time=t, amplitude=a)` where `t` is a vector of arrival
times and `a` is a vector of complex amplitudes of the arrivals. The arrivals
are sorted in ascending order of arrival times.

# Examples:
```julia-repl
julia> x = chirp(1000, 5000, 0.1, 40960; window=(tukey, 0.05))
julia> x4 = resample(x, 4)
julia> y4 = samerateas(x4, zeros(32768))
julia> y4[128:127+length(x4)] = real(x4)          # time 0.000775𝓈, index 32.75
julia> y4[254:253+length(x4)] += -0.8 * real(x4)  # time 0.001544𝓈, index 64.25
julia> y4[513:512+length(x4)] += 0.6 * real(x4)   # time 0.003125𝓈, index 129.0
julia> y = resample(y4, 1//4)
julia> y .+= 0.1 * randn(length(y))
julia> findsignal(x, y, 3; coarse=true)
(time = Float32[0.000781, 0.001538, 0.003125], amplitude = ComplexF64[...])
julia> findsignal(x, y, 3)
(time = Float32[33, 64, 129], [0.000775, 0.001545, 0.003124], amplitude = ComplexF64[...])
```
"""
# function findsignal2(r, s, n=1; prominence=0.2, coarse=false)
#   # coarse arrival time estimation
#   r = analytic(r)
#   r = r / std(r)
#   s = analytic(s)
#   mfo = mfilter(r, s) / length(r)
#   absmfo = abs.(samples(mfo))
#   p, _ = findmaxima(absmfo)
#   peakproms!(p, absmfo; minprom=prominence*maximum(absmfo))
#   length(p) > length(s)/10 && return (time=Float64[], amplitude=ComplexF64[])
#   h = absmfo[p]
#   ndx = sortperm(h; rev=true)
#   length(ndx) > n && (ndx = ndx[1:n])
#   p = p[ndx]
#   if coarse
#     t = time(Float64.(p), s)
#     ndx = sortperm(t)
#     return (time=t[ndx], amplitude=samples(mfo[p[ndx]]))
#   end
#   # iterative fine arrival time estimation
#   margin = 5   # arrival time may vary up to margin from coarse estimates
#   i::Int = minimum(p)
#   n = maximum(p) - i + length(r) + 2 * margin
#   n = nextfastfft(n)
#   i = max(1, i - margin)
#   N = n
#   i + N - 1 > length(s) && (N = length(s) - i + 1)
#   X = fft(vcat(samples(r), zeros(n-length(r))))
#   soln = let p=p, f=fftfreq(n)
#     function reconstruct(v)
#       ii = @view v[1:length(p)]
#       aa = @views complex.(v[length(p)+1:2*length(p)], v[2*length(p)+1:3*length(p)])
#       Z = mapreduce(+, zip(ii, aa)) do (i, a)
#         a .* X .* cis.(-2π .* i .* f)
#       end
#       @view real(ifft(Z))[1:N]
#     end
#     v0 = [p .- i; real.(mfo[p]); imag.(mfo[p])]
#     optimize(v -> sum(abs2, reconstruct(v) .- s[i:i+N-1]), v0)
#   end
#   v = minimizer(soln)
#   pp = v[1:length(p)] .+ i

# #   try
# #   t = time(pp, s)
# #   a = complex.(v[length(p)+1:2*length(p)], v[2*length(p)+1:3*length(p)])
# #   ndx = sortperm(t)
# #   (time=t[ndx], amplitude=a[ndx], mfo=mfo ? m : empty(m))
# #   catch err
# #     # println(err)
# #   end
#   t = pp
#   a = complex.(v[length(p)+1:2*length(p)], v[2*length(p)+1:3*length(p)])
#   ndx = sortperm(t)
#   (sample=t[ndx], amplitude=a[ndx])
# end

# function findsignal(r, s, n=1; prominence=0.0, finetune=2, mingap=1, mfo=false)
#     # coarse arrival time estimation
#     r = analytic(r)
#     s = analytic(s)
#     m = mfilter(r, s)
#     m ./= (std(r) * length(r))
#     T = eltype(m)
#     m̄ = abs.(samples(m))
#     p = argmaxima(m̄, mingap)
#     prominence > 0 && peakproms!(p, m̄; minprom=prominence*maximum(m̄))
#     if length(p) > length(s)/10
#       return (time=Float64[], amplitude=T[], mfo=mfo ? m : empty(m))
#     end
#     h = m̄[p]
#     ndx = partialsortperm(h, 1:n; rev=true)
#     p = p[ndx]
#     if finetune == 0
#       t = time(Float64.(p), s)
#       ndx = sortperm(t)
#       return (time=t[ndx], amplitude=samples(m[p[ndx]]), mfo=mfo ? m : empty(m))
#     end
#     # iterative fine arrival time estimation
#     i::Int = minimum(p)
#     N::Int = maximum(p) - i + length(r) + 2 * finetune
#     N = nextfastfft(N)
#     i = max(1, i - finetune)
#     X = fft(vcat(samples(r), zeros(N-length(r))))
#     S = fft(vcat(samples(s[i:min(i+N-1,end)]), zeros(max(i+N-1-length(s),0))))
#     soln = let P=length(p), f=fftfreq(N)
#       optimize([p .- i; real.(m[p]); imag.(m[p])], BFGS(); autodiff=:forward) do v
#         ii = @view v[1:P]
#         aa = @views complex.(v[P+1:2P], v[2P+1:3P])
#         X̄ = mapreduce(+, zip(ii, aa)) do (i, a)
#           a .* X .* cis.(-2π .* i .* f)
#         end
#         sum(abs2, X̄ .- S)
#       end
#     end
#     v = minimizer(soln)
#     pp = v[1:length(p)] .+ i
#     t = time(pp, s)
#     a = complex.(v[length(p)+1:2*length(p)], v[2*length(p)+1:3*length(p)])
#     ndx = sortperm(t)
#     (time=t[ndx], amplitude=a[ndx], mfo=mfo ? m : empty(m))
#   end

#~ find extrema
extrema_in_file(aufname; res_dir=nothing) = extrema_in_file(aufname, res_dir)
function extrema_in_file(aufname, res_dir=nothing)
    if isdir(aufname)
        return extrema_in_file.( readdir.(aufname; join=true) |> skiphiddenfiles, Ref(res_dir) )
    end

    try
        if length(aufname) < 4 
            return
        elseif aufname[end-3:end] == "flac"
            data, fs = flac2signal(aufname)
        elseif aufname[end-2:end] == "mat"
            data, fs, _,_, timestamp = readAudio(aufname)
        else
            return
        end
        @info aufname
        mkpath(res_dir)

        # data, fs, _,_, timestamp = readAudio(aufname)
        # data, fs = flac2signal(aufname)
        res = (extrema_and_indices(data), energy(data), aufname)
        res_dict = Dict("extrema_and_indices" => res[1], "energy" => res[2], "aufname" => aufname)
        if res_dir != nothing
            open(joinpath(res_dir, basename(aufname)*".json"), "w") do f
                JSON.print(f,res_dict,4)
                # write(f, res_dict)
            end
            save(joinpath(res_dir, basename(aufname)*".jld2"), "res", res)
        else
            JSON.print(stdout,res_dict,4)
        end

    catch err
        @error(aufname)
        println(err)
        return
    end
end


function previousfastfft(n, lessby=0)
    if nextfastfft(n-lessby) == nextfastfft(n)
        return previousfastfft(n, lessby+1)
    else
        return nextfastfft(n-lessby)
    end
end

function truncate_fft_end(data)
    newdata = @view(data[begin:previousfastfft(size(data,1)), :])
    size(newdata,1) < size(data,1) && @info("truncated data for FFT")
    newdata
end

windowing(data, i, window, col=1:size(data,2)) = data[window[begin]+i : window[end]+i, col]

function add(a,b)
    a_or_b = size(a,1) < size(b,1)
    len = a_or_b ? size(a,1) : size(b,1)
    a_or_b ? a .+ b[1:len,:] : a[1:len,:] .+ b
end

extrema_abs(args; dims=1, kwargs...) = extrema(args; dims=dims, kwargs...) .|> x->maximum(abs.(x))

p2p_db(x) = 20 * log10( -reduce(-, extrema(x)) )
p2p_db(arr::Array{T, 2} where T) = mapslices2(p2p_db, arr)

norm_max(args...; norm_func=x->maximum(abs.(x); dims=1), kwargs...) = args[1]./norm_func(args[1])


"""
    finddelay2(x, y)

Estimate the delay of x with respect to y by locating the peak of their
cross-correlation.

The output delay will be positive when x is delayed with respect y, negative if
advanced, 0 otherwise.

# Example
```jldoctest
julia> finddelay2([0, 0, 1, 2, 3], [1, 2, 3])
2

julia> finddelay2([1, 2, 3], [0, 0, 1, 2, 3])
-2
```
"""
function finddelay2(x_o::AbstractVector{<: Real}, y_o::AbstractVector{<: Real}; norm_func=x->x)
    x = norm_func(x_o)
    y = norm_func(y_o)

    s = xcorr(y, x, padmode=:none)
    max_corr = maximum(abs, s)
    max_idxs = findall(x -> abs(x) == max_corr, s)

    center_idx = length(x)
    # Delay is position of peak cross-correlation relative to center.
    # If the maximum cross-correlation is not unique, use the position
    # closest to the center.
    d_ind = argmin(abs.(center_idx .- max_idxs))
    d = center_idx - max_idxs[d_ind]
    return d, max_corr, s
end

# a = [
#     begin
#         snip = windowing(data_filt, d["res_impulsetrain"].pind_good[i], window_impulsive)
#         norm_snip = norm_max(snip)
#         extrema_abs([mfilter(norm_snip[:,1], norm_snip[:,3]) mfilter(norm_snip[:,1], norm_snip[:,2]) mfilter(norm_snip[:,3], norm_snip[:,2])]|>collect)
#     end
    
#     for i in 1:length(d["res_impulsetrain"].pind_good_inS)
# ]