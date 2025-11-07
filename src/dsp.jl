using DSP
import SignalAnalysis.findsignal

function filter_simple(data, band_pass, band_stop=nothing; fs=1, butterworth_size=butterworth_size, mapslices2=mapslices2, kwargs...)
    data_filt = data;
    if !iszero(band_pass[1]) || !isinf(band_pass[2])
        filter_type = nothing
        if iszero(band_pass[1]) 
            filter_type = Lowpass(band_pass[2]/fs*2)
        elseif isinf(band_pass[2]) 
            filter_type = Highpass(band_pass[1]/fs*2)
        else
            filter_type = Bandpass(band_pass[1]/fs*2, band_pass[2]/fs*2)
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
                    filter_type = Bandstop(bs[1]/fs*2, bs[2]/fs*2)
                    filter_weight = filter_weight * digitalfilter(filter_type, Butterworth(butterworth_size))
                end
            end

        end
        # data_filt = mapslices( x -> filtfilt( filter_weight, x), data, dims=1)
        if size(data,2) > 1
            data_filt = mapslices2( x -> filtfilt( filter_weight, x), data; kwargs...)
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


function mapslices2(func, arr::AbstractArray{T, 2} where T; kwargs...)
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

        extrema_indices[j] = func(@view(arr[:,j]); kwargs...) #(min_val, min_idx, max_val, max_idx)
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

function mapblocks(func, arr::AbstractArray{T} where T; dims=ndims(arr), kwargs...)
    # out_arr = Array{Any}(undef,size(arr, dims))
    # for (i, block) in enumerate(eachslice(arr; dims=dims))
    #     out_arr[i] = func(block)
    # end
    out_arr = Array{typeof(first(arr))}(undef, size(func(selectdim(arr, dims, 1)))..., size(arr, dims))
    outdim = ndims(out_arr)
    Threads.@threads for i in 1:size(out_arr, outdim)
        @inbounds selectdim(out_arr, outdim, i) .= func(selectdim(arr, dims, i))
    end
    # cat(out_arr, dims=dims)
    return out_arr
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

# peak to peak dB calculation
p2p_db(x) = 20 * log10( -reduce(-, extrema(x)) )
p2p_db(arr::Array{T, 2} where T) = mapslices2(p2p_db, arr)

# example normalization function
# nfunc(x) = sqrt(sum(abs2.(x)))
# nfunc(x) = energy(x)
# nfunc(x) = 1
# nfunc(x) = maximum(abs.(x); dims=1)
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
function finddelay2(x_o::AbstractVector{<: Real}, y_o::AbstractVector{<: Real};
    norm_func=x->x, flag_norm_rms=false)

    x = norm_func(x_o)
    y = norm_func(y_o)

    s = xcorr(y, x, padmode=:none)
    # @debug argmax(abs.(s))
    max_corr = maximum(abs, s)
    # @debug max_corr
    max_idxs = findall(x -> abs(x) == max_corr, s)

    center_idx = length(x)
    # Delay is position of peak cross-correlation relative to center.
    # If the maximum cross-correlation is not unique, use the position
    # closest to the center.
    d_ind = argmin(abs.(center_idx .- max_idxs))
    d = center_idx - max_idxs[d_ind]

    max_corr = s[max_idxs[d_ind]]

    if flag_norm_rms
        norm_factor = sqrt(sum(abs2, x) * sum(abs2, y))
        max_corr /= norm_factor
        s ./= norm_factor
    end
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

# using DSP
using SignalBase
function psd2(data; fs=1.0, nfft=512, noverlap=div(nfft,2),
    window=hamming(nfft), xscale=:auto, yrange=50)
    p=[];pow=Array{Float64}(undef, Int(nfft/2+1), size(data,2));
    for i = 1:size(data,2)
        p = welch_pgram(data[:,i], nfft, noverlap; fs=inHz(fs), window=window)
        pow[:,i] = 10*log10.(p.power)
    end

    return pow, p.freq
end

function correct_ambient(data_fs; ch_list=nothing, correction=nothing, kwargs...)
    @info ch_list, correction
    if data_fs isa Tuple
        data,fs = data_fs
    end
    if isnothing(ch_list) && isnothing(correction)
        ch_list, correction, freqss, pow, pp, ch_noisylist, data = ambientnoise_correction(data_fs; kwargs...)
    end
    if !isempty(ch_list)
        data_corrected = Float64.(data)
        data_corrected[:,ch_list] = data_corrected[:,ch_list] .* 10 .^ (correction./20)
        return data_corrected
    end
    return data
end

function ambientnoise_correction(data_fs; res_dir=nothing, ch_list=:auto, ch_db=:auto, threshold_fft_error=:auto, threshold_fft_error_multiple=.35, rx_vect=nothing)
	# data, fs = readAudio(aufname)
	data, fs = data_fs
	!isnothing(rx_vect) && (data = @view data[:,get_relevant_channels(rx_vect)])
	# @info size(data)
	pow, freqss = psd2(data; fs=fs)
	correction=nothing

	pp=copy(pow)
	pow_median = nothing
	if ch_list == :auto
		pow_median = mapslices(median,pow;dims=2)
		pow_diff = pow_median .- pow
		threshold_fft_error==:auto && (threshold_fft_error = size(pow,1) * threshold_fft_error_multiple; @info "threshold_fft_error: $threshold_fft_error")
		count_result = count(pow_diff .> 3; dims=1)[:]
		ch_list = findall(>(threshold_fft_error), count_result)
		@info "corrected channels: $ch_list"
		@debug count_result

		# plot(pow_diff)|>display
		count_noisyresult = count(pow_diff .< -3; dims=1)[:]
		ch_noisylist = findall(>(threshold_fft_error), count_noisyresult)
		@info "noisy list: $ch_noisylist"
	else
		ch_list = nothing
	end

	if !isnothing(ch_list)
		if ch_db != :auto
			pp[:, ch_list] = pow[:, ch_list] .+ ch_db
			correction = ch_db
		else
			correction = median(pow_median .- pow[:,ch_list]; dims=1)
			# @info size(pow[:, ch_list]), size(correction)
		end
		pp[:, ch_list] = pow[:, ch_list] .+ correction
		@info "corrections: $correction"
	end
	return ch_list, correction, freqss, pow, pp, ch_noisylist, data
end
ambientnoise_correction(aufname::String; kwargs...) = ambientnoise_correction(readAudio(aufname); kwargs...)

"""
Example:
```
 ch_list, correction, freqss, pow, pp, ch_noisylist = find_correction(in_dir; res_dir=res_dir, rx_vect=rx_vect, ch_db=9)
```
"""
function find_correction(in_dir; func_filter= y-> joinpath(y, "acoustic", filter(x->startswith(x,"Ambient") && endswith(x,".ogg"), readdir(joinpath(y,"acoustic")) )[1]) , kwargs...)
    aufname = func_filter(in_dir)
    @info aufname
    output = ambientnoise_correction(aufname; kwargs...)
    return output
end

    

filter_band(x, band_pass=[0, Inf]; fs=fs) = mapslices(extrema, filter_simple(x, band_pass; fs=fs, mapslices2=mapslices, dims=1);dims=1)
# filter_bandfft()


diff_func(a::AbstractVector; kwargs...) = diff_func(a; dims=1, kwargs...)

"""
    diff_func(A::AbstractVector)
    diff_func(A::AbstractArray; dims::Integer)

Finite diff_funcerence operator on a vector or a multidimensional array `A`. In the
latter case the dimension to operate on needs to be specified with the `dims`
keyword argument.

!!! compat "Julia 1.1"
    `diff_func` for arrays with dimension higher than 2 requires at least Julia 1.1.

# Examples
```jldoctest
julia> a = [2 4; 6 16]
2×2 Matrix{Int64}:
 2   4
 6  16

julia> diff_func(a, dims=2)
2×1 Matrix{Int64}:
  2
 10

julia> diff_func(vec(a))
3-element Vector{Int64}:
  4
 -2
 12
```
"""
function diff_func(a::AbstractArray{T,N}; dims::Integer, kwargs...) where {T,N}
    require_one_based_indexing(a)
    1 <= dims <= N || throw(ArgumentError("dimension $dims out of range (1:$N)"))

    r = axes(a)
    r0 = ntuple(i -> i == dims ? UnitRange(1, last(r[i]) - 1) : UnitRange(r[i]), N)
    r1 = ntuple(i -> i == dims ? UnitRange(2, last(r[i])) : UnitRange(r[i]), N)

    return view(a, r1...) .- view(a, r0...)
end
function diff_func(r::AbstractRange{T}; dims::Integer=1, func=(r,i)->r[i+1]-r[i]) where {T}
    dims == 1 || throw(ArgumentError("dimension $dims out of range (1:1)"))
    return [@inbounds func(r,i) for i in firstindex(r):lastindex(r)-1]
end


"""
    similarity_matrix(signals::Vector{Vector{Float64}})

Computes the normalized cross-correlation similarity score between all pairs of signals.
Returns a matrix where element (i, j) is the maximum normalized cross-correlation between signals[i] and signals[j].
"""
function similarity_matrix(signals::Vector{Vector{Float64}})
    n = length(signals)
    sim_matrix = zeros(Float64, n, n)

    for i in 1:n
        for j in i:n
            sig1 = signals[i]
            sig2 = signals[j]
            # Normalize both signals to zero mean and unit variance
            sig1_norm = (sig1 .- mean(sig1)) ./ std(sig1)
            sig2_norm = (sig2 .- mean(sig2)) ./ std(sig2)

            # Compute full cross-correlation
            cc = xcorr(sig1_norm, sig2_norm)#, mode = :full)
            # Normalize by product of norms to get correlation coefficient
            norm_factor = sqrt(sum(sig1_norm .^ 2) * sum(sig2_norm .^ 2))
            ncc = cc ./ norm_factor

            # Max correlation is our similarity score
            max_sim = maximum(ncc)
            sim_matrix[i, j] = max_sim
            sim_matrix[j, i] = max_sim  # symmetric
        end
    end

    return sim_matrix
end

# # Example usage:
# s1 = randn(100)
# s2 = 2 .* s1 .+ 0.1 .* randn(100)  # Highly similar, different scale
# s3 = reverse(s1)                   # Less similar
# s4 = randn(80)                     # Random, shorter length

# signals = [s1, s2, s3, s4]
# sim_matrix = similarity_matrix(signals)
# println("Similarity Matrix:\n", round.(sim_matrix; digits=3))

# compute normal rfft
function compute_rfft(snip, fs=1.0; type=:amplitude, plot=plot) 
	fft_val = rfft(snip, 1) .|> abs
	freqss =  fftfreq2(size(snip,1),fs)  #0:(fs/size(snip,1)):fs÷2

	type == :log && (fft_val = 20 .* log10.(fft_val))
	@debug (size(snip), size(freqss), size(fft_val), typeof(freqss), typeof(fft_val))
	fft_val, freqss
end

"""
sig2rgb(sig; fs=1.0, rgb_bands=[[1000, 70_000], [70_000, 120_000], [120_000, 170_000]], kwargs...)

Convert a time-domain signal into an RGBA representation by aggregating FFT magnitudes
over three user-defined frequency bands (R, G, B) and deriving an alpha channel from
the aggregated energy.

Arguments
- sig: Input signal. Typically a 1-D array of samples, but any shape accepted by
    compute_rfft is allowed (compute_rfft must return magnitude data with frequency
    bins on the first dimension and frames/columns on the second).
- fs::Real: Sampling frequency (Hz). Passed to compute_rfft. Default: 1.0.
- rgb_bands::AbstractVector{<:AbstractVector}: A length-3 collection of 2-element
    ranges [low, high] (same units as fs) that define the frequency bands mapped to
    the R, G and B channels respectively. Default is [[1000, 70000], [70000, 120000], [120000, 170000]].
- kwargs...: Additional keyword arguments forwarded to compute_rfft (e.g., window, nfft).

Returns
- (rgb_array, rgba_view)
    - rgb_array::Array{Float32,2}: A 4×Nframes array where rows 1..3 are the R,G,B channel
        magnitudes (per-frame normalized to [0,1]) and row 4 is the alpha channel computed
        as the per-frame total energy normalized by the maximum total energy across frames.
    - rgba_view: A ColorTypes-compatible view created with colorview(RGBA, rgb_array)
        suitable for visualization or image I/O.

Behavior
- The function calls compute_rfft(sig, fs; kwargs...) and expects two return values:
    (fft_val, freqss). fft_val should be an array of non-negative magnitudes with
    freq bins along the first dimension and frames along the second; freqss should
    be a vector of frequency values corresponding to the rows of fft_val.
- For each band in rgb_bands, the function selects frequency bins where freqss ∈ [low, high)
    and computes the mean magnitude across those bins for each frame to produce each color channel.
- The alpha channel is computed as the framewise sum of R+G+B, then normalized by the
    maximum sum across all frames to lie in [0,1].
- Each color channel R,G,B is normalized per-frame by the maximum among the three channels
    for that frame to scale values to [0,1].

Notes / Caveats
- If a specified band contains no matching frequency bins (e.g., band outside the FFT range),
    the mean over an empty slice will produce NaNs. Ensure rgb_bands intersect freqss.
- If the per-frame maximum across R/G/B is zero (silent frame), per-frame normalization
    may produce NaN or Inf; callers may wish to filter or clamp such frames.
- Units of rgb_bands must match units of freqss (typically Hz).
- compute_rfft must produce magnitude values (not complex spectra); if it returns complex
    results the user should ensure magnitudes are passed (or compute_rfft should do so).

Example
- Typical call:
    sig2rgb(signal, fs=48000, rgb_bands=[[1000,7000],[7000,12000],[12000,20000]])
    returns a 4×N array and a corresponding color view for display.
"""
function sig2rgb(sig; fs=1.0, rgb_bands=[[1000, 70_000], [70_000, 120_000], [120_000, 170_000]], kwargs...)
    fft_val, freqss = compute_rfft(sig, fs; kwargs...)

    # Map frequency bands to RGB channels
    rgb = zeros(Float32, 4, size(fft_val,2))
    for (i, band) in enumerate(rgb_bands)
        # Find frequencies within the band
        mask = (freqss .>= band[1]) .& (freqss .< band[2])
        # Compute average magnitude in the band
        rgb[i, :] = mean(fft_val[mask,:], dims=1)
    end

    sum_rgb = sum(rgb; dims=1)
    rgb[4,:] = sum_rgb./maximum(sum_rgb)

    # Normalize each channel to [0, 1]
    rgb[1:3, :] .= @view(rgb[1:3, :]) ./ maximum(@view(rgb[1:3, :]); dims=1)


    # convert to color
    rgb, colorview(RGBA, rgb)

    # return rgb
end

function get_color_clicks(win_anal, savefname;
	clips=clips, res=res, fs=fs, 
	rgb_bands=[[10_000, 60_000], [60_000, 110_000], [110_000, 160_000]],
	rgbs_alpha_offset=0.0)
	if dirname(savefname) |> isdir == false
		mkpath(dirname(savefname))
	end
	rgbs, rgba = sig2rgb(clips[:,win_anal]; fs=fs, rgb_bands=rgb_bands)
	rgbs[4,:] .+= rgbs_alpha_offset
	p = plot(res.res_impulsetrain.pind_good_inS[win_anal], rgbs[4,:]; 
		color=rgbs|> eachcol .|> x-> RGBA(x...), 
		# color=rgba,
		seriestype=:scatter, 
		hover = string.(win_anal) .* ", " .* string.(round.(res.res_impulsetrain.pind_good_inS[win_anal]; digits=3)) .*"s",
		bg=:black,markerstrokewidth = 0,
		size=(1000,600))
	savefig(savefname)
	return p
end


function extract_clips_single_channel(clips_fixed, ref_channel=0)
    if ref_channel != 0
        return clips = hcat(map(x-> x[:,ref_channel], clips_fixed)...)
    else
        return clips = hcat(map(x-> x[:, argmax(maximum(abs.(x); dims=1))[2] ], clips_fixed)...)
    end
end

# function filter_extract(data::AbstractMatrix, fs::Real; 
#         band_pass=[0, Inf], band_stop=nothing, ref_channel=0, window_impulsive=(-100:100), pind_good=Int[], kwargs...)
#     data_filt = filter_simple(data, band_pass, band_stop; fs=fs)

    
#         data_filt = filter_simple(data, band_pass, band_stop; fs=fs, mapslices2=mapslices2, dims=1, kwargs...)
#     clips_fixed = extract_clips(data_filt, window_extract .+ res.res_impulsetrain.pind_good, NaN; flag_matrix=false)
#     clips = extract_clips_single_channel(clips_fixed, ref_channel)
#     return clips, data_filt
# end

# function filter_extract(data::AbstractMatrix, fs::Real; band_pass=[0, Inf], band_stop=nothing, ref_channel=0, window_impulsive=(-100:100), pind_good=Int[], kwargs...)
#     data_filt = filter_simple(data, band_pass, band_stop; fs=fs, mapslices2=mapslices2, dims=1, kwargs...)
#     clips_fixed = extract_clips(data_filt, window_extract .+ res.res_impulsetrain.pind_good, NaN; flag_matrix=false)
#     clips = extract_clips_single_channel(clips_fixed, ref_channel)
#     return clips, data_filt
# end

# rgbs, a = sig2rgb(clips[:,:]; fs=fs, 
#     rgb_bands=[[10_000, 60_000], [60_000, 110_000], [110_000, 160_000]])
# save("temp/test.png", a)

# Helper function for resizing
function nn_resize(S, F, T)
    fsrc, tsrc = size(S)
    out = Matrix{Float32}(undef, F, T)
    for j in 1:T
        tj = clamp(round(Int, (j-1)/(T-1) * (tsrc-1) + 1), 1, tsrc)
        for i in 1:F
            fi = clamp(round(Int, (i-1)/(F-1) * (fsrc-1) + 1), 1, fsrc)
            out[i,j] = S[fi, tj]
        end
    end
    out
end

function multispec_rgb_image(sig::AbstractVector, fs::Real;
    nffts = [128, 512, 2048],
    window=hann,
    norm=:per_channel,
    dynrange=80.0,
    gamma=1/2.2,
    savepath::Union{Nothing,String}=nothing
)
    # Compute spectrograms for each nfft
    S = Matrix{Float64}[]
    for nfft in nffts
        hop = nfft ÷ 4
        s = stft(sig, nfft, hop; window=window) .|> abs
        push!(S, s)
    end

    # Resize and normalize each spectrogram
    target_f = maximum(size(s,1) for s in S)
    target_t = maximum(size(s,2) for s in S)
    rgb = Array{Float32,3}(undef, target_f, target_t, 3)
    for c in 1:3
        # Resize
        rgb[:,:,c] = nn_resize(S[c], target_f, target_t)
        # Normalize
        mx = maximum(rgb[:,:,c])
        if mx > 0
            rgb[:,:,c] ./= mx
        end
    end

    # Convert to RGB image and flip vertically
    img = colorview(RGB, permutedims(rgb, (3,1,2))[:, end:-1:1, :])

    # Optionally save
    if !isnothing(savepath)
        save(savepath, img)
    end

    return img, rgb
end

# function autocorrelation_analysis(clips; threshold_autocor=4, threshold_n_autocor=3, nfunc=x->sqrt(sum(abs2.(x))), plot_dir::Union{Nothing,String}=nothing)
#     autocor_clips = Int[]
#     for i in eachindex(clips)
#         ac = xcorr(clips[:,i], clips[:,i]; padmode=:none)
#         ac = ac ./ maximum(abs.(ac))
#         ac = @view ac[length(clips[:,i]):end]
#         peaks, _ = findmaxima(ac)
#         peakproms!(peaks, ac; minprom=0.1)
#         n_peaks = length(peaks)
#         if n_peaks >= threshold_n_autocor && maximum(ac[peaks]) >= threshold_autocor/10
#             push!(autocor_clips, i)
#         end

#         if !isnothing(plot_dir)
#             plot(ac; title="Clip $i Autocorrelation (n_peaks=$n_peaks)", xlabel="Lag", ylabel="Normalized Amplitude")
#             scatter!(peaks, ac[peaks]; color=:red, label="Peaks")
#             savefig(joinpath(plot_dir, "clip_$(lpad(i,3,'0'))_autocorrelation.png"))
#             close("all")
#         end
#     end
#     return autocor_clips
# end
function autocor_analysis(clips, fs, ref_channel, train_start_ind)
    if size(clips[i], 2) > 1
        clips[i] = clips[i][:,ref_channel]
    end
    a = plot(signal(clips[i],fs); title=string(i))
    b = specgram(clips[i]; fs=fs, colorbar=nothing, nfft=128, crange=80)
    c = specgram(clips[i]; fs=fs, colorbar=nothing, nfft=round(Int,fs*.01)|>nextfastfft )
    selections = train_start_ind[i]:train_start_ind[i+1]-1
    snip = clips_fixed[:, selections]
    d = plot_time_fft(snip, fs; legend_position=:outerbottom, labels=reshape(string.(selections),1,length(selections)))

    correls = map(snip|>eachcol) do ref
        map(x-> mfilter(norm_max(ref; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy
    end

    e = Plots.bar(selections, map(x-> mfilter(norm_max(x; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy)
    Plots.hline!([threshold_autocor], label="Threshold", color=:red)

    for ftype = output_types
            savefig(joinpath(clips_plot_dir, "clip_$(i).$ftype"))
            if flag_extra_plot
            Plots.plot(d, e; layout=(2,1), legend_position=:outerbottom)#, size=(1200,800), title=string(i))
            savefig(joinpath(clips_plot_dir, "all_$(i).$ftype"))
            end
    end
    # savefig(joinpath(clips_plot_dir, "clip_$(i).html"))

    if count(map(x-> mfilter(norm_max(x; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy .> threshold_autocor) > threshold_n_autocor
        push!(autocor_clips, i)
    end
end


function analyze_clips(
    target,#::AbstractString,
    summary_fname::AbstractString,
    result_directory::AbstractString,
    impulsive_band_pass::AbstractVector,
    window_extract,
    threshold_autocor::Real=4,
    threshold_n_autocor::Int=3
    ;
    func = x->x,
    plot_dir_prefix::AbstractString="temp/clips_train_",
    output_types = [],#["html"],
    flag_extra_plot = false,
    fname2dt_func = DEFAULT_fname2timestamp_func,
    ref_channel=ref_channel
)
    if isfile(target)
        aufname = target
    else
        # Find closest row and get audio file
        result = find_closest_row(summary_fname, target)
        aufname = result.filepath
    end
    @info "Processing audio file: $aufname"

    if isfile(result_directory)
        respath = result_directory
    else
    # Find result path
        respath = readdir(result_directory; join=true) |>
            filter(isdir) .|> readdirjoin .|>
            filter(endswith(".jld2")) .|>
            filter(contains(splitext(basename(aufname))[1])) |>
            filter(!isempty) |> first
    end

    @info "loading result from: $respath"
    res = load(respath)
    res = dict2namedtuple(res)

    if res.res_impulsetrain.train_start |> isempty
        @warn "No impulsive train detected in the result."
        return Int[], "", threshold_autocor, threshold_n_autocor
    end

    # Read and filter audio
    data, fs, _, _, timestamp = readAudio(aufname; fname2timestamp_func=fname2dt_func)
    data_filt = filter_simple(data, impulsive_band_pass; fs=fs)

    # Extract clips
    clips_fixed = extract_clips(data_filt, window_extract .+ res.res_impulsetrain.pind_good, NaN; flag_matrix=true)
    clips = extract_clips(data_filt, [res.res_impulsetrain.train_start res.res_impulsetrain.train_end], NaN; flag_matrix=false)

    # Prepare output directory
    clips_plot_dir = plot_dir_prefix * Dates.format(timestamp, "yyyymmdd_HHMMSS")
    mkpath(clips_plot_dir)

    train_start_ind = [res.res_impulsetrain.train_start_ind... length(res.res_impulsetrain.pind_good)+1]
    nfunc(x) = sqrt(sum(abs2.(x)))
    autocor_clips = Int[]

    for i in eachindex(clips)
        out = func(clips[i])
        selections = train_start_ind[i]:train_start_ind[i+1]-1
        snip = clips_fixed[:, selections]

    end
    
    for ftype = output_types
        make_clip_index_html(clips_plot_dir; outname=basename(clips_plot_dir)*"_$(ftype)clips.html", output_type=ftype, prefix="clip_")
        if flag_extra_plot
            make_clip_index_html(clips_plot_dir; outname=basename(clips_plot_dir)*"_$(ftype)_all.html", output_type=ftype, prefix="all_")
        end
    end
    # make_clip_index_html(clips_plot_dir; outname=basename(clips_plot_dir)*".html", output_types)
    return out
end



function stretch(a, newmin::Real=0.0, newmax::Real=1.0; dim::Union{Nothing,Int}=nothing)
    """
    Linearly stretch/rescale array `a` to the interval [newmin, newmax].
    If `dim` is `nothing` (default) the whole array is rescaled.
    If `dim` is an integer, rescaling is performed independently for each slice along that dimension.

    Constant slices (min == max) become filled with the midpoint (newmin+newmax)/2.
    Returned array is Float64 for numeric inputs (to avoid integer truncation).
    """
    if newmax == newmin
        error("newmin and newmax must differ")
    end

    if dim === nothing
        af = Float64.(a)
        amin = minimum(af)
        amax = maximum(af)
        if amax == amin
            return fill((newmin + newmax)/2, size(a))
        end
        return (af .- amin) .* ((newmax - newmin) / (amax - amin)) .+ newmin
    else
        f = x -> begin
            xf = Float64.(x)
            amin = minimum(xf); amax = maximum(xf)
            if amax == amin
                fill((newmin + newmax)/2, size(xf))
            else
                (xf .- amin) .* ((newmax - newmin) / (amax - amin)) .+ newmin
            end
        end
        return mapslices(f, a; dims=dim)
    end
end



function rgb_from_value(val, vmin, vmax; flag_clamp::Bool=true)
    """
    Map a scalar `val` in [vmin, vmax] to an RGB triple (each in 0..1).
    The mapping uses an HSV hue sweep from red (0°) at vmin through green/yellow to blue (240°) at vmax.
    If `flag_clamp` is true (default) values outside [vmin, vmax] are clamped.
    Returns a 3-tuple (r,g,b) of Float64.
    """
    if vmax == vmin
        return (0.5, 0.5, 0.5)
    end

    t = (val - vmin) / (vmax - vmin)
    if flag_clamp
        t = clamp(t, 0.0, 1.0)
    end

    # Map t in [0,1] to hue in [0° -> 240°] (red -> blue)
    h = t * (240.0/360.0)  # normalized hue in [0,1]

    s = 1.0
    v = 1.0

    # HSV to RGB
    if s == 0.0
        return (v, v, v)
    end
    h6 = h * 6.0
    sector = floor(Int, h6) % 6
    f = h6 - floor(h6)
    p = v * (1 - s)
    q = v * (1 - s * f)
    tcol = v * (1 - s * (1 - f))

    if sector == 0
        return (v, tcol, p)
    elseif sector == 1
        return (q, v, p)
    elseif sector == 2
        return (p, v, tcol)
    elseif sector == 3
        return (p, q, v)
    elseif sector == 4
        return (tcol, p, v)
    else # sector == 5
        return (v, p, q)
    end
end