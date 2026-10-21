using Optim
using SignalAnalysis.Units
include("dsp.jl")

default_tdoa2dir_solver = NelderMead
default_window = -5000:10000
# default_getTDOA_func = x->x; #get_tdoa_raw_MaxEnergyRefChannel
if (!@isdefined ref_channel) 
    ref_channel = 1
end
if !@isdefined default_c
    default_c =1540 #speed of sound
end

# include("config.jl")
get_tdoa_raw(data, peaks, ref_channel; window = default_window) = get_tdoa_raw(data, peaks; window = window , ref_channel=ref_channel)

function get_tdoa_raw(data, peaks; window = default_window , ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		for ch in 1:size(data,2)
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			@debug ([window.+peaks[i], ch ], [window.+peaks[i], ref_channel ])
			win = window.+peaks[i]
			if win[1]<1
				win = 1:win[end]
			elseif win[end]>size(data,1)
				win = win[1]:size(data,1)
			end
			tdoa[i,ch] = finddelay(data[win, ch ], data[win, ref_channel ])
		end
	end
	return tdoa
end

function get_tdoa_raw_MaxEnergyRefChannel(data, peaks; window = default_window , ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		win = window.+peaks[i]
		if win[1]<1
			win = 1:win[end]
		elseif win[end]>size(data,1)
			win = win[1]:size(data,1)
		end
		
		ref_channel = energy(data[win,:]) |> argmax
		for ch in 1:size(data,2)
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			tdoa[i,ch] = finddelay(data[win, ch ], data[win, ref_channel ])
		end
	end
	return tdoa
end

function get_tdoa_raw_MaxEnergyRefChannel_resample(data, peaks; window = default_window , ref_channel=ref_channel, resample_ratio=10)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		win = window.+peaks[i]
		if win[1]<1
			win = 1:win[end]
		elseif win[end]>size(data,1)
			win = win[1]:size(data,1)
		end
		
		ref_channel = energy(data[win,:]) |> argmax
		data_new = resample(data[win,:], resample_ratio; dims=1)
		for ch in 1:size(data,2)
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			tdoa[i,ch] = finddelay(data_new[:, ch ], data_new[:, ref_channel ])
		end
	end
	return tdoa ./ resample_ratio
end

function get_tdoa_raw_MaxPeakRefChannel(data, peaks; window = default_window , ref_channel=ref_channel, ref_window=-30:33)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	ref_signals = Array{eltype(data)}(undef,64, length(peaks))
	for i in eachindex(peaks)
		win = window.+peaks[i]
		win_ref = ref_window.+peaks[i]
		if win[1]<1
			win = 1:win[end]
			win_ref = 1:win_ref[end]
		elseif win[end]>size(data,1)
			win = win[1]:size(data,1)
			win_ref = win_ref[1]:size(data,1)
		end
		
		ref_channel = argmax(extrema(data[win,:], dims=1)[:] .|> x->maximum(abs.(x)))
		ref_signal = data[win_ref, ref_channel]
		ref_signals[:,i] = ref_signal
		# @debug size(ref_signal)
		for ch in 1:size(data,2)
			# @debug (i,ch)
			# if ch==ref_channel
			# 	tdoa[i,ch]=0
			# 	continue
			# end
			tdoa[i,ch] = finddelay(data[win, ch ], ref_signal)
		end
	end
	return tdoa#, ref_signals
end

function get_tdoa_raw_flexi(data, windows = default_window ; window = default_window, ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(windows),size(data,2))
	for i in eachindex(windows)
		window = windows[i]
		for ch in 1:size(data,2)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			tdoa[i,ch] = finddelay(data[window, ch ], data[window, ref_channel ])
		end
	end
	return tdoa
end

get_tdoa_envelope(data, peaks, ref_channel; window = default_window) =  get_tdoa_envelope(data, peaks; window = window , ref_channel=ref_channel)
function get_tdoa_envelope(data, peaks; window = default_window, ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		for ch in 1:size(data,2)
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			@debug window.+peaks[i]
			tdoa[i,ch] = finddelay(data[window.+peaks[i], ch ] |> hilbert .|> abs, data[window.+peaks[i], ref_channel ] |> hilbert .|> abs)
		end
	end
	tdoa
end

envelope_lpfiltered(data; lpf_maxFreq=10_000, fs=500_000) = filter_simple(abs.(hilbert(data)), [0 lpf_maxFreq]; fs=fs)

get_tdoa_envelope_filtered(data, peaks, ref_channel; window = default_window) =  get_tdoa_envelope(data, peaks; window = window , ref_channel=ref_channel)
function get_tdoa_envelope_filtered(data, peaks; window = default_window, ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		for ch in 1:size(data,2)
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			tdoa[i,ch] = finddelay(data[window.+peaks[i], ch ] |> envelope_lpfiltered, data[window.+peaks[i], ref_channel ] |> envelope_lpfiltered)
		end
	end
	tdoa
end

get_tdoa_max(data, peaks, ref_channel; window = default_window) = get_tdoa_max(data, peaks; window = window , ref_channel=ref_channel)
function get_tdoa_max(data, peaks; window = default_window , ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		maxi, max_ind = findmax( data[window.+peaks[i], : ], dims=1)
		# mini, min_ind = findmin( data[window.+peaks[i], : ], dims=1)

		tdoa[i,:] = getindex.(max_ind,1)[:] .- getindex(max_ind[ref_channel],1)
	end
	return tdoa
end

get_tdoa_min(data, peaks, ref_channel; window = default_window) = get_tdoa_min(data, peaks; window = window , ref_channel=ref_channel)
function get_tdoa_min(data, peaks; window = default_window , ref_channel=ref_channel, func=findmax)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		# maxi, max_ind = func( data[window.+peaks[i], : ], dims=1)
		maxi, max_ind = findmin( data[window.+peaks[i], : ], dims=1)

		tdoa[i,:] = getindex.(max_ind,1)[:] .- getindex(max_ind[ref_channel],1)
	end
	return tdoa
end


#~ tdoa to direction
# tdoa: Nxm, rx: complex
# angles: [azimuth elevation] in radians
# rx: vector 3xN
## N is number of ping, m is number of sensors
function cost_tdoa2ang(tdoa,ang,rx_vect,fs; c=default_c, ref_channel=ref_channel)
	# ref_channel = 1
	sv = [sin(ang[1])*cos(ang[2]) sin(ang[2]) cos(ang[1])*cos(ang[2])]
	dt = -sv * rx_vect ./c .* fs
	tdoa = tdoa .- tdoa[ref_channel]
	dt = dt' .- dt[ref_channel]# .+ tdoa[1]
	sum( (dt-tdoa).^2)
end

mod_sign(a,b) = sign(a)*mod(abs(a),b)

function tdoa2dir(tdoas,rx_vect,fs; 
	solver_func=default_tdoa2dir_solver, cost_tdoa2ang=cost_tdoa2ang, return_residual=false)
	
	initial_vals = zeros(size(tdoas,1),2) .+0.0
	for i in 1:size(tdoas,1)
		r2 = optimize(x -> cost_tdoa2ang(tdoas[i,:], x, rx_vect, fs), initial_vals[i,:], solver_func())#BFGS() #SimulatedAnnealing())#
		initial_vals[i,:] = Optim.minimizer(r2)
	end

	if return_residual
		return mod_sign.(initial_vals, 2*pi), cost_tdoa2ang.([tdoas[i,:] for i in axes(tdoas,1)], collect(eachrow(initial_vals)), Ref(rx_vect), Ref(fs))
	else
		return mod_sign.(initial_vals, 2*pi)
	end
end

function detection2angle(data, pind_good, rx_vect; fs=fs, window=default_window, ref_channel=ref_channel, channels_relevant=get_relevant_channels(rx_vect),
	getTDOA_func=default_getTDOA_func, solver_func=default_tdoa2dir_solver, cost_tdoa2ang=cost_tdoa2ang, return_residual=false)
    
	tdoa = getTDOA_func( @view(data[:,channels_relevant]), pind_good ; window=window, ref_channel=ref_channel)
    tdoa2dir(tdoa, rx_vect, fs; solver_func=solver_func, cost_tdoa2ang=cost_tdoa2ang, return_residual=return_residual), tdoa
	# nothing, tdoa
end

function detection2angle_envelope(data, pind_good, rx_vect; window=default_window , ref_channel=ref_channel, solver_func=default_tdoa2dir_solver)
    tdoa = get_tdoa_envelope(data, pind_good; window=window, ref_channel=ref_channel)
    tdoa2dir(tdoa, rx_vect), tdoa
end


function angle2px(ang, fov_angle=fov_angle, imsize=imsize)
	@debug angle2pxd( rad2deg.(ang), fov_angle, imsize)
	return angle2pxd( rad2deg.(ang), fov_angle, imsize)
end

function angle2pxd(angd, fov_angled=fov_angle, imsize=imsize)
	return angd ./ [fov_angled[1]/2 -fov_angled[2]/2] .* [imsize[2]/2 imsize[1]/2] .+ [imsize[2]/2 imsize[1]/2]
end



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

function detect_angle(aufname_data_fs::Tuple, rx_vect; 
	res_dir=nothing, ref_channel=1, dist=dist_impulsive, window=window_impulsive, threshold=threshold_impulsive,
	tdoa2dir=tdoa2dir, cost_tdoa2ang=cost_tdoa2ang, solver_func=default_tdoa2dir_solver, getTDOA_func=default_getTDOA_func)

    aufname, data, fs = aufname_data_fs

    detections = detect_impulse(aufname_data_fs, res_dir;
        ref_channel=ref_channel, dist=dist, window=window, threshold=threshold)
    tdoa = getTDOA_func(data, detections.pind_good ; window=window, ref_channel=ref_channel)
    ang = tdoa2dir(tdoa, rx_vect; solver_func=solver_func, cost_tdoa2ang=cost_tdoa2ang)
end

#~ V2
#~ tdoa to direction
# tdoa: Nxm, rx: complex
# angles: [azimuth elevation] in radians
# rx: vector 3xN
## N is number of ping, m is number of sensors
function cost_tdoa2ang2(tdoa,ang,rx_vect,fs; c=default_c, ref_channel=ref_channel)
	# ref_channel = 1
	ang = ang[1]
	sv = [sin(ang) cos(ang) 0]
	dt = sv * rx_vect ./c .* fs
	tdoa = tdoa .- tdoa[ref_channel]
	dt = dt' .- dt[ref_channel]# .+ tdoa[1]
	sum( (dt-tdoa).^2)
end

function tdoa2dir2(tdoas,rx_vect,fs; solver_func=default_tdoa2dir_solver)
	initial_vals = zeros(size(tdoas,1),1) .+0.0
	for i in 1:size(tdoas,1)
		r2 = optimize(x -> cost_tdoa2ang2(tdoas[i,:], x, rx_vect, fs), [initial_vals[i]], solver_func())#BFGS() #SimulatedAnnealing())#
		initial_vals[i] = Optim.minimizer(r2)[1]
	end
	mod_sign.(initial_vals, 2*pi)
end

function detection2angle2(data, pind_good, rx_vect; fs=fs, window=default_window, ref_channel=ref_channel, getTDOA_func=default_getTDOA_func, solver_func=default_tdoa2dir_solver)
    tdoa = getTDOA_func(data, pind_good ; window=window, ref_channel=ref_channel)
    tdoa2dir2(tdoa, rx_vect, fs), tdoa
	# nothing, tdoa
end

function cost_tdoa2ang3(tdoa,ang,rx_vect,fs; c=default_c, ref_channel=ref_channel)
	# ref_channel = 1
	# ang = ang[1]
	sv = [sin(ang[1])*cos(ang[2]) cos(ang[1])*cos(ang[2]) sin(ang[2])]
	dt = sv * rx_vect ./c .* fs
	tdoa = tdoa .- tdoa[ref_channel]
	dt = dt' .- dt[ref_channel]# .+ tdoa[1]
	sum( (dt-tdoa).^2)
end

get_tdoa_minmax(data, peaks, ref_channel=nothing; window = default_window) = get_tdoa_max(data, peaks; window = window , ref_channel=ref_channel)
function get_tdoa_minmax(data, peaks; window = default_window , ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		ex = extrema_and_indices(data[window.+peaks[i], : ])
		max_or_min = ex .|> x -> x[2] > x[1]
		@debug max_or_min
		if sum(max_or_min) > length(ex)/2
			tdoa[i,:] =  ex .|> x -> x[4]
		else
			tdoa[i,:] =  ex .|> x -> x[3]
		end
	end
	return tdoa
end

function get_tdoa_findTrigger(data, peaks; window = default_window , ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		win = window.+peaks[i]
		if win[1]<1
			win = 1:win[end]
		elseif win[end]>size(data,1)
			win = win[1]:size(data,1)
		end
		
		ref_channel = energy(data[win,:]) |> argmax
		for ch in 1:size(data,2)
			@debug (i,ch)
			# if ch==ref_channel
			# 	tdoa[i,ch]=0
			# 	continue
			# end
			tdoa[i,ch] = findTrigger(data[win, ch ], 1);
		end
	end
	return tdoa
end

function get_tdoa_findTrigger_waveletdenoise(data, peaks; window = default_window , ref_channel=ref_channel)
	tdoa = Array{Int}(undef,length(peaks),size(data,2))
	for i in eachindex(peaks)
		win = window.+peaks[i]
		if win[1]<1
			win = 1:win[end]
		elseif win[end]>size(data,1)
			win = win[1]:size(data,1)
		end
		
		ref_channel = energy(data[win,:]) |> argmax
		Threads.@threads for ch in 1:size(data,2)
			@debug (i,ch)
			@debug win
			# if ch==ref_channel
			# 	tdoa[i,ch]=0
			# 	continue
			# end
			tdoa[i,ch] = findTrigger(denoise(data[win, ch ], TI=true), 1);
		end
	end
	return tdoa
end

#~ 
ang2vect(ang) = [sin(ang[1])*cos(ang[2]) sin(ang[2]) cos(ang[1])*cos(ang[2])]

function compensate_imu(angs, datetimes, df)
	indx = searchsortedfirst.(Ref(ahrs.time), ZonedDateTime(timestamp) .+ map( x -> Nanosecond(round(x * 1e9)), res.res_impulsetrain.pind_good_inS))
	ang_compensated = map( x -> x = x > (-2*pi +pi/4) ? x : x+2pi, ang .- deg2rad.(ahrs[indx,:yaw]));
	
end

function steering2(rxpos::AbstractMatrix, c, θ::AbstractMatrix)
	nsensors = size(rxpos, 2)
	ndir = size(θ, 1)
	pos0 = sum(rxpos; dims=2) / nsensors
	rxpos = rxpos .- pos0
	size(rxpos, 1) == 1 && (rxpos = vcat(rxpos, zeros(2, nsensors)))
	size(rxpos, 1) == 2 && (rxpos = vcat(rxpos, zeros(1, nsensors)))
	γ = θ[:,1]
	ϕ = θ[:,2]
	sd = zeros(nsensors, ndir)
	for k ∈ 1:ndir
	  for j ∈ 1:nsensors
		sd[j,k] = -rxpos[:,j]' * [cos(ϕ[k])*sin(γ[k]), sin(ϕ[k]), -cos(ϕ[k])*cos(γ[k])] / c
	  end
	end
	sd
  end
  
  steering2(rxpos::AbstractVector, c, θ) = steering(collect(rxpos'), c, θ)
  

function beamform_output(s, sd; fs=framerate(s), output_func=p2p_db)
	ndir = size(sd, 2)
	bfo = zeros(eltype(s), ndir)
	Threads.@threads for k ∈ 1:ndir
	  temp_bfo = zeros(eltype(s), nframes(s))
	  for j ∈ 1:nchannels(s)
		@views temp_bfo .+= padded(s[:,j], 0; delay=round(Int, -sd[j,k]*fs))
	  end
	  bfo[k] = output_func(temp_bfo)
	end
	bfo
	# signal(bfo, fs)
  end



#~ sub-sample TDOA: GCC / SRP-PHAT and a closed-form direction solver
# plan: reviews/plans/trace-the-function-run-contiguous-folder-eager-puppy.md
# TDOA convention as finddelay(data[:,ch], data[:,ref]): positive = channel arrives later than ref
using FFTW, LinearAlgebra, Statistics

_clip_window(win, n) = max(first(win), 1):min(last(win), n)

# largest physically possible lag (samples) for each sensor pair, plus a margin
tau_max_samples(rx_vect, fs; c=default_c, margin=2) = distance_pairwise(rx_vect) ./ c .* fs .+ margin

# geometry of the channels in `data` (rx_vect may still hold the NaN columns of unused channels)
_rx_for(data, rx_vect) = size(rx_vect, 2) == size(data, 2) ? rx_vect : rx_vect[:, get_relevant_channels(rx_vect)]

# band for the GCC mask, top capped below Nyquist
default_gcc_band(band, fs) = (b = something(band, impulsive_band_pass); [b[1], min(b[2], 0.98 * fs / 2)])

# parabolic interpolation around index k; returns (fractional index, interpolated peak value)
function subsample_peak(s, k)
	(k <= firstindex(s) || k >= lastindex(s)) && return float(k), s[k]
	y0, y1, y2 = s[k-1], s[k], s[k+1]
	den = y0 - 2y1 + y2
	iszero(den) && return float(k), y1
	δ = clamp(0.5 * (y0 - y2) / den, -0.5, 0.5)
	return k + δ, y1 - 0.25 * (y0 - y2) * δ
end

# mean |rfft|^2 of consecutive length-L segments of n, zero-padded to N (same scaling as the click snippet)
function noise_psd(n, L, N; taper=nothing)
	nseg = length(n) ÷ L
	P = zeros(N ÷ 2 + 1)
	for s in 1:nseg
		seg = Float64.(n[(s-1)*L+1:s*L])
		isnothing(taper) || (seg .*= taper)
		P .+= abs2.(rfft([seg; zeros(N - L)]))
	end
	return P ./ max(nseg, 1)
end

# frequency weighting of the cross-spectrum R = X conj(Y)
# :scot uses auto-spectra smoothed over `smooth_bins` bins (unsmoothed, it equals :phat for a single click)
# :ml weights each bin by its coherence from a noise spectrum (falls back to :phat without noise)
function gcc_weight(weighting, X, Y, R; L, N, noise_x=nothing, noise_y=nothing, taper=nothing, smooth_bins=5)
	ε = 1e-12 * maximum(abs, R) + floatmin(Float64)
	weighting == :none && return 1.0
	weighting == :phat && return 1 ./ (abs.(R) .+ ε)
	if weighting == :scot
		smooth(P) = (h = smooth_bins ÷ 2; [mean(P[max(1, i - h):min(end, i + h)]) for i in eachindex(P)])
		return 1 ./ (sqrt.(smooth(abs2.(X)) .* smooth(abs2.(Y))) .+ ε)
	end
	if weighting == :ml
		(isnothing(noise_x) || isnothing(noise_y)) && return 1 ./ (abs.(R) .+ ε)
		Nx = noise_psd(noise_x, L, N; taper=taper)
		Ny = noise_psd(noise_y, L, N; taper=taper)
		Sx = max.(abs2.(X) .- Nx, 0.0)
		Sy = max.(abs2.(Y) .- Ny, 0.0)
		γ2 = min.(Sx .* Sy ./ ((Sx .+ Nx) .* (Sy .+ Ny) .+ ε^2), 0.99)
		return γ2 ./ ((1 .- γ2) .* (abs.(R) .+ ε))
	end
	error("unknown weighting: $weighting")
end

# weighted cross-correlation of x against y over lags |τ| <= maxlag (samples), upsampled by `interp`
# returns (lags, r); r is normalized so that 1 = all weighted bins in phase; peak at τ > 0: x later than y
function gcc_curve(x, y; fs, maxlag, weighting=:phat, band=nothing, interp=8, noise_x=nothing, noise_y=nothing, taper=nothing)
	L = length(x)
	N = nextfastfft(2L)
	X = rfft([Float64.(x); zeros(N - L)])
	Y = rfft([Float64.(y); zeros(N - L)])
	R = X .* conj.(Y)
	R = R .* gcc_weight(weighting, X, Y, R; L=L, N=N, noise_x=noise_x, noise_y=noise_y, taper=taper)
	if !isnothing(band)
		f = rfftfreq(N, fs)
		R[(f .< band[1]) .| (f .> band[2])] .= 0
	end
	Nup = N * interp
	Rup = zeros(ComplexF64, Nup ÷ 2 + 1)
	Rup[1:length(R)] .= R
	r = irfft(Rup, Nup) ./ max(2 * sum(abs, R) / Nup, floatmin(Float64))
	M = min(ceil(Int, maxlag * interp), Nup ÷ 2 - 1)
	return (-M:M) ./ interp, r[mod.(-M:M, Nup) .+ 1]
end

# linear interpolation of a gcc_curve result at lag τ (clamped to its lag range)
function curve_at((lags, r), τ)
	x = clamp((τ - first(lags)) / step(lags) + 1, 1, length(r))
	k = min(floor(Int, x), length(r) - 1)
	return r[k] + (x - k) * (r[k+1] - r[k])
end

# delay (samples) of x relative to y and its normalized peak height; with `guide` (samples), take the
# local maximum nearest to it among those reaching guide_ratio of the highest one
function gcc_delay(x, y; guide=nothing, guide_ratio=0.7, kwargs...)
	lags, r = gcc_curve(x, y; kwargs...)
	k = argmax(r)
	if !isnothing(guide)
		cand = [i for i in 2:length(r)-1 if r[i] >= r[i-1] && r[i] >= r[i+1] && r[i] >= guide_ratio * r[k]]
		isempty(cand) || (k = cand[argmin(abs.(lags[cand] .- guide))])
	end
	kf, peak = subsample_peak(r, k)
	return first(lags) + (kf - 1) * step(lags), peak
end

# per-channel delays (sum zero) from pairwise delays τ[p] = t[j] - t[k] of pairs[p] = (j, k):
# least squares when all pairs are given, else reference pairs (j, ref) with t[ref] = 0
# also returns the rms pair misfit (for 3 sensors the closure |τ12+τ23+τ31| is 3x this)
function pairs_to_delays(τ, pairs, m)
	t = zeros(m)
	if length(pairs) == m * (m - 1) ÷ 2
		for (p, (j, k)) in enumerate(pairs)
			t[j] += τ[p]
			t[k] -= τ[p]
		end
		t ./= m
	else
		for (p, (j, _)) in enumerate(pairs)
			t[j] = τ[p]
		end
	end
	return t, sqrt(sum(p -> abs2(t[pairs[p][1]] - t[pairs[p][2]] - τ[p]), eachindex(pairs)) / length(pairs))
end

# GCC TDOA: frequency weighting (:phat, :scot, :ml or :none), band mask, physical lag limit and sub-sample peak
# all_pairs: estimate every pair and reconcile them by least squares; guide=:envelope: on near-equal peaks,
# take the one nearest the envelope-correlation delay; noise for :ml is taken noise_len windows before the click
function get_tdoa_gcc(data, peaks; window=default_window, ref_channel=ref_channel, fs, rx_vect=rx_vect, c=default_c,
	weighting=:phat, band=nothing, interp=8, all_pairs=true, guide=nothing, taper=0.0, noise_len=4, return_quality=false)

	m = size(data, 2)
	τmax = tau_max_samples(_rx_for(data, rx_vect), fs; c=c)
	band = default_gcc_band(band, fs)
	pairs = all_pairs ? [(j, k) for j in 1:m for k in j+1:m] : [(j, ref_channel) for j in 1:m if j != ref_channel]
	tdoa = zeros(length(peaks), m)
	peak = zeros(length(peaks), length(pairs))
	pair_rms = zeros(length(peaks))
	Threads.@threads for i in eachindex(peaks)
		win = _clip_window(window .+ peaks[i], size(data, 1))
		L = length(win)
		w = taper > 0 ? tukey(L, taper) : nothing
		snip = Float64.(data[win, :])
		isnothing(w) || (snip .*= w)
		noise = nothing
		if weighting == :ml
			nwin = _clip_window((first(win) - noise_len * L):(first(win) - 1), size(data, 1))
			length(nwin) >= L && (noise = data[nwin, :])
		end
		env = guide == :envelope ? reduce(hcat, [abs.(hilbert(snip[:, j])) for j in 1:m]) : nothing
		τ = zeros(length(pairs))
		for (p, (j, k)) in enumerate(pairs)
			g = isnothing(env) ? nothing :
				gcc_delay(env[:, j] .- mean(env[:, j]), env[:, k] .- mean(env[:, k]); fs=fs, maxlag=τmax[j, k], weighting=:none, interp=interp)[1]
			τ[p], peak[i, p] = gcc_delay(snip[:, j], snip[:, k]; fs=fs, maxlag=τmax[j, k], weighting=weighting, band=band, interp=interp,
				guide=g, taper=w, noise_x=isnothing(noise) ? nothing : noise[:, j], noise_y=isnothing(noise) ? nothing : noise[:, k])
		end
		t, pair_rms[i] = pairs_to_delays(τ, pairs, m)
		tdoa[i, :] = t .- t[ref_channel]
	end
	return return_quality ? (tdoa, (; peak, pair_rms, pairs)) : tdoa
end

# SRP-PHAT: the direction (u, v) = (sin(az)cos(el), sin(el)) maximizing the summed pairwise GCC curves,
# grid_deg coarse grid then a local grid `refine` times finer; returns the TDOAs predicted at that direction
function get_tdoa_srp(data, peaks; window=default_window, ref_channel=ref_channel, fs, rx_vect=rx_vect, c=default_c,
	weighting=:phat, band=nothing, interp=8, grid_deg=1.0, refine=10, taper=0.0, return_quality=false)

	m = size(data, 2)
	rx = _rx_for(data, rx_vect)
	τmax = tau_max_samples(rx, fs; c=c)
	band = default_gcc_band(band, fs)
	pairs = [(j, k) for j in 1:m for k in j+1:m]
	arrival(u, v) = -(rx' * [u, v, sqrt(max(1 - u^2 - v^2, 0.0))]) ./ c .* fs   # as in cost_tdoa2ang
	pairdelays(u, v) = (t = arrival(u, v); [t[j] - t[k] for (j, k) in pairs])
	h = sind(grid_deg)
	grid = [(u, v) for u in -1:h:1 for v in -1:h:1 if u^2 + v^2 <= 1]
	D = reduce(hcat, [pairdelays(u, v) for (u, v) in grid])
	fine = range(-h, h; length=2refine + 1)

	tdoa = zeros(length(peaks), m)
	srp = zeros(length(peaks))
	uv = zeros(length(peaks), 2)
	Threads.@threads for i in eachindex(peaks)
		win = _clip_window(window .+ peaks[i], size(data, 1))
		snip = Float64.(data[win, :])
		taper > 0 && (snip .*= tukey(length(win), taper))
		curves = [gcc_curve(snip[:, j], snip[:, k]; fs=fs, maxlag=τmax[j, k], weighting=weighting, band=band, interp=interp) for (j, k) in pairs]
		score(d) = sum(p -> curve_at(curves[p], d[p]), eachindex(pairs))
		g = argmax(gi -> score(view(D, :, gi)), eachindex(grid))
		best, (bu, bv) = score(view(D, :, g)), grid[g]
		u0, v0 = bu, bv
		for du in fine, dv in fine
			u, v = u0 + du, v0 + dv
			u^2 + v^2 <= 1 || continue
			s = score(pairdelays(u, v))
			s > best && ((best, bu, bv) = (s, u, v))
		end
		t = arrival(bu, bv)
		tdoa[i, :] = t .- t[ref_channel]
		srp[i] = best / length(pairs)
		uv[i, :] = [bu, bv]
	end
	return return_quality ? (tdoa, (; srp, uv)) : tdoa
end

# unclamped direction cosines (u, v) = (sin(az)cos(el), sin(el)) solving the TDOAs of a planar array
# (TDOAs are linear in u, v, see cost_tdoa2ang); exact for 3 sensors, least squares for more
function tdoa2uv(tdoas, rx_vect, fs; c=default_c, ref_channel=ref_channel)
	others = setdiff(axes(rx_vect, 2), ref_channel)
	A = -(rx_vect[1:2, others] .- rx_vect[1:2, ref_channel])' ./ c .* fs
	return reduce(vcat, [(A \ (tdoas[i, others] .- tdoas[i, ref_channel]))' for i in axes(tdoas, 1)])
end

# closed-form tdoa2dir for planar arrays (same call and output); infeasible TDOAs (u²+v² > 1) are
# clamped to the unit circle, which shows as a non-zero residual; non-planar arrays fall back to tdoa2dir
function tdoa2dir_ls(tdoas, rx_vect, fs; c=default_c, ref_channel=ref_channel, return_residual=false, cost_tdoa2ang=cost_tdoa2ang, kwargs...)
	z = rx_vect[3, :]
	all(abs.(z .- z[1]) .< 1e-9) ||
		return tdoa2dir(tdoas, rx_vect, fs; cost_tdoa2ang=cost_tdoa2ang, return_residual=return_residual, kwargs...)
	uv = tdoa2uv(tdoas, rx_vect, fs; c=c, ref_channel=ref_channel)
	angs = zeros(size(tdoas, 1), 2)
	for i in axes(uv, 1)
		u, v = uv[i, 1], uv[i, 2]
		r = hypot(u, v)
		r > 1 && (u /= r; v /= r)
		el = asin(clamp(v, -1, 1))
		angs[i, :] = [cos(el) > 1e-12 ? asin(clamp(u / cos(el), -1, 1)) : 0.0, el]
	end
	angs = mod_sign.(angs, 2 * pi)
	return_residual || return angs
	return angs, cost_tdoa2ang.([tdoas[i, :] for i in axes(tdoas, 1)], collect(eachrow(angs)), Ref(rx_vect), Ref(fs))
end

#~ fixed copies of existing TDOA functions (the originals are left untouched)
# get_tdoa_minmax with the intended polarity vote, |max| vs |min| (get_tdoa_minmax compares max with min,
# which is always true, so it is plain argmax), and window clipping at the file edges
function get_tdoa_minmax2(data, peaks; window=default_window, ref_channel=ref_channel)
	tdoa = Array{Int}(undef, length(peaks), size(data, 2))
	for i in eachindex(peaks)
		ex = extrema_and_indices(data[_clip_window(window .+ peaks[i], size(data, 1)), :])
		use_max = count(x -> x[2] > -x[1], ex) > length(ex) / 2
		tdoa[i, :] = ex .|> x -> use_max ? x[4] : x[3]
	end
	return tdoa
end

# get_tdoa_envelope_filtered with the real fs passed to envelope_lpfiltered (whose default assumes 500 kHz)
# and window clipping at the file edges
function get_tdoa_envelope_filtered2(data, peaks; window=default_window, ref_channel=ref_channel, fs, lpf_maxFreq=10_000)
	tdoa = Array{Int}(undef, length(peaks), size(data, 2))
	for i in eachindex(peaks)
		win = _clip_window(window .+ peaks[i], size(data, 1))
		env_ref = envelope_lpfiltered(data[win, ref_channel]; lpf_maxFreq=lpf_maxFreq, fs=fs)
		for ch in 1:size(data, 2)
			tdoa[i, ch] = ch == ref_channel ? 0 : finddelay(envelope_lpfiltered(data[win, ch]; lpf_maxFreq=lpf_maxFreq, fs=fs), env_ref)
		end
	end
	return tdoa
end


default_getTDOA_func = get_tdoa_minmax #get_tdoa_raw_MaxEnergyRefChannel #get_tdoa_raw #get_tdoa_envelope#get_tdoa_raw_MaxEnergyRefChannel
include("config.jl")

@info "localization.jl loaded"

