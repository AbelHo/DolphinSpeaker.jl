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
		@debug size(ref_signal)
		for ch in 1:size(data,2)
			@debug (i,ch)
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

function detection2angle(data, pind_good, rx_vect; fs=fs, window=default_window, ref_channel=ref_channel, channels_relevant=1:size(rx_vect,2),
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



default_getTDOA_func = get_tdoa_minmax #get_tdoa_raw_MaxEnergyRefChannel #get_tdoa_raw #get_tdoa_envelope#get_tdoa_raw_MaxEnergyRefChannel
include("config.jl")

@info "localization.jl loaded"

