using Pkg
Pkg.activate(".")
using FileIO, FLAC, DSP, SignalAnalysis, SignalAnalysis.Units
using Optim
using VideoIO, FFMPEG
using Plots
# plotlyjs()
include("pinger_calibration.jl")
include("../../../research/analysis_interface/Scripts/readImages.jl")

# aufname = "/Users/abel/Documents/data/aspod/field/bahamas_2022/Aud_2022-06-25_10.33.28.wav"
# vidfname = "/Users/abel/Documents/data/aspod/field/bahamas_2022/Vid_2022-06-25_10.33.28.mkv"
# #"/Volumes/PortableSSD/aspod/aspod4_2/2022-09-11_BoatTest_wav/Aud_B001C0044_20220912040444_0001_20131219_021518.wav"
# #"/Users/abel/Downloads/temp/data/aspod/maui_2022/2022.03.04/Aud_A003C0154_20220304130251_0001_20131219_100827.wav" #"/Users/abel/Downloads/temp/data/aspod/test/a4_2/2022-09-11_BoatTest_flac/Aud_B001C0044_20220912040444_0001_20131219_021518.flac"
# rx = 0.14722/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) )
# rx_vect = [real(rx); imag(rx); zeros(1,3)]
# #=
#   0.07361    5.2046e-18  -0.07361
#  -0.0424988  0.0849975   -0.0424988
#   0.0        0.0          0.0
# =#

# c = 1540 # m/s speed of sound
# dist = 15000#pinger
# rx_vect = [real(rx); imag(rx); zeros(1,3)]
# Plots.scatter(rx_vect[1,:],rx_vect[2,:],rx_vect[3,:])
# xlabel!("X");ylabel!("Y")

# vid = VideoIO.openvideo(vidfname)
# data, fs = load(aufname)
# println("Duration: " * string(size(data,1)/fs) *"seconds")
# pind, ppeak = findPings(data|>hilbert.|>abs; ref_channel=3, dist=dist)
# # p = findPings(data|>hilbert.|>abs; ref_channel=3, dist=15000)

# thresh = 0.1 # 0.1#pinger/clickler #0.2
# threshold_indices = findall(>(thresh), ppeak)
# pind_good = pind[threshold_indices]
# pind_good_inS = (pind_good.-1) ./fs
# # overthresh=filter(x -> x>thresh, ppeak)

# res_dir = "/Users/abel/Documents/data_res/aspod/real"
# isdir(res_dir) || mkpath(res_dir)
# audacity_label(pind_good./fs, joinpath(res_dir, splitext(aufname)[1]*".txt" |> basename))

# tdoa = get_tdoa_envelope(data, pind_good )
# tdoa_raw = get_tdoa_raw(data, pind_good )
# ang = tdoa2dir(tdoa, rx_vect)
# # plot( pind_good/fs, (mod.(directions,2*pi) .|> rad2deg) .* sign.(directions) )

##=#

# if isempty(eventsFilePath) 
# 	peaks, eventTimings, eventSignals = findPings(hydroAudioFile; save_fn=nothing,
#    dist=Int(400000*0.9), channels=5)
# else
# 	eventsFile = load(eventsFilePath)
# 	peaks = eventsFile["peaks"]
# 	eventTimings = eventsFile["peak_times"]
# 	eventSignals = eventsFile["peak_vals"]
# end;

function get_tdoa_raw(data, peaks; window = -5000:10000, ref_channel=3)
	tdoa = Array{Int}(undef,length(peaks),size(data,2)-1)
	for i in eachindex(peaks)
		for ch in 1:size(data,2)-1
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			tdoa[i,ch] = finddelay(data[window.+peaks[i], ch ], data[window.+peaks[i], ref_channel ])
		end
	end
	tdoa
end

function get_tdoa_envelope(data, peaks; window = -5000:10000, ref_channel=3)
	tdoa = Array{Int}(undef,length(peaks),size(data,2)-1)
	for i in eachindex(peaks)
		for ch in 1:size(data,2)-1
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			tdoa[i,ch] = finddelay(data[window.+peaks[i], ch ] |> hilbert .|> abs, data[window.+peaks[i], ref_channel ] |> hilbert .|> abs)
		end
	end
	tdoa
end


function get_tdoa_raw(data, peaks; window = -5000:10000, ref_channel=3)
	tdoa = Array{Int}(undef,length(peaks),size(data,2)-1)
	for i in eachindex(peaks)
		for ch in 1:size(data,2)-1
			@debug (i,ch)
			if ch==ref_channel
				tdoa[i,ch]=0
				continue
			end
			tdoa[i,ch] = finddelay(data[window.+peaks[i], ch ], data[window.+peaks[i], ref_channel ])
		end
	end
	return tdoa
end


function ding()
	#find Envelope TDoA
	ping_thresh = 400
	
	# data=signal(loadDataBin(hydroAudioFile), fs);
	# data2 = loadDataBin(hydroAudioFile);
	ref_channel = 1 
	channels = 5
	
	tdoa = Array{Int}(undef,length(peaks[1]),channels-1)
	for i in 1:length(peaks[1])
		for ch in 1:channels-1
			# if ch==ref_channel
			# 	tdoa[i,ch]=0
			# 	continue
			# end
			# tdoa[i,ch] = finddelay(data[window[1]+eventTimings[i]s : window[2]+eventTimings[i]s, ch ] |> hilbert .|> abs, data[window[1]+eventTimings[i]s : window[2]+eventTimings[i]s, ref_channel ] |> hilbert .|> abs)
			# tdoa[i,ch] = finddelay(data2[window2.+peaks[1][i], ch ] |> hilbert .|> abs, data2[window2.+peaks[1][i], ref_channel ] |> hilbert .|> abs)
			
			res = findfirst(  (data[window[1]+eventTimings[i]s : window[2]+eventTimings[i]s, ch ] |> hilbert .|> abs ) .> ping_thresh)
			if isnothing(res)
				res = findfirst(  (data[window[1]+eventTimings[i]s : window[2]+eventTimings[i]s, ch ] |> hilbert .|> abs ) .> ping_thresh/2 )
				isnothing(res) ? tdoa[i,ch] = -99999999 : tdoa[i,ch]=res
				
			else
				tdoa[i,ch] = res
			end
		end
	end
end



#~ tdoa to direction
# tdoa: Nxm, rx: complex
# angles: [azimuth elevation] in radians
# rx: vector 3xN
## N is number of ping, m is number of sensors
function cost_ray(angles,rx, ref_chan,fs,tdoa)
	cost = 0
	c= 1500
#     i=1
#     angles=Array{Float64}(undef,size(tdoa,1),2)
	for i in 1:size(tdoa,1)
		ang_vect = [cos(angles[2,i])*sin(angles[1,i])  sin(angles[2,i])  cos(angles[2,i])*cos(angles[1,i])]
		delays = ang_vect * rx /c*fs 

		delays = delays .- delays[ref_chan]
		cost += (tdoa[i,:].-tdoa[i,ref_chan] -delays).^2 |> sum
	end

	cost
end

function cost_ray_1var(ang_rx, ref_chan,fs,tdoa)
	#     ref_chan=1
	#     fs=400000
		
	cost_ray(ang_rx[1:2,1:size(tdoa,1)],ang_rx[:,size(tdoa,1)+1:end], ref_chan,fs,tdoa)
end

function cost_tdoa2ang(tdoa,ang,rx_vect,fs; c=1540, ref_channel=1)
	# ref_channel = 1
	sv = [sin(ang[1])*cos(ang[2]) sin(ang[2]) cos(ang[1])*cos(ang[2])]
	dt = -sv * rx_vect ./c .* fs
	tdoa = tdoa .- tdoa[ref_channel]
	dt = dt' .- dt[ref_channel]# .+ tdoa[1]
	sum( (dt-tdoa).^2)
end

mod_sign(a,b) = sign(a)*mod(abs(a),b)
function tdoa2dir(tdoas,rx_vect,fs; solver_func=NelderMead)
	initial_vals = zeros(size(tdoas,1),2) .+0.0
	for i in 1:size(tdoas,1)
		r2 = optimize(x -> cost_tdoa2ang(tdoas[i,1:3], x, rx_vect, fs), initial_vals[i,:], solver_func())#BFGS() #SimulatedAnnealing())#
		initial_vals[i,:] = Optim.minimizer(r2)
	end
	mod_sign.(initial_vals, 2*pi)
end

function vline2(xvals, previus_plot, max_y=3000)
	previus_plot
	for ind=1:length(xvals)
		Plots.plot!([xvals[ind]; xvals[ind]], [0; max_y*(1-0.05*ind)], color=palette(:default)[ind] , width=3)
	end
	
	# for x in xvals
	# 	a=Plots.plot!([x; x], [0; max_y])
	# end
	return previus_plot
end

function plotTDOA(data, i, eventTimings, tdoa; window = [-5000, 10000], func=plotlyjs, func_on_data=x->x)
	# i=selected_events[event]
	# using Plots
	func
	# specgram(data[ 0.0s+eventTimings[i]s : 2.0s+eventTimings[i]s ,1])
	# print(tdoa[i,:])
	# d = data[ window[1]+eventTimings[i]s : window[2]+eventTimings[i]s,1:4] |> hilbert .|> abs
	d = data[ window[1]+eventTimings[i] : window[end]+eventTimings[i],1:size(tdoa,2)] |> hilbert .|> abs
	a = Plots.plot(d)
	vline2( tdoa[i,:].-window[1] , a, maximum(d))
	# xlims!( window[1], window[end])
	# Plots.xlims!( minimum(tdoa[i,:])*1000 -0.5, maximum(tdoa[i,:])*1000+0.5)
	# Plots.xlims!( minimum(tdoa[i,:])/fs*1000 -0.5, maximum(tdoa[i,:])/fs*1000+0.5)
	# Plots.plot!([20; 20],[0; 3000])
	title!(string(i)*" tdoa:"*string(tdoa[i,:]))
end

function plotTDOA_raw(data, i, eventTimings, tdoa; window = [-5000, 10000], func=plotlyjs)
	# i=selected_events[event]
	# using Plots
	# plotlyjs()
	func()
	# specgram(data[ 0.0s+eventTimings[i]s : 2.0s+eventTimings[i]s ,1])
	# print(tdoa[i,:])"
	# d = data[ window[1]+eventTimings[i]s : window[2]+eventTimings[i]s,1:4] |> hilbert .|> abs
	d = data[ window[1]+eventTimings[i] : window[end]+eventTimings[i],1:size(tdoa,2)]
	a = Plots.plot(d)
	@debug tdoa[i,:]
	@debug window[1]
	vline2( tdoa[i,:].-window[1] , a, maximum(d))
	# xlims!( window[1], window[end])

	# xlims!( -window[1]-50, -window[1]+100)
	# Plots.xlims!( minimum(tdoa[i,:])*1000 -0.5, maximum(tdoa[i,:])*1000+0.5)
	# Plots.xlims!( minimum(tdoa[i,:])/fs*1000 -0.5, maximum(tdoa[i,:])/fs*1000+0.5)
	# Plots.plot!([20; 20],[0; 3000])
	title!(string(i)*" tdoa:"*string(tdoa[i,:]))
end


function plotSpecgram(data, i, eventTimings, tdoa; window = [-5000, 10000], func=plotlyjs)

	func()
	a = specgram(data[ window[1]+eventTimings[i] : window[end]+eventTimings[i], 1]; colorbar=nothing, nfft=512, fs=400000)#Plots.plot(d)
	maxi = 30
	ylims!(0,maxi)
	# vline2( (tdoa[i,:].-minimum(tdoa[i,:])) .* maxi , a, maxi)

	title!(string(i)*" tdoa:"*string(tdoa[i,:]))
end



function angle2px(ang, fov_angle=[62.61721188568244,35.793211268714096,71.6855447884958], imsize=(2160, 3840))
	angle2pxd( rad2deg.(ang), fov_angle, imsize)
end

function angle2pxd(angd, fov_angled=[62.61721188568244,35.793211268714096,71.6855447884958], imsize=(2160, 3840))
	angd ./ [fov_angled[1]/2 -fov_angled[2]/2] .* [imsize[2]/2 imsize[1]/2] .+ [imsize[2]/2 imsize[1]/2]
end

function plotImg(img, angs2, tdoa_slider)
	horizontal_angle = 62.61721188568244
	vertical_angle = 35.793211268714096
	# diagonal_angle = 71.6855447884958

	gr()
	plot(img)
	
	# ang = [ angs2[tdoa_slider,1]|>rad2deg, ]
	# ang = ((angs2[:,tdoa_slider] .|> rad2deg ) .+90 )./180 .*size(img)[[2,1]]
	ang = angle2px(angs2[tdoa_slider,:]')
	scatter!( [ang[1]], [ang[2]]  ) #([1000; 4000], [2000; 500])
	
	# ang = angle2pixel([angs2[:,tdoa_slider]], cameraCalibrationFile)
	# scatter!( [ang[1]], [ang[2]]  ) #([1000; 4000], [2000; 500])
	
	# plot!(title=
	# 	string(tdoa_slider)*" tdoa:"*string(tdoa[tdoa_slider,:]) *"\n"*
	# 	string(angs2[tdoa_slider,:] .|> rad2deg),
	# 	legend=false)
end



# (tdoa[oklist,1:3]./fs, rx_vect)
# gr()
# plot(img)
# ang = ((angs2[:,tdoa_slider] .|> rad2deg ) .+90 )./180 .*size(img)[[2,1]]
# scatter!( [ang[1]], [ang[2]]  ) #([1000; 4000], [2000; 500])

# ang = angle2pixel([angs2[:,tdoa_slider]], cameraCalibrationFile)
# scatter!( [ang[1]], [ang[2]]  ) #([1000; 4000], [2000; 500])
# plot!(title=string(angs2[:,tdoa_slider] .|> rad2deg), legend=false)

# angs2[:,tdoa_slider] .|> rad2deg

# ((angs2[:,tdoa_slider] .|> rad2deg ) .+90 )./180 .*size(img)[[2,1]]
# plot(angs2' .|> rad2deg)



function get_fps(file::AbstractString, streamno::Integer = 0)
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
    fps_strs = FFMPEG.exe(
        `-v 0 -of compact=p=0 -select_streams v:0 -show_entries stream=r_frame_rate $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug fps_strs
	try
		fps = split(fps_strs[1], '=')[2]
		if occursin("No such file or directory", fps)
			error("Could not find file $file")
		elseif occursin("N/A", fps)
			return nothing
		end
		return reduce(//, parse.(Int, split(fps,'/')) )
		# return round(reduce(/, parse.(Float64, split(fps,'/')) ), digits=3)
	catch err
		@debug(err)
		return NaN
	end
end

function get_framerate(file::AbstractString, streamno::Integer = 0, video_or_audio="v")
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
	if video_or_audio == "v"
		entries = "r_frame_rate"
	elseif video_or_audio == "a"
		entries = "sample_rate"
	end
	
    fps_strs = FFMPEG.exe(
        `-v 0 -of compact=p=0 -select_streams $video_or_audio:$streamno -show_entries stream="$entries" $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug fps_strs
	try
		fps = split(fps_strs[1], '=')[2]
		if occursin("No such file or directory", fps)
			error("Could not find file $file")
		elseif occursin("N/A", fps)
			return nothing
		end
		return reduce(//, parse.(Int, split(fps,'/')) )
		# return round(reduce(/, parse.(Float64, split(fps,'/')) ), digits=3)
	catch err
		@debug(err)
		return NaN
	end
end

function get_whatever(file::AbstractString, streamno::Integer = 0, video_or_audio="v"; entries_custom=nothing)
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
	if video_or_audio == "v"
		entries = "r_frame_rate"
	elseif video_or_audio == "a"
		entries = "sample_rate"
	end
	if !isnothing(entries)
		entries = entries_custom
	end
	
    fps_strs = FFMPEG.exe(
        `-v 0 -of compact=p=0 -select_streams $video_or_audio:$streamno -show_entries stream="$entries" $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug fps_strs
	try
		fps = split(fps_strs[1], '=')[2]
		if occursin("No such file or directory", fps)
			error("Could not find file $file")
		elseif occursin("N/A", fps)
			return nothing
		end
		if occursin('/', fps)
			return reduce(//, parse.(Int, split(fps,'/')) )
		else
			return parse.(Float64, fps) 
		end
		# return round(reduce(/, parse.(Float64, split(fps,'/')) ), digits=3)
	catch err
		@debug(err)
		return NaN
	end
end

"""
    get_number_frames(file [, streamno])
Query the the container `file` for the number of frames in video stream
`streamno` if applicable, instead returning `nothing` if the container does not
report the number of frames. Will not decode the video to count the number of
frames in a video.
"""
function get_number_frames(file::AbstractString, streamno::Integer = 0)
    streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
    frame_strs = FFMPEG.exe(
		`-v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets $file`, #-hide_banner
        # `-v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets -of csv=p=0 $file`,
        command = FFMPEG.ffprobe,
        collect = true,
    )
	@debug frame_strs
	frame_str = frame_strs[1]
	# num_frames = parse(Int, split(frame_str,'=')[end])
    if occursin("No such file or directory", frame_str)
        error("Could not find file $file")
    elseif occursin("N/A", frame_str)
        return NaN
    end
	
	try
		frame_str = frame_strs[2]
	    return parse(Int, split(frame_str,'=')[end])
	catch err
		@debug (err)
		return NaN
	end
end

function get_duration(file::AbstractString, streamno::Integer = 0)
	try
		streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))

		
		frame_strs = FFMPEG.exe(
			`-show_entries format=duration -v quiet -of csv="p=0" $file`,
			# `-v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets -of csv=p=0 $file`,
			command = FFMPEG.ffprobe,
			collect = true,
		)
		@debug frame_strs
		frame_str = frame_strs[1]
		# num_frames = parse(Int, split(frame_str,'=')[end])
		if occursin("No such file or directory", frame_str)
			error("Could not find file $file")
		elseif occursin("N/A", frame_str)
			@debug "manually calculate duration"
			return get_number_frames(file) / get_fps(file) |> Float64
		end
	
	# try
		frame_str = frame_str
	    return parse(Float64, split(frame_str,'=')[end])
	catch err
		@debug (err)
		return NaN
	end
end

# function get_duration_smart(fname; vidtype=r".mkv|.MP4|.avi|.mp4", autype=r".wav|.mat|.flac.mp3|.aac")
# 	occursin(vidtype, fname) && return get_number_frames(fname) / get_fps(fname) |> Float64
# 	occursin(autype, fname)  && return get_duration2(fname)
# end

# function get_duration2(vidfname)
# 	try
# 		get_number_frames(vidfname) / get_fps(vidfname) |> Float64
# 	catch err
# 		@debug(err)
# 		# # try
# 		# # 	get_duration2(vidfname)
# 		# # catch err
# 		# 	@debug(err)
# 			NaN
# 		# end
# 	end
# end


function get_videos_audiodata(vidfname)
	strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -f s16le -acodec pcm_s16le -loglevel error -`)# .|> Int16;
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# vid_audioFS =  get_framerate(vidfname, 0, "a")
	# length(strs)/vid_audioFS/2 #vid_auDur
	# get_duration(vidfname) #vid_Dur

	# # method 1, fast and creates array
	# vid_audiodata = Array{Int16}(undef, Int(length(strs)/2)) # 8bits to 16bits per frame
	# map!( x -> strs[2x[1]-1] + Int16(256)*strs[2x[1]] , vid_audiodata, 1:length(vid_audiodata)) # convert to 8bits little endian to Int16 merging each 2 bytes to 1 frame
	# vid_audiodata, get_framerate(vidfname, 0, "a")
	
	# # method 2, slow but easy to read, create array
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# method 3, fastest, doesnt create array, does the job
	vid_audiodata = reinterpret(Int16, strs)#, get_framerate(vidfname, 0, "a")
	if get_whatever(vidfname, 0, "a"; entries_custom="channels")>1
		vid_audiodata = reshape(vid_audiodata, get_whatever(vidfname, 0, "a"; entries_custom="channels")|>Int, :)'
	end
	vid_audiodata, get_framerate(vidfname, 0, "a")
end

function get_videos_audiodata_direct(vidfname)
	strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -f s16le -acodec pcm_s16le -loglevel error -`)# .|> Int16;
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# vid_audioFS =  get_framerate(vidfname, 0, "a")
	# length(strs)/vid_audioFS/2 #vid_auDur
	# get_duration(vidfname) #vid_Dur

	# # method 1, fast and creates array
	# vid_audiodata = Array{Int16}(undef, Int(length(strs)/2)) # 8bits to 16bits per frame
	# map!( x -> strs[2x[1]-1] + Int16(256)*strs[2x[1]] , vid_audiodata, 1:length(vid_audiodata)) # convert to 8bits little endian to Int16 merging each 2 bytes to 1 frame
	# vid_audiodata, get_framerate(vidfname, 0, "a")
	
	# # method 2, slow but easy to read, create array
	# ltoh.(reinterpret(Int16, strs)), get_framerate(vidfname, 0, "a")
	# method 3, fastest, doesnt create array, does the job
	reinterpret(Int16, strs), get_framerate(vidfname, 0, "a")
end


@info "DONE!"

# frame_strs = FFMPEG.exe(`-i $vidfname -f s16le -acodec pcm_s16le -loglevel error -`, command=FFMPEG.ffmpeg, collect=true,)
# strs = read(`ffmpeg -i $vidfname -f s16le -acodec pcm_s16le -loglevel error -`);

# strs = @ffmpeg_env read(`$ffmpeg -i $vidfname -f s16le -acodec pcm_s16le -loglevel error -`)# .|> Int16;
# vid_audioFS =  get_framerate(vidfname, 0, "a")
# length(strs)/vid_audioFS/2 #vid_auDur
# get_duration(vidfname) #vid_Dur

# vid_audiodata = Array{Int16}(undef, Int(length(strs)/2))
# map!( x -> strs[2x[1]-1] + Int16(256)*strs[2x[1]] , vid_audiodata, 1:length(vid_audiodata))
# # map!( x -> (2x[1]-1,2x[1]) , enumerate(vid_audiodata))
# # a = map( x -> strs[2x[1]-1] + strs[2x[1]] , enumerate(vid_audiodata))

# @btime vid_audiodata = Array{Int16}(undef, Int(length(strs)/2)); map!( x -> strs[2x[1]-1] + Int16(256)*strs[2x[1]] , vid_audiodata, 1:length(vid_audiodata))


# readdir(folname; join=true) |> skiphiddenfiles .|> [_ get_duration(_)]
# a = readdir("/Volumes/My Passport/Bahamas_2022/2022.06.25/0001"; join=true) |> skiphiddenfiles .|> get_duration

# folname = "/Volumes/My Passport/Bahamas_2022/2022.06.25/0002"
# flist = readdir(folname; join=true) |> skiphiddenfiles
# vidtype=r".mkv|.MP4|.avi|.mp4"; autype=r".wav|.mat|.flac.mp3"
# vids = filter( x -> occursin(vidtype, x), flist)
# aus = filter( x -> occursin(autype, x), flist)

# trigger_times = [basename.(vids) findVidAudioBlip.(vids; plot_window_inS=nothing) basename.(aus) findAudioBlip.(aus; plot_window_inS=nothing)]