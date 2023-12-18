using Plots
using ProgressMeter
import GLMakie
using Pipe: @pipe
# include("map.jl")
using Base.Threads

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
	a = plot(d)
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


function plotImg(img, angs2, tdoa_slider; angle2px=(x,y)->(x,y), fov_angle=fov_angle, plotfunc=gr, markersize=10, kwargs...)
	# horizontal_angle = 62.61721188568244
	# vertical_angle = 35.793211268714096
	# diagonal_angle = 71.6855447884958

	plotfunc()
	Plots.plot(img)
	
	# ang = [ angs2[tdoa_slider,1]|>rad2deg, ]
	# ang = ((angs2[:,tdoa_slider] .|> rad2deg ) .+90 )./180 .*size(img)[[2,1]]
	px = angle2px(angs2[tdoa_slider,:]', fov_angle)
	@debug tdoa_slider, px#, angs2[tdoa_slider,:]'
	Plots.scatter!( [px[1]], [px[2]]; markersize=markersize, kwargs...) #([1000; 4000], [2000; 500])
	# @info [tdoa_slider ang]
	
	# ang = angle2pixel([angs2[:,tdoa_slider]], cameraCalibrationFile)
	# scatter!( [ang[1]], [ang[2]]  ) #([1000; 4000], [2000; 500])
	
	# plot!(title=
	# 	string(tdoa_slider)*" tdoa:"*string(tdoa[tdoa_slider,:]) *"\n"*
	# 	string(angs2[tdoa_slider,:] .|> rad2deg),
	# 	legend=false)
end

function plot_one_event(event_plots_dir, vidfname, pind_good_inS, data, pind_threshold_indices, tdoa, window, ang;
	func2=plotTDOA_raw, plotfunc=gr, i=1, kwargs...)
	vid = VideoIO.openvideo(vidfname)
    imsize = raw_frame_size(vid)
	p = plotImg( readImage(vid, pind_good_inS[i]), ang, i; angle2px=angle2px, plotfunc=plotfunc, kwargs...)
        p = Plots.plot!(title=
            string(i)*"_"*string(pind_good_inS[i])*"s tdoa:"*string(tdoa[i,:]) *"\n"*
            string(@pipe ang[i,:] .|> rad2deg .|> round(_; digits=3))*"\n"*
            string(angle2px(ang[i,:]') .|> round .|> Int),
            legend=false)
        p = Plots.xlims!( 1, imsize[1])
        p = Plots.ylims!( 1, imsize[2])

        # savefig(p, joinpath(event_plots_dir, string(i)*".png"))

        p2 = func2(data, i, pind_threshold_indices, tdoa; window=window, func=plotfunc)
        # p2 = plot!(title=
        #     string(i)*"_"*string(pind_good_inS[i])*"s tdoa:"*string(tdoa[i,:]),# *"\n"*
        #     # string(ang[i,:] .|> rad2deg),
        #     legend=false)
		# @info p
		# @info p2
        pnew = Plots.plot(p, p2, layout=@layout [a b])
		# @info pnew
		return pnew
end


# event_plots_dir = joinpath(res_dir, basename(vidfname)*"_clicks"*"_"*"_t"*string(thresh)*"_d"*string(dist))
# plot_all_clicks(event_plots_dir*"_tdoa-raw", vidfname, pind_good_inS, data, pind[threshold_indices], tdoa_raw, window, ang)
# plot_all_clicks(joinpath(res_dir, splitext(basename(aufname))[1] ), vidfname, detected_tonal_inS, data_filt, (detected_tonal_inS.*fs).|>round.|>Int, hcat(tdoa...)', 1:1000, vcat(pxs...) .|> deg2rad)
function plot_all_clicks(event_plots_dir, vidfname, pind_good_inS, data, pind_threshold_indices, tdoa, window, ang, ppeak=ones(length(pind_good_inS));
	func2=plotTDOA_raw, plotfunc=gr,
	plotsize=(600,400), kwargs...)

	vid = VideoIO.openvideo(vidfname)
    imsize = raw_frame_size(vid)
    mkpath(event_plots_dir)
	plotfunc()
	
	ppeak_norm = ppeak ./ maximum(ppeak)
	# if Threads.nthreads()>1
	# 	close(vid)
	# 	vid = vidfname
	# end

	@info "Writing each detection to image..." 
    @time @showprogress "Writing each detection to image..." for i in eachindex(pind_good_inS)
	# Threads.@threads for i in eachindex(pind_good_inS)
	# @time @threads for i in eachindex(pind_good_inS)
        # @show angle2px(ang[i,:]') .|> round .|> Int |> string
        p = plotImg( readImage(vid, pind_good_inS[i]), ang, i; angle2px=angle2px, plotfunc=plotfunc, alpha=ppeak_norm[i], kwargs...)
        p = Plots.plot!(title=
            string(i)*"_"*string(pind_good_inS[i])*"s tdoa:"*string(tdoa[i,:]) *"\n"*
            string(@pipe ang[i,:] .|> rad2deg .|> round(_; digits=3))*"\n"*
            string(angle2px(ang[i,:]') .|> round .|> Int),
            legend=false)
        p = Plots.xlims!( 1, imsize[1])
        p = Plots.ylims!( 1, imsize[2])

        # savefig(p, joinpath(event_plots_dir, string(i)*".png"))

        p2 = func2(data, i, pind_threshold_indices, tdoa; window=window, func=plotfunc)
        # p2 = plot!(title=
        #     string(i)*"_"*string(pind_good_inS[i])*"s tdoa:"*string(tdoa[i,:]),# *"\n"*
        #     # string(ang[i,:] .|> rad2deg),
        #     legend=false)
		# @info p
		# @info p2
        pnew = Plots.plot(p, p2, layout=@layout [a b]; size=plotsize)
		if false #plotfunc == plotlyjs
			open(joinpath(event_plots_dir, string(i)*".html"), "w") do io
				PlotlyBase.to_html(io, pnew)
			end
		else
	        Plots.savefig(pnew, joinpath(event_plots_dir, string(i)*".png"))
		end
    end
	# if vid isa String
	# 	return
	# end
	close(vid)
end



# function plot_vidNaudio(img, data; resolution=(1920,1080))

    # indata = @view data[1:400000,1]
    
    # resolution=(1920,1080)
    # fig = Figure(;resolution=result_resolution)
    # ax = Axis(fig[9,16])
function plot_summary!(fig, indata, fs, img; resolution=result_resolution, nfft=4000, freq_crop=1:401, title="")
    # fig = Figure(;resolution=resolution)
    empty!(fig);
    aspect_ratio = reduce(//, resolution);
    # ax = Axis(fig[denominator(aspect_ratio), numerator(aspect_ratio)])
    # ga = fig[denominator(aspect_ratio), numerator(aspect_ratio)] = GridLayout();

	@debug (size(indata), typeof(indata), fs, fig);
    lines!(Axis(fig[end,1:end-1]), (1:length(indata))./fs, indata);
    tightlimits!(Axis(fig[end,1:end-1])) #GLMakie.xlims!(0,length(indata)/fs);
    # GLMakie.ylims!(-500, 500)
    # ylim
    # resize_to_layout!(fig)
    
    spec = spectrogram(indata, nfft; nfft=nfft,fs=fs);
    a = GLMakie.heatmap( fig[1:end,end], @view(spec.freq[freq_crop]), spec.time, pow2db.(@view spec.power[freq_crop,:]), xticklabelrotation=pi/2);
    a.axis.xticks = (0:20_000:40_000, string.(Int.(0:20_000:40_000)./1000) );

    # a = GLMakie.heatmap(fig[1:end-1,1:end-1], @view(img'[:,end:-1:1]) )
	
	b = image(fig[1:end-1,1:end-1], @view(img'[:,end:-1:1]));
	b.axis.title = title;#"time: frame: click: whistle:"
	
	# GLMakie.vlines!([0,0.04])
    # a.axis.xticks = nothing
    return Makie.colorbuffer(fig.scene);
end

plot_summary!(args...; kwargs...) = plot_summary!(args;kwargs)

function plot_summary_plots!(fig, indata, fs, img; resolution=result_resolution, nfft=4000, freq_crop=1:401, title="")
	spec = spectrogram(indata, nfft; nfft=nfft,fs=fs);
	p1 = Plots.plot(img; ticks = false)#, xlims=(0,resolution[1]), ylims=(0,resolution[2]));
	p2 = Plots.plot((1:length(indata))./fs, indata; label=nothing, xlims=(0,length(indata)/fs));
	p3 = Plots.heatmap(@view(spec.freq[freq_crop]), spec.time, pow2db.(@view spec.power[freq_crop,:])'; colorbar=nothing,
		xticks=(0:20_000:40_000, string.(Int.(0:20_000:40_000)./1000)));
	Plots.plot(p1,p2,p3, size=result_resolution, layout =
	Plots.@layout[ [a{0.9h}; b] c{0.058823529411764705w}] ) #@layout[ [a{9/10h};b]{16/17w} c]
	Plots.title!(title)
end

function plot_summary_plots_img!(fig, indata, fs, img; resolution=result_resolution, nfft=4000, freq_crop=1:401, title="")
	# @info (args, kwargs)
	p = plot_summary_plots!(fig, indata, fs, img; resolution=resolution, nfft=nfft, freq_crop=freq_crop, title=title)
	iobuf = IOBuffer()
	png(p, iobuf)
	load(iobuf)
end

# function plot_summary!( args; kwargs...)
#     plot_summary!( args...; resolution=resolution, nfft=nfft, freq_crop=nfft)
# end
# function plot_nextVideo!(vid, data)
function plot_nextVideo!(img, counter, extra_arg) 
    # pind_vidframes_list, p_pixels_list, colour, ptsize = extra_arg
    extra_arg1, extra_arg2 = extra_arg
    overlay_points!(img, counter, extra_arg1)
    data, fs = extra_arg2
    
    fig = Figure(;resolution=resolution)
    plot_summary!(fig, indata, fs, img)
end

function get_pixel(img)
	fig = Figure()
	ax1 = Axis(fig[1,1])
	image!(ax1,img)
	pts = []

	register_interaction!(ax1, :my_interaction) do event::GLMakie.MouseEvent, axis
		if event.type === MouseEventTypes.leftclick
			println("You clicked on the axis at datapos $(event.data)")
			push!(pts, event.data)
		end
	end
end

function get_pixelLoc(pts, fig, vidfname, pind_good_inS, data, pind_threshold_indices, tdoa, window, ang)
    vid = VideoIO.openvideo(vidfname)
    imsize = raw_frame_size(vid)
    # mkpath(event_plots_dir)


	# fig = Figure()
	ax1 = Axis(fig[1,1])
	image!(ax1,img)
	# pts = []
	i = 1

	register_interaction!(ax1, :my_interaction) do event::GLMakie.MouseEvent, axis
		if event.type === MouseEventTypes.leftclick
			println("You clicked on the axis at datapos $(event.data)")
			push!(pts, event.data)
			
			i += 1
			image!(ax1, readImage(vid, pind_good_inS[i]))
			ax1.title = string(i)
		end
	end
	# while i <= length(d["pind_good_inS"])
	# 	sleep(5)
	# 	# print(int2str(i))
	# end
	println("done checking")

	return pts
end

function plot_ang(res_new, ang, aufname="", res_dir=nothing; type="Amplitude", ylabel="Azimuth(°)", label="clicks")
	maxi = res_new.ppeak |> maximum
	
	if lowercase(type)==lowercase("Amplitude")
		transparency = res_new.ppeak./maxi
	elseif lowercase(type)==lowercase("Power")
		transparency = (res_new.ppeak.^2)./(maxi^2)
	else
		@error "Type unknown, use default Amplitutde"
		type="Amplitude"
	end

	p = scatter(res_new.pind_good_inS, ang.|>rad2deg; alpha=transparency,
    xlabel="Time(s)", ylabel=ylabel, label=label,
    title="Detection "*type*"(Transparency)",
	markershape=:auto)

	if !isnothing(res_dir)
		mkpath(res_dir)
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_$label"*type*".html"))
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_$label"*type*".png"))
	end
	p
end
function plot_ang!(res_new, ang, aufname="", res_dir=nothing; type="Amplitude")
	maxi = res_new.ppeak |> maximum
	
	if lowercase(type)==lowercase("Amplitude")
		transparency = res_new.ppeak./maxi
	elseif lowercase(type)==lowercase("Power")
		transparency = (res_new.ppeak.^2)./(maxi^2)
	else
		@error "Type unknown, use default Amplitutde"
		type="Amplitude"
	end

	p = scatter!(res_new.pind_good_inS, ang.|>rad2deg; alpha=transparency,
    xlabel="Time(s)", ylabel="Azimuth(°)", label="clicks",
    title="Detection "*type*"(Transparency)",
	markershape=:auto)

	if !isnothing(res_dir)
		mkpath(res_dir)
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_click"*type*".html"))
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_click"*type*".png"))
	end
	p
end

function plot_azimuth(res_new, ang, aufname="", res_dir=nothing; type="Amplitude")
	maxi = res_new.ppeak |> maximum
	
	if lowercase(type)==lowercase("Amplitude")
		transparency = res_new.ppeak./maxi
	elseif lowercase(type)==lowercase("Power")
		transparency = (res_new.ppeak.^2)./(maxi^2)
	else
		@error "Type unknown, use default Amplitutde"
		type="Amplitude"
	end

	p = scatter(res_new.pind_good_inS, ang[:,2].|>rad2deg; alpha=transparency,
    xlabel="Time(s)", ylabel="Azimuth(°)", label="clicks",
    title="Detection "*type*"(Transparency)")

	if !isnothing(res_dir)
		mkpath(res_dir)
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_click"*type*".html"))
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_click"*type*".png"))
	end
	p
end

function plot_elevation(res_new, ang, aufname="", res_dir=nothing; type="Amplitude")
	maxi = res_new.ppeak |> maximum
	
	if lowercase(type)==lowercase("Amplitude")
		transparency = res_new.ppeak./maxi
	elseif lowercase(type)==lowercase("Power")
		transparency = (res_new.ppeak.^2)./(maxi^2)
	else
		@error "Type unknown, use default Amplitutde"
		type="Amplitude"
	end

	p=scatter!(res_new.pind_good_inS, ang[:,1].|>rad2deg; alpha=transparency,
    xlabel="Time(s)", ylabel="Azimuth(°)", label="clicks",
    title="Detection "*type*"(Transparency)")

	if !isnothing(res_dir)
		mkpath(res_dir)
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_click"*type*".html"))
		savefig(joinpath(res_dir, splitext(basename(aufname))[1] *"_click"*type*".png"))
	end
	p
end

# function func_splatkwargs3(args...; func=x->x, kwargs...)
# 	@info size(args)
# 	p = plot!()
# 	kwarg=[]
# 	kwarg_splat=[]
# 	for key in keys(kwargs)
# 		if kwargs[key] isa Array
# 			push!(kwarg, key=>kwargs[key])
# 		else
# 			push!(kwarg_splat, key=>kwargs[key])
# 		end
# 	end
# 	for i in 1:size(args[1],1) #eachindex(args[1])
# 		@info i
# 		@info func
# 		@info [x[1]=>x[2][i] for x in kwarg]
# 		# @info kwarg_splat

# 		p=func(args...[i,:]; [x[1]=>x[2][i] for x in kwarg]..., kwarg_splat...)
# 	end
# 	p
# end

function plot_linesfrom3(origin, pts; kwargs...)
	p = plot!()
	kwarg=[]
	kwarg_splat=[]
	for key in keys(kwargs)
		if kwargs[key] isa Array
			push!(kwarg, key=>kwargs[key])
		else
			push!(kwarg_splat, key=>kwargs[key])
		end
	end
	for i in eachindex(pts)
		p=plot!([origin; pts[i]]; [x[1]=>x[2][i] for x in kwarg]..., kwarg_splat...)
	end
	p
end


function plot_linesfrom(origin, pts, alphas=repeat([1],size(pts,1)); kwargs...)
	p = plot!()
	for i in eachindex(pts)
		p=plot!([origin; pts[i]]; alpha=alphas[i], kwargs...)
	end
	p
end

function plot_spectro(res_new) 
	heatmap(res_new.res_tonal.time_index, res_new.res_tonal.freqss, res_new.res_tonal.mag_max; size=(1500,800)); title!(res_new.res_tonal.percent_quiet|>string)
end
# save_fig(p, res_new, res_dir) = savefig(p, joinpath(res_dir,"tonal_bw"*string(res_new.res_tonal.freqss|>extrema)*"_pq"*(res_new.res_tonal.percent_quiet|>string)*".html"))
# "/Users/abel/Documents/data_res/aspod/Hawaii_2022-09/punnet3/tonal_bw(400.0, 63900.0)_pq0.html"

function plot_map_detections_bearing(node_loc,radius, res_new, ang2, ang2_corrected)
	p=plot_nodes(node_loc; radius=radius)

    #    alpha=(res_new.ppeak ./ maximum(res_new.ppeak)), 
    #    color=palette_continuous(res_new.pind_good_inS;regular=false) )

	# plot per click train
	colmap = palette_continuous(res_new.train_start;regular=false)
	for i ∈ 1:length(res_new.train_start_ind)
		if i < length(res_new.train_start_ind) 
			win = res_new.train_start_ind[i]:res_new.train_start_ind[i+1]
		else
			win = res_new.train_start_ind[i]:length(res_new.pind_good_inS)
		end
		p=plot_nodes(node_loc; radius=radius)
		
		plot!([repeat([node_loc[2]], 1,size(ang2_corrected[win],1)) ; pts_from_location.(Ref(node_loc[2]), ang2_corrected[win]', Ref(1000))];
		color=colmap[i], alpha=0.05); 
		# plot_linesfrom3(node_loc[2], pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], Ref(1000)); 
		# alpha=res_new.ppeak[win] ./ maximum(res_new.ppeak) .* 0.25, 
		# color=colmap[i])#palette_continuous(res_new.pind_good_inS;regular=false))
		scatter!(pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], 600 .- (res_new.pind_good_inS[win])), 
		alpha=res_new.ppeak[win] ./ maximum(res_new.ppeak) .* 0.25, 
		color=palette_continuous(res_new.pind_good_inS[win];regular=false))
		# color=colmap[i])

		p=title!(string(i) *"__"* string(res_new.pind_good_inS[[win[1]; win[end]]]) *"s")
		display(p)
		# sleep(1)
	end
	return p
end

function plot_map_detections_bearing_gif(node_loc,radius, res_new, ang2, ang2_corrected;
		overlay=false, fps=25)
	

    #    alpha=(res_new.ppeak ./ maximum(res_new.ppeak)), 
    #    color=palette_continuous(res_new.pind_good_inS;regular=false) )

	# plot per click train
	gr()
	colmap = palette_continuous(0:1/fps:res_new.pind_good_inS[end];regular=false)
	i=1
	p=plot_nodes(node_loc; radius=radius)
	anim = @animate for t ∈ 0:1/fps:res_new.pind_good_inS[end] #1:length(res_new.train_start_ind)
		win = findall(x -> x>t && x<t+1/fps, res_new.pind_good_inS)
		# if i < length(res_new.train_start_ind) 
		# 	win = res_new.train_start_ind[i]:res_new.train_start_ind[i+1]
		# else
		# 	win = res_new.train_start_ind[i]:length(res_new.pind_good_inS)
		# end
		overlay || plot_nodes(node_loc; radius=radius);
		# plot!([repeat([node_loc[2]], 1,size(ang2_corrected[win],1)) ; pts_from_location.(Ref(node_loc[2]), ang2_corrected[win]', Ref(1000))];
		# color=colmap[i], alpha=0.05); 
		plot_linesfrom3(node_loc[2], pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], Ref(1000)); 
		alpha=res_new.ppeak[win] ./ maximum(res_new.ppeak), 
		color=colmap[i])#palette_continuous(res_new.pind_good_inS;regular=false))
		scatter!(pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], 600 .- (res_new.pind_good_inS[win])), 
		alpha=res_new.ppeak[win] ./ maximum(res_new.ppeak), 
		# color=palette_continuous(res_new.pind_good_inS[win];regular=false))
		color=colmap[i])
		i += 1
		p=title!(string(t) *"s")
		# display(p)
		# sleep(1)
		# @info(t)
		mod(t,1) == 0 && @info( string(t) *"/"* string(res_new.pind_good_inS[end]))
	end
	return anim
end

function plot_map_detections_bearing_gif_hold(node_loc,radius, res_new, ang2, ang2_corrected;
	overlay=false, fps=25, t_fade_lag=2)


	#    alpha=(res_new.ppeak ./ maximum(res_new.ppeak)), 
	#    color=palette_continuous(res_new.pind_good_inS;regular=false) )

	# plot per click train
	gr()
	colmap = palette_continuous(0:1/fps:res_new.pind_good_inS[end];regular=false)
	i=1
	p=plot_nodes(node_loc; radius=radius)
	anim = @animate for t ∈ 0:1/fps:res_new.pind_good_inS[end] #1:length(res_new.train_start_ind)
		# win = findall(x -> x>t && x<t+1/fps, res_new.pind_good_inS)
		win = findall(x -> x>t-t_fade_lag && x<t+1/fps, res_new.pind_good_inS)
		# if i < length(res_new.train_start_ind) 
		# 	win = res_new.train_start_ind[i]:res_new.train_start_ind[i+1]
		# else
		# 	win = res_new.train_start_ind[i]:length(res_new.pind_good_inS)
		# end
		overlay || plot_nodes(node_loc; radius=radius);
		# plot!([repeat([node_loc[2]], 1,size(ang2_corrected[win],1)) ; pts_from_location.(Ref(node_loc[2]), ang2_corrected[win]', Ref(1000))];
		# color=colmap[i], alpha=0.05); 
		if isempty(win)
			time_normalized = []
		else
			time_normalized = res_new.pind_good_inS[win] .- minimum(res_new.pind_good_inS[win])
			time_normalized = time_normalized ./ maximum(time_normalized)
		# end

			plot_linesfrom3(node_loc[2], pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], Ref(radius*3)); 
			alpha=time_normalized , #res_new.ppeak[win] ./ maximum(res_new.ppeak), 
			color=palette_continuous(res_new.ppeak[win] ./ maximum(res_new.ppeak);regular=false, colortype=:red))
			# color=colmap[i])#
			# color=palette_continuous(res_new.pind_good_inS[win];regular=false))

			# @debug pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], 300 .+ (res_new.pind_good_inS[win]))
			# @debug ang2_corrected[win]
			# @debug res_new.ppeak[win] ./ maximum(res_new.ppeak) .* 100
			plot_arrow!.(pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], 300 .+ (res_new.pind_good_inS[win])), 
				ang2_corrected[win], res_new.ppeak[win] ./ maximum(res_new.ppeak) .* 100; 
				color=:green
				)

			# scatter!(pts_from_location.(Ref(node_loc[2]), ang2_corrected[win], 600 .- (res_new.pind_good_inS[win])), 
			# alpha=res_new.ppeak[win] ./ maximum(res_new.ppeak[win]), 
			# color=palette_continuous(res_new.pind_good_inS[win];regular=false))
			# # color=colmap[i])
		end

		i += 1
		p=title!(string(round(t;digits=3)) *"s")
		# display(p)
		# sleep(1)
		# @info(t)
		mod(t,1) == 0 && @info( string(t) *"/"* string(res_new.pind_good_inS[end]))
	end
	return anim
end

# fname = "/Volumes/One Touch/res/Hawaii_2022-09/punnet_yellow/4/2022-09-25/counts.csv"
# using CSV, TimeZones, DataFrames
function plot_detection_summary(fname; res_dir=nothing, filetype=".png", plotly_flag=false)
	
	detections = CSV.read(fname, DataFrame)
	detections.datetime = ZonedDateTime.(detections.datetime .|> String)

	#~ plotly
	if plotly_flag
		p = PlotlyJS.Plot(detections.datetime .|> DateTime, [detections.num_noise detections.num_tonal],
		# labels=Dict("num_noise"=>"noise", "num_tonal"=>"tonal"),
		PlotlyJS.Layout(
		title="Detections ("* string(detections.datetime[1] |> Date) *")",
		xaxis_title="DateTime",
		yaxis_title="Number of Detections"
		# legend_title="Legend Title"
		)
		# ; labels=["noise";"tonal"]
		)

		if !isnothing(res_dir)
			open(joinpath(res_dir, dirname(fname)|>basename) * filetype, "w") do io
				PlotlyBase.to_html(io, p)
			end
		end
		return p
	end


	plot(detections.datetime, [detections.num_tonal]; 
		labels="number of tonals",
	# plot(detections.datetime, [detections.num_tonal detections.num_noise]; 
	# 	labels=["number of tonals" "number of noise"],
		title="Detections ("* string(detections.datetime[1] |> Date) *")",
		xlabel="DateTime",
		ylabel="Number of Detections"
		)

	if !isnothing(res_dir)
		savefig(joinpath(res_dir, dirname(fname)|>basename) * filetype)
	end
end


#~ others
# plot_detection_summary("/Volumes/One Touch/res/Hawaii_2022-09/punnet_yellow/4/2022-09-17/counts.csv"; res_dir="/Volumes/One Touch/res/Hawaii_2022-09/punnet_yellow/4/summary")
# plot_detection_summary("/Volumes/One Touch/res/Hawaii_2022-09/punnet_yellow/4/2022-09-27/counts.csv"; res_dir="/Volumes/One Touch/res/Hawaii_2022-09/punnet_yellow/4/summary", filetype=".html", plotly_flag=true)

function detectionsfiles2plot(fname; res_dir=nothing)
    df = CSV.read(fname, DataFrame)
    df.datetime = ZonedDateTime.(String.(df.datetime)) .|> DateTime
    # plot(df.datetime, [df.num_tonal df.num_noise]; label=["tonal" "noise"])

    p = PlotlyJS.plot([
        PlotlyJS.scatter(df, x=:datetime, y=:num_tonal; name="tonal"),
        PlotlyJS.scatter(df, x=:datetime, y=:num_noise; name="noise")
        ],
        PlotlyJS.Layout(
            title="Detections ("* string(df.datetime[1] |> Date) *")",
            xaxis_title="DateTime",
            yaxis_title="Number of Detections"
            # label=["noise" "tonal"]
            # legend_title="Legend Title"
        )
    )

    mkpath(res_dir)
    open(joinpath(res_dir,basename(dirname(fname))*".html"), "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
    # open(joinpath(res_dir,basename(dirname(fname))*".png"), "w") do io
    #     PlotlyBase.to_image(io, p.plot)
    # end
    PlotlyJS.savefig(p.plot, joinpath(res_dir,basename(dirname(fname))*".png"))

end

# 	a=PlotlyJS.plot(x=df.datetime .|> DateTime, y=df.num_noise, name="noise")
# 	b=PlotlyJS.scatter(x=df.datetime .|> DateTime, y=df.num_tonal, name="tonal")
# 	plot([a,b])

# 	PlotlyJS.Plot(detections.datetime .|> DateTime, [detections.num_noise detections.num_tonal],
# 	labels=Dict("num_noise"=>"noise", "num_tonal"=>"tonal"),
# 	PlotlyJS.Layout(
#     title="Plot Title",
#     xaxis_title="X Axis Title",
#     yaxis_title="Y Axis Title",
#     legend_title="Legend Title")
# 	# ; labels=["noise";"tonal"]
# 	)

# 	# PlotlyJS.add_vrect!(fig, "2022-09-27T07:07:26.201-1000", "2022-09-27T17:07:26.201-1000", fillcolor="LightSalmon", opacity=0.5,
#         #    layer="below", line_width=0)

# 	PlotlyJS.Plot(detection, df.datetime .|> DateTime, [df.num_noise df.num_tonal],
# 	labels=Dict("num_noise"=>"noise", "num_tonal"=>"tonal"),
# 	PlotlyJS.Layout(
#     title="Plot Title",
#     xaxis_title="X Axis Title",
#     yaxis_title="Y Axis Title",
#     legend_title="Legend Title")
# 	# ; labels=["noise";"tonal"]
# 	)

