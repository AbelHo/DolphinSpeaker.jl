import Plots
using ProgressMeter
using Pipe: @pipe
# include("map.jl")
using Base.Threads
# import GLMakie
using PlotlyBase

"""
	quiver_error_threshold(x, y, px, py; q=0.9, vid_width=nothing, vid_height=nothing, kwargs...)

Compute Euclidean distances between annotated positions (`x`,`y`) and estimated
pixel positions (`px`,`py`), compute the `q` quantile threshold of those
distances, select points with distance <= threshold and return a quiver plot of
the error vectors for those points.

Returns a NamedTuple with fields `:plot`, `:inds`, and `:threshold`.
"""
function quiver_error_threshold(x::AbstractVector, y::AbstractVector,
		px::AbstractVector, py::AbstractVector; q::Real=0.9,
		width=nothing, height=nothing, flag_display::Bool=true, kwargs...)

	# ensure vectors have same length
	n = length(x)
	@assert length(y) == n && length(px) == n && length(py) == n "Input vectors must have the same length"

	# Euclidean distances per point
	dist = sqrt.((x .- px).^2 .+ (y .- py).^2)

	# quantile threshold
	threshold = Statistics.quantile(dist, q)
	@info "quiver_error_threshold: n=$(n), threshold=$(threshold)"

	inds = findall(dist .<= threshold)

	# build quiver plot for selected indices
	dx = x[inds] .- px[inds]
	dy = y[inds] .- py[inds]

	p = Plots.scatter(px[inds], py[inds]; aspect_ratio=:equal)
	Plots.scatter!(p, x[inds], y[inds])
	
	Plots.quiver!(p, px[inds], py[inds], quiver=(dx, dy); aspect_ratio=:equal, color=:blue,
		label = "error vectors ($(Int(round(q*100)))% quantile)",
		xlabel = "Pixel px", ylabel = "Pixel py",
		title = "Annotated vs Estimated Positions ($(Int(round(q*100)))% quantile)", kwargs...)

	# apply optional limits if provided
	if width !== nothing
		try
			Plots.xlims!(p, 0, width)
		catch _
		end
	end
	if height !== nothing
		try
			Plots.ylims!(p, 0, height)
		catch _
		end
	end
	flag_display && Plots.display(p)
	return (plot = p, inds = inds, threshold = threshold)
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
"""
	Draw vertical lines at positions in `x` with the default color changing scheme.
	
	# Arguments
	- `x`: A vector of x-coordinates where vertical lines should be drawn.
	- `kwargs...`: Additional keyword arguments to customize the appearance of the lines.
	
	# Returns
	- The modified plot object with vertical lines added.
"""
function vlinec!(x; kwargs...)
	for i in eachindex(x)
		vline!([x[i]]; color=palette(:default)[i], kwargs...)
	end
	return plot!()
end

"""
	Draw horizontal lines at positions in `y` with the default color changing scheme.
	
	# Arguments
	- `y`: A vector of y-coordinates where horizontal lines should be drawn.
	- `kwargs...`: Additional keyword arguments to customize the appearance of the lines.
	
	# Returns
	- The modified plot object with horizontal lines added.
"""
function hlinec!(y; kwargs...)
	for i in eachindex(y)
		hline!([y[i]]; color=palette(:default)[i], kwargs...)
	end
	return plot!()
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
        pnew = Plots.plot(p, p2, layout=Plots.@layout [a b])
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
        pnew = Plots.plot(p, p2, layout=Plots.@layout [a b]; size=plotsize)
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

# function get_pixel(img)
# 	fig = Figure()
# 	ax1 = Axis(fig[1,1])
# 	image!(ax1,img)
# 	pts = []

# 	register_interaction!(ax1, :my_interaction) do event::GLMakie.MouseEvent, axis
# 		if event.type === MouseEventTypes.leftclick
# 			println("You clicked on the axis at datapos $(event.data)")
# 			push!(pts, event.data)
# 		end
# 	end
# 	return pts
# end

# function get_pixelLoc(pts, fig, vidfname, pind_good_inS)#, data, pind_threshold_indices, tdoa, window, ang)
#     vid = VideoIO.openvideo(vidfname)
#     imsize = raw_frame_size(vid)
#     # mkpath(event_plots_dir)


# 	# fig = Figure()
# 	ax1 = Axis(fig[1,1])
# 	image!(ax1,img)
# 	# pts = []
# 	i = 1

# 	register_interaction!(ax1, :my_interaction) do event::GLMakie.MouseEvent, axis
# 		if event.type === MouseEventTypes.leftclick
# 			println("You clicked on the axis at datapos $(event.data)")
# 			push!(pts, event.data)
			
# 			i += 1
# 			@info i
# 			image!(ax1, readImage(vid, pind_good_inS[i]))
# 			ax1.title = string(i)
# 		end
# 	end
# 	while i <= length(pind_good_inS)
# 		sleep(1)
# 		# print(int2str(i))
# 	end
# 	println("done checking")

# 	return pts
# end

# function get_pixelLoc2(vidfname, timestamps)
#     # Open the video
#     vid = VideoIO.openvideo(vidfname)

#     # Initialize an array to store the clicked coordinates
#     pts = []

#     # Create a figure and axis
#     fig = Figure()
#     ax1 = Axis(fig[1,1])
# 	image(@view(img[end:-1:1, :])')

#     # Register a mouse click interaction
#     register_interaction!(ax1, :my_interaction) do event::GLMakie.MouseEvent, axis
#         if event.type === MouseEventTypes.leftclick
#             println("You clicked on the axis at datapos $(event.data)")
# 			global i
# 			@info i
# 			image!(ax1, readImage(vid, timestamps[i]))
# 			ax1.title = string(i)

# 			i += 1
#             push!(pts, event.data)
#         end
#     end

#     # Loop over the timestamp
# 	i = 1; pts = [];
#     while i < length(timestamps) i in 1:length(timestamps)
#         # Read the frame at the current timestamp
#         img = readImage(vid, timestamps[i])

#         # Display the frame
#         # image!(ax1, img); #
# 		image(ax1,@view(img[end:-1:1, :])')
#         ax1.title = string(i)

#         # Wait for a click
#         # while length(pts) < i
#             sleep(0.1)
#         # end
#     end

#     return pts
# end


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
	anim = Plots.@animate for t ∈ 0:1/fps:res_new.pind_good_inS[end] #1:length(res_new.train_start_ind)
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
	anim = Plots.@animate for t ∈ 0:1/fps:res_new.pind_good_inS[end] #1:length(res_new.train_start_ind)
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
	try
		detections.datetime = ZonedDateTime.(detections.datetime .|> string)
	catch error
		@error "no zoned date time"
	end

	mkpath(res_dir)
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

function detectionsfiles2plot(fname; res_dir=nothing, plottype=PlotlyJS.bar)
    df = CSV.read(fname, DataFrame)
    try
		df.datetime = ZonedDateTime.(String.(df.datetime)) .|> DateTime
	catch error
		@error "no zoned date time"
	end
    # plot(df.datetime, [df.num_tonal df.num_noise]; label=["tonal" "noise"])

    p = PlotlyJS.plot([
        plottype(df, x=:datetime, y=:num_tonal; name="tonal"),
        plottype(df, x=:datetime, y=:num_noise; name="noise"),
        plottype(df, x=:datetime, y=:num_impulsetrain; name="click_train"),
        plottype(df, x=:datetime, y=:num_impulseINtrain; name="click")
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
	return p
end

function detectionsfiles2plot2(fname; res_dir=nothing, plottype=PlotlyJS.bar, 
	add_daynight=true, latitude=21.3, longitude=-157.8, # Hawaii coordinates by default
	trace_list = [:num_tonal=>"tonal", :num_noise=>"noise", :num_impulsetrain=>"click_train", :num_impulseINtrain=>"click",
	:count=>"count"]
	)

	dts_colname = :datetime


	if typeof(fname) == DataFrame
		df = fname
		dts_colname = Symbol( hasproperty(df, dts_colname) ? dts_colname : hasproperty(df, Symbol("timestamp")) ? "timestamp" : hasproperty(df, Symbol("time")) ? "time" : hasproperty(df, Symbol("date_time")) ? "date_time" : "dts" )
		fname = Dates.format(df[1,dts_colname], "yyyy-mm-dd_HHMMSS") * "/"
		@info "autoname: $fname"
    else
        df = CSV.read(fname, DataFrame)
		dts_colname = Symbol( hasproperty(df, dts_colname) ? dts_colname : hasproperty(df, Symbol("timestamp")) ? "timestamp" : hasproperty(df, Symbol("time")) ? "time" : hasproperty(df, Symbol("date_time")) ? "date_time" : "dts" )
	end

    try
        df[!, dts_colname] = ZonedDateTime.(String.(df[!,dts_colname])) .|> DateTime
    catch error
        @error "no zoned date time"
    end
    @debug "after zoned dt"
    # Create the base plot
    layout = PlotlyJS.Layout(
        title="Detections ("* string(df[1, dts_colname] |> Date) *")",
        xaxis_title="DateTime",
        yaxis_title="Number of Detections",
        # Add transparency to better see the shading
		plot_bgcolor="rgba(229, 236, 246, 1)" # background darker default color
    )
    
    p = PlotlyJS.plot(layout)
    
    # Add day/night shading if requested
    if add_daynight
        # Get unique dates in the data
        dates = unique(Date.(df[!, dts_colname]))
        
        # Add night rectangles (can use SunCalc.jl or Astro.jl for accurate sunrise/sunset)
        shapes = PlotlyJS.Shape[]
        
        for date in dates
			@debug "date: $date"
            # Approximate sunrise (6am) and sunset (6pm)
            # For more accuracy, use sunrise/sunset calculation based on lat/long
            sun_time = get_sun_times(date, latitude, longitude)
			@debug sun_time
			sunrise = sun_time.sunrise
            sunset = sun_time.sunset
			# sunrise = DateTime(date) + Hour(6)
            # sunset = DateTime(date) + Hour(18)
            @debug "sunrise: $sunrise, sunset: $sunset"
            # Day rectangle from surise to sunset
            push!(shapes, PlotlyJS.rect(
                x0 = sunrise,
                x1 = sunset,
                y0 = 0,
                y1 = 1,
                yref = "paper",
                fillcolor = "rgba(255, 255, 255, 0.4)",
                line_width = 0,
                layer = "below"
            ))
            
            # # Night rectangle from sunset to midnight
            # push!(shapes, PlotlyJS.rect(
            #     x0 = sunset,
            #     x1 = DateTime(date + Day(1)),
            #     y0 = 0,
            #     y1 = 1,
            #     yref = "paper",
            #     fillcolor = "rgba(55, 55, 150, 0.2)",
            #     line_width = 0,
            #     layer = "below"
            # ))
        end
        
        # Add shapes to layout
        # p.layout["shapes"] = shapes
		layout[:shapes] = shapes
		p = PlotlyJS.plot(layout)
    end
    @debug "HERE!!!!!!!!!!!!!!!!!!!!!"
    # Add data traces
	for (col, name) in trace_list
		if hasproperty(df, col)
			PlotlyJS.add_trace!(p, plottype(df, x=dts_colname, y=col; name=name))
		else
			@warn "Column $col not found in DataFrame"
		end
	end
    # PlotlyJS.add_trace!(p, plottype(df, x=:datetime, y=:num_tonal; name="tonal"))
    # PlotlyJS.add_trace!(p, plottype(df, x=:datetime, y=:num_noise; name="noise"))
    # PlotlyJS.add_trace!(p, plottype(df, x=:datetime, y=:num_impulsetrain; name="click_train"))
    # PlotlyJS.add_trace!(p, plottype(df, x=:datetime, y=:num_impulseINtrain; name="click"))

	isnothing(res_dir) && (return p)
    # Save output
    mkpath(res_dir)
    PlotlyJS.savefig(p, joinpath(res_dir,basename(dirname(fname))*".png"))
	PlotlyJS.savefig(p, joinpath(res_dir,basename(dirname(fname))*".html"))
    # open(joinpath(res_dir,basename(dirname(fname))*".html"), "w") do io
    #     PlotlyBase.to_html(io, p)
    # end
    return p
end

#~ get sun timings
using SunCalc, TimeZones, TimeZoneFinder
function get_sun_times(date, latitude, longitude; localtime=true)
    # times = SunCalc.getTimes(date, latitude, longitude)
    times = SunCalc.getTimes(date, latitude, longitude)
	@debug times
    if localtime
        return NamedTuple{keys(times)}(
            TimeZones.ZonedDateTime.(values(times),
                timezone_at(latitude, longitude)|>Ref; from_utc=true)
                .|> DateTime
            )

        # TimeZones.ZonedDateTime.([times[:sunrise], times[:sunset]],
        #         timezone_at(latitude, longitude)|>Ref; from_utc=true) .|> DateTime
    end
    
    return times #times[:sunrise], times[:sunset]
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


function plot_fft(snip, fs=1.0; type=:amplitude, plot=plot) 
	fft_val = rfft(snip, 1) .|> abs
	freqss =  fftfreq2(size(snip,1),fs)  #0:(fs/size(snip,1)):fs÷2

	type == :log && (fft_val = 20 .* log10.(fft_val))
	@debug (size(snip), size(freqss), size(fft_val), typeof(freqss), typeof(fft_val))
	plot(freqss, fft_val)
end

plot_fft!(args...; kwargs...) = plot_fft(args...; plot=plot!,kwargs...)

lay = Plots.@layout [a b];
function plot_time_fft(snip, fs=1.0; layout=lay, kwargs_plotfft=(), kwargs...)
	# a = plot(signal(snip,fs))
	a = plot(snip)
	b = plot_fft(snip,fs; kwargs_plotfft...)
	plot(a,b; layout=layout, kwargs...)
end

# norm_max(args...; norm_func=x->maximum(abs.(x); dims=1), kwargs...) = args[1]./norm_func(args[1])

plot_norm(args...; norm_func=x->maximum(abs.(x); dims=1), kwargs...) = plot(args[1]./norm_func(args[1]), args[2:end]...; kwargs...)
plot_norm!(args...; norm_func=x->maximum(abs.(x); dims=1), kwargs...) = plot!(args[1]./norm_func(args[1]), args[2:end]...; kwargs...)

function psd_plot(args...; kwargs...)
    pow, freq = psd2(args...; kwargs...)
    plot(freq, pow, xlabel="Frequency (Hz)", ylabel="Power Spectral Density (dB/Hz)")#, xscale=:log10)
end

function psd_plot!(args...; kwargs...)
    pow, freq = psd2(args...; kwargs...)
    plot!(freq, pow, xlabel="Frequency (Hz)", ylabel="Power Spectral Density (dB/Hz)")#, xscale=:log10)
end

function psd_plot_file(aufname; res_fol=missing, ch_list=nothing)
	data, fs = readAudio(aufname)
	if isnothing(ch_list) 
		ch_list = 1:size(data,2)
	end
	p = psd_plot(@view(data[:,ch_list]); fs=fs, nfft=1024*4)
	title!(basename(aufname))
	if !ismissing(res_fol) 
		savefig(joinpath(res_fol, (splitext(aufname)[1]|>basename) * ".html") )
	end
	return p
end


"""
Require: include("dsp.jl")
usage:
plot_ambient(aufname; ch_db=9, res_dir="/Users/abel/Documents/data_res/concretecho/Ambient/amb_balance", rx_vect=rx_vect)
"""
function plot_ambient(data_fs; res_dir=nothing, kwargs...)
	ch_list, correction, freqss, pow, pp = ambientnoise_correction(data_fs; kwargs...)

	p1=plot(freqss,pow); p2=plot(freqss,pp);
	p = plot(p1,p2; layout=Plots.@layout([a;b]))
	title!(basename(aufname))
	if !isnothing(res_dir) 
		savefig(p, joinpath(res_dir, (splitext(aufname)[1]|>basename) * ".html") )
		savefig(p, joinpath(res_dir, (splitext(aufname)[1]|>basename) * ".png") )
	end
	return p
end

plot_ambient(aufname::String; kwargs...) = plot_ambient(readAudio(aufname); kwargs...)


# 3D waterfall plot
"""
	waterfall_3d(clips, selections)

Generate a 3D waterfall plot visualizing audio clips.

# Arguments
- `clips::Matrix`: Matrix where rows represent samples and columns represent clips.
- `selections::Vector`: Vector of indices or identifiers for the clips to display.

# Returns
A 3D surface plot with sample indices on the x-axis, selection identifiers on the y-axis,
and amplitude values on the z-axis.

# Description
Creates a 3D visualization of multiple audio clips arranged in a waterfall-like structure,
allowing for visual comparison of signal patterns across different selections.

# Example
"""
function waterfall_3d(clips, xlabels=1:size(clips,1))#, x = repeat(1:size(clips, 1), 1, size(clips, 2)))
    n_samples = size(clips, 1)
	n_clips = size(clips,2)
	# Create meshgrid-like arrays
    x = repeat(xlabels, 1, n_clips)
    y = repeat(reshape(1:n_clips, 1, :), n_samples, 1)
    z = clips[:, 1:n_clips]
    
    Plots.surface(x, y, z, 
           title="3D Waterfall Plot",
           xlabel="Sample Index", 
           ylabel="repetition",
           zlabel="Amplitude")
end

# function waterfall_3d(clips, selections=1:size(clips, 2))
#     n_samples = size(clips, 1)
#     n_clips = length(selections)
    
#     # Create meshgrid-like arrays
#     x = repeat(1:n_samples, 1, n_clips)
#     y = repeat(reshape(selections, 1, :), n_samples, 1)
#     z = clips[:, 1:length(selections)]
    
#     Plots.surface(x, y, z, 
#            title="3D Waterfall Plot",
#            xlabel="Sample Index", 
#            ylabel="Selection",
#            zlabel="Amplitude")
# end

# Apply to your data:
# waterfall_3d(clips, selections)
# """
# respath = "/media/spin/anas2/data_res/dolphin/calf/Single_ball/res_20250910/20250227_101058/20250227_10.10.58_log_t251.28621036482684_d200__cps60.0.jld2"
# train_indx = [ x:( i+1<length(res.res_impulsetrain.train_start_ind) ? res.res_impulsetrain.train_start_ind[i+1]-1 : length(res.res_impulsetrain.pind_good)) for (i,x) in enumerate(res.res_impulsetrain.train_start_ind)]
# res_dir = joinpath(result_directory, "20250227_10.10.58")
# plot_color_clicks.(train_indx, 
#     "$res_dir/color_click/click_" .* [ "$(res.res_impulsetrain.pind_good_inS[train_indx[i][1]])_$(train_indx[i][1])" for i in 1:length(train_indx)] 
#         .* "..html";
#     rgbs_alpha_offset=0.2
#     )

# """
function plot_color_clicks(win_anal, savefname;
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
	!isempty(savefname) && savefig(savefname)
	return p
end


"""
    multispec_rgb(x, fs; nffts=[128,512,2048], hops = nothing,
                  db=true, dynrange=80.0, ref=1.0, gamma=1/2.2,
                  power=true, pad=true, norm=:per_channel)

Compute three spectrograms with different `nfft` values and pack them into RGB channels.

Returns:
  specs  :: NTuple{3,Matrix{Float64}}   (each power/amp spectrogram, size: fbin × frame)
  freqs  :: NTuple{3,Vector{Float64}}
  times  :: NTuple{3,Vector{Float64}}
  rgb    :: Array{Float32,3}  (F × T × 3, nearest-upscaled & 0–1)

Arguments:
  x        : 1-D signal
  fs       : sampling rate
  nffts    : vector of 3 FFT sizes (R,G,B)
  hops     : hop sizes (defaults to nfft ÷ 4 each)
  db       : convert to dB if true
  dynrange : clamp lower (max - dynrange) in dB scaling
  ref      : reference for dB (10*log10(P/ref))
  gamma    : gamma correction applied after 0–1 scaling
  power    : if true use |X|^2, else |X|
  pad      : reflect pad to fit last frame
  norm     : :per_channel (independent scaling) or :global

Example:
    specs, freqs, times, rgb = multispec_rgb(sig, fs; nffts=[128,512,2048])
"""
function multispec_rgb(x::AbstractVector, fs;
    nffts = [128,512,2048],
    hops::Union{Nothing,AbstractVector}=nothing,
    db=true, dynrange=80.0, ref=1.0, gamma=1/2.2,
    power=true, pad=true, norm=:per_channel)

    @assert length(nffts)==3 "Need exactly 3 nfft values"
    hops === nothing && (hops = nffts .÷ 4)

    # Hann window helper
    hann(n) = 0.5 .- 0.5*cos.(2π*(0:n-1)/(n-1))

    function one_spec(x, fs, nfft, hop)
        w = hann(nfft)
        L = length(x)
        if pad && L < nfft
            xpad = vcat(x, zeros(eltype(x), nfft-L))
            L = length(xpad); xuse = xpad
        else
            xuse = x
        end
        nframes = 1 + max(0, (L - nfft) ÷ hop)
        nfreq = nfft ÷ 2 + 1
        S = Matrix{Float64}(undef, nfreq, nframes)
        for k in 0:nframes-1
            i1 = k*hop + 1
            i2 = i1 + nfft - 1
            if i2 > L
                if pad
                    seg = similar(xuse, nfft)
                    nremain = L - i1 + 1
                    seg[1:nremain] .= @view xuse[i1:end]
                    seg[nremain+1:end] .= 0
                else
                    break
                end
            else
                seg = @view xuse[i1:i2]
            end
            frame = seg .* w
            spec = rfft(frame)
            mag = abs.(spec)
            power && (mag .*= mag)
            S[:, k+1] = mag
        end
        # Compute times at frame centers
        times = ( (0:size(S,2)-1).*hop .+ (nfft/2) ) ./ fs
        freqs = (0:nfreq-1) .* (fs/nfft)
        return S, freqs, times
    end

	raw_specs = Vector{Matrix{Float64}}(undef, 3)
	freqlist  = Vector{Vector{Float64}}(undef, 3)
	timelist  = Vector{Vector{Float64}}(undef, 3)

	for (i,(nfft,hop)) in enumerate(zip(nffts,hops))
		S,freqs,times = one_spec(x, fs, nfft, hop)
		raw_specs[i] = S
		freqlist[i]  = freqs
		timelist[i]  = times
	end

	raw_specs_tuple = (raw_specs[1], raw_specs[2], raw_specs[3])
	freqlist_tuple  = (freqlist[1], freqlist[2], freqlist[3])
	timelist_tuple  = (timelist[1], timelist[2], timelist[3])

    # dB & scaling
	proc_specs = map(raw_specs_tuple) do S
        if db
            # SdB = 10 .* log10.(S .+ eps()) .- 10*log10(ref)
			SdB = log10.(S)
            mx = maximum(SdB)
            SdB_clamped = clamp.(SdB, mx - dynrange, mx)
            A = (SdB_clamped .- (mx - dynrange)) ./ dynrange
        else
            # amplitude/power direct normalization
            mx = maximum(S)
            A = mx>0 ? S./mx : S
        end
        gamma == 1 ? A : A .^ gamma
    end

    # Combine with nearest-neighbor resizing to largest dimensions
    target_f = maximum(size(S,1) for S in proc_specs)
    target_t = maximum(size(S,2) for S in proc_specs)
    rgb = Array{Float32,3}(undef, target_f, target_t, 3)

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

    for c in 1:3
        rgb[:,:,c] = nn_resize(proc_specs[c], target_f, target_t)
    end

    if norm == :global
        mx = maximum(rgb)
        mx>0 && (rgb ./= mx)
    end

	return raw_specs_tuple, freqlist_tuple, timelist_tuple, rgb
end

"""
    multispec_rgb_plot(x, fs; kwargs...)

Convenience wrapper: runs `multispec_rgb` and returns an RGB image matrix
(Height × Width × 3, Float32 in 0–1).
"""
function multispec_rgb_plot(x, fs; kwargs...)
    _, _, _, rgb = multispec_rgb(x, fs; kwargs...)
    return rgb
end

# # Example (uncomment to test):
# using FileIO, ImageCore

# sig = randn(10_000)
# sig = chirp(1000, 2000, 3, 44100) |> real
# # specs, freqs, times, rgb = multispec_rgb(sig, 48_000; nffts=[128,256,2048])
# # imm = colorview(RGB, permutedims(rgb, (3,1,2)))
# # save("temp/multispec.png", imm)

# # specs .|> size
# # freqs .|> size
# # times .|> size
# # rgb |> size

# # rgb[:,:,1] |> Plots.heatmap
# # rgb[:,:,2] |> Plots.heatmap
# # rgb[:,:,3] |> Plots.heatmap
# # Plots.heatmap(specs[1])
# # Plots.heatmap(specs[2])
# # Plots.heatmap(specs[3])
# # # specgram(clips[i]; fs=fs, colorbar=nothing)


# specgram(sig)

# fs=44100
# nffts = [128,512,2048]; hops = zeros(Int,length(nffts))
# S=[]; F=[]; T=[];
# for (nfft, hop) in zip(nffts, hops)
# 	# @info typeof.([sig, nfft, hop])
#     s = stft(sig, nfft, hop; window=hann) .|> abs
#     push!(S, s)
#     # push!(F, f)
#     # push!(T, t)
# end

# target_f = maximum(size(s,1) for s in S)
# target_t = maximum(size(s,2) for s in S)
# rgb = Array{Float32,3}(undef, target_f, target_t, 3)
# rgb .= 0
# function nn_resize(S, F, T)
# 	fsrc, tsrc = size(S)
# 	out = Matrix{Float32}(undef, F, T)
# 	for j in 1:T
# 		tj = clamp(round(Int, (j-1)/(T-1) * (tsrc-1) + 1), 1, tsrc)
# 		for i in 1:F
# 			fi = clamp(round(Int, (i-1)/(F-1) * (fsrc-1) + 1), 1, fsrc)
# 			# @info i, j, fi, tj
# 			out[i,j] = S[fi, tj]
# 		end
# 	end
# 	out
# end

# for c in 1:3
# 	rgb[:,:,c] = nn_resize(S[c], target_f, target_t)
# 	rgb[:,:,c] ./= maximum(rgb[:,:,c])
# end

# # rgb to RGB image
# img = colorview(RGB, permutedims(rgb, (3,1,2))[:, end:-1:1, :])
# save("temp/multispec_img_hann-green.png", img)

# rgb[:,:,1] |> Plots.heatmap
# savefig("temp/multispec_1.html")
# rgb[:,:,2] |> Plots.heatmap
# savefig("temp/multispec_2.html")
# rgb[:,:,3] |> Plots.heatmap
# savefig("temp/multispec_3.html")

# # save img

# specgram(sig; fs=fs, nfft=128, noverlap=0, colorbar=nothing
# 	, downsample=nothing, pooling=nothing)
# savefig("temp/multispec_specgram_128.html")