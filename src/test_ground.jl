# ang = tdoa2dir([d["tdoa"];], rx_vect, fs)

using FileIO
using VideoIO
using DSP
include("audio.jl")
include("localization.jl")
include("plotting.jl")
include("readImages.jl")

d = load("/Users/abel/Documents/data_res/aspod/bahamas_2023/aspod/2022-05-25/0003/Aud_20131218_213638_t819.1750000000001_d200.jld2")
aufname = "/Users/abel/Documents/data/aspod/field/bahamas_2023/aspod/2022-05-25/0003/Aud_20131218_213638.wav"
vidfname = "/Users/abel/Documents/data/aspod/field/bahamas_2023/aspod/2022-05-25/0003/Vid_20131218_213637.mkv"
res_dir = "/Users/abel/Documents/data_res/aspod/bahamas_2023/aspod/2022-05-25/0003"

data, fs = readAudio(aufname)
filter_weight = digitalfilter(Bandpass(band_pass[1],band_pass[2], fs=fs), Butterworth(butterworth_size))
data_filt = mapslices( x -> filtfilt( filter_weight, x), data[:,1:end-1], dims=1)

ang = tdoa2dir([d["tdoa"];], rx_vect, fs)

plot_one_event(joinpath(res_dir,splitext(basename(vidfname))[1]), 
vidfname, d["pind_good_inS"], data_filt, d["pind_good"], 
d["tdoa"], -100:500, ang; func2=plotTDOA_raw, plotfunc=gr, i=ind)


# plot_all_clicks(joinpath(res_dir,splitext(basename(vidfname))[1]), 
# vidfname, d["pind_good_inS"], data_filt, d["pind_good"], 
# d["tdoa"], -100:500, ang; func2=plotTDOA_raw, plotfunc=gr)


fig = Figure()
# ax1 = Axis(fig[1,1])
pts = []
get_pixelLoc(pts, fig, 
vidfname, d["pind_good_inS"], data_filt, d["pind_good"], 
d["tdoa"], -100:500, ang)


points = Observable(Point2f[])

scene = Scene(camera = campixel!)
linesegments!(scene, points, color = :black)
scatter!(scene, points, color = :gray)

pts = []
on(events(scene).mousebutton) do event
    if event.button == Mouse.left
        if event.action == Mouse.press# || event.action == Mouse.release
            mp = events(scene).mouseposition[]
            push!(pts, mp)
            @info mp
            # notify(points)
        end
    end
end

scene


fig = Figure()
ax1 = Axis(fig[1,1])
image!(ax1,img)
pts = []

register_interaction!(ax1, :my_interaction) do event::MouseEvent, axis
    if event.type === MouseEventTypes.leftclick
        println("You clicked on the axis at datapos $(event.data)")
        push!(pts, event.data)
    end
end


