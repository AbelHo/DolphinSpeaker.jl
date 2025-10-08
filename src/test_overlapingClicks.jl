include("detector.jl")
include("plotting.jl"); 
plotlyjs()
include("tabulate_data.jl")
include("test_make_video.jl")
include("audio.jl")
using ImageShow

aufname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/20250227_10.16.59_log.flac"
respath = "/media/spin/anas2/data_res/dolphin/calf/new_2025_freeplay_report/20250227_10.16.59_log_t292.4385100072437_d200__cps60.0.jld2" #"/media/spin/anas2/data_res/dolphin/calf/new_2025_freeplay_report/20250227_10.16.59_log_t292.43851000724715_d200__cps60.0.jld2"
vidfname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/20250227_10.16.59_log.mp4"
res_dir = "temp"

# new set
aufname = "/media/spin/anas2/data/calf/upload/Single_ball_S1/20250221_10.53.41_log.flac"
respath = "/media/spin/anas2/data_res/dolphin/calf/res_20250904/ball2/20250221_10.53.41_log_t313.82475408659667_d200__cps60.0.jld2"
vidfname = "/media/spin/anas2/data/calf/upload/Single_ball_S1/20250221_10.53.41_log.mkv"
res_dir = dirname(respath)

isfile(aufname)
isfile(respath)
isfile(vidfname)
res_dir = joinpath(res_dir, splitext(basename(aufname))[1])

autocor_clips, plot_dir = analyze_and_plot_clips(
    aufname,
    "",
    respath,
    impulsive_band_pass,
    [-25 102];
    output_types = ["png", "html"],
    plot_dir_prefix = joinpath(res_dir,"clips/OVERLAP_clips_train2_"),
    flag_extra_plot = true
)

ftype = "html"
clips_plot_dir = "temp/clips/OVERLAP_clips_train2_20250227_101659"
make_clip_index_html(clips_plot_dir; outname=basename(clips_plot_dir)*"_$(ftype)_all.html", output_type=ftype, prefix="all_")

data, fs, _, _, timestamp = readAudio(aufname)
data_filt = filter_simple(data, impulsive_band_pass; fs=fs)

res = load(respath)
res = dict2namedtuple(res)

# Extract clips
nfunc(x) = sqrt(sum(abs2.(x)))
window_extract = [-25 102]
clips_fixed = extract_clips(data_filt, window_extract .+ res.res_impulsetrain.pind_good, NaN; flag_matrix=false)
data=nothing;data_filt=nothing; GC.gc()

# clips = hcat(map(x-> x[:,ref_channel], clips_fixed)...)
# select loudest channel
clips = hcat(map(x-> x[:, argmax(maximum(abs.(x); dims=1))[2] ], clips_fixed)...)
win_anal = 521:590 #52681:52710 #624:635 #594:655 #9771:9800
plot_time_fft(clips[:,win_anal],fs; size=(1000,600))
plot_time_fft(clips[:,win_anal]|>norm_max,fs; size=(1000,600))
xcs = finddelay2.(clips[:,win_anal[1]] |> Ref, clips[1:end,win_anal]|>eachcol; norm_func=x->norm_max(x; norm_func=nfunc) )
plot( map(x->x[2], xcs); size=(1000,600), label="max corr");
plot!(clips[:,win_anal] |> eachcol .|> energy |> norm_max, label="energy")
res.res_impulsetrain.pind_good[win_anal] |> diff |> plot
time_diff = detect_impulsegroup(res.res_impulsetrain.pind_good[win_anal])
Plots.histogram(filter(>(0), time_diff), bins=0:20:10000)

td = filter(>(0), time_diff)# |> plot #lots.histogram
td2 = time_diff |> copy
# change all values <1 and >600 into zeros
td2[td2 .< 1] .= 0
td2[td2 .> 1200] .= 0

Plots.histogram(filter(>(0), time_diff), bins=0:20:10000)

# add color based on frequencies
rgbs, rgba = sig2rgb(clips[:,:]; fs=fs, 
    rgb_bands=[[10_000, 60_000], [60_000, 110_000], [110_000, 160_000]])
rgbs[4,:] .+= 0.2
plot( res.res_impulsetrain.pind_good_inS[win_anal], ones(win_anal|>length,); 
    color=a, seriestype=:scatter, bg=:black,
    hover = string.(win_anal))

plot( res.res_impulsetrain.pind_good_inS[win_anal], ones(win_anal|>length,); 
    color=rgba[win_anal], seriestype=:scatter, bg=:black,
    hover = string.(win_anal))

# custom plot segment
win_anal = 656:1001 #1:length(res.res_impulsetrain.pind_good_inS)    #4974:6090 #5000:5400
rgbs, rgba = sig2rgb(clips[:,win_anal]; fs=fs, 
    rgb_bands=[[10_000, 60_000], [60_000, 110_000], [110_000, 160_000]])
rgbs[4,:] .+= 0.2
# plot(res.res_impulsetrain.pind_good_inS[win_anal], rgbs[4,:]; 
#     color=rgba, seriestype=:scatter, bg=:black,
#     hover = string.(win_anal),
#     size=(1000,600))
plot(res.res_impulsetrain.pind_good_inS[win_anal], rgbs[4,:]; 
    color=rgbs|> eachcol .|> x-> RGBA(x...), 
    # color=rgba,
    seriestype=:scatter, 
    hover = string.(win_anal),
    bg=:black,markerstrokewidth = 0,
    size=(1000,600))
savefig("temp/color_click/temp_$(res.res_impulsetrain.pind_good_inS[win_anal[1]]).html")

# run each click train
train_indx = [ x:( i+1<length(res.res_impulsetrain.train_start_ind) ? res.res_impulsetrain.train_start_ind[i+1]-1 : length(res.res_impulsetrain.pind_good)) for (i,x) in enumerate(res.res_impulsetrain.train_start_ind)]

"temp/color_click/temp_$(res.res_impulsetrain.pind_good_inS[train_indx]).html"



plot_color_clicks.(train_indx, 
    "$res_dir/color_click/click_" .* [ "$(res.res_impulsetrain.pind_good_inS[train_indx[i][1]])_$(train_indx[i][1])" for i in 1:length(train_indx)] 
        .* "..html";
    rgbs_alpha_offset=0.2
    )


d = load("/media/spin/anas2/data_res/dolphin/calf/new_2025_freeplay_report/20250227_10.16.59_log_t292.4385100072437_d200__cps60.0_angles.jld2") #20250227_10.16.59_log_t292.43851000724715_d200__cps60.0_angles.jld2")
d = dict2namedtuple(d);

# vidfname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/20250227_10.16.59_log.mp4"
vid = VideoIO.openvideo(vidfname)
img = get_image(vid, 0)

indx = findall(==(5), d.pixel_related_impulsive[1])
# display_corners!(img, d.pixel_related_impulsive[2][indx,:] |> eachrow; dot_size=25, colors=[1,0,0], opacity=0.4)
# img

Plots.scatter(d.pixel_related_impulsive[2][indx,1], d.pixel_related_impulsive[2][indx,2]; 
    color=rgba[indx], markerstrokewidth=0, legend=false); hline!([size(img,1)]); vline!([size(img,2)])

mkpath("temp/vid_rope2")
fps = get_fps(vidfname)
# for ind = minimum(d.pixel_related_impulsive[1]):maximum(d.pixel_related_impulsive[1])
ind = 206; min_alpha = 1
    indx = findall(==(ind), d.pixel_related_impulsive[1]);
    Plots.scatter(d.pixel_related_impulsive[2][indx,1], d.pixel_related_impulsive[2][indx,2]; 
        color=rgbs[:,indx] .+ [zeros(3,length(indx)); ones(1,indx|>length)*min_alpha] |> eachcol .|> x-> RGBA(x...), 
        markerstrokewidth=0, legend=false, title="Vid Frame: $ind, $(round(ind/fps; digits=3))s");
        hline!([0, size(img,1)]); vline!([0, size(img,2)]) #|> display
    savefig("temp/vid_rope2/frame_$ind.png")
# end

display_corners3!(img, d.pixel_related_impulsive[2][indx,:], 50, [1,0,0], 1.0)
img

using Plots, ColorTypes

# Example data
x = rand(10)
y = rand(10)
colors = [RGB(rand(), rand(), rand()) for _ in 1:10]  # List of RGB colors for each point

Plots.scatter(x, y; color=rgba[indx] .|> RGB, markerstrokewidth=0, legend=false)


#~ boris behavior labels
using CSV, DataFrames, CategoricalArrays
include("video.jl")
behavior_label_fname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/20250227_10.16.59_log_overlaid_t9_d200.mkv_normalized-audio_boris.csv"
df = CSV.read(behavior_label_fname, DataFrame)
filter(x->x.Behavior == "Biting", df)
categories = df.Behavior |> CategoricalArray 
df.categoryID = categories .|> levelcode
[categories|>unique df.categoryID|>unique]

audacity_label([df[:,"Start (s)"] df[:,"Stop (s)"]], splitext(aufname)[1]*"_behavior-audacity.txt", 
    df.Behavior .* "(" .* df.Subject .*")")

audacity_label( filter(x->x.Behavior == "Biting", df), splitext(aufname)[1]*"_behavior-audacity_biting.txt")


# make clips
vidfname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/20250227_10.16.59_log.mp4"
res_dir = "/media/spin/anas2/data_res/dolphin/calf/new_2025_freeplay_report/video/tarsier/6s"
make_video_clips_ffmpeg(df, vidfname; outdir=res_dir, min_duration=6.0)

res_dir = "/media/spin/anas2/data_res/dolphin/calf/new_2025_freeplay_report/video/tarsier/6s_all"
split_video_by_duration(vidfname, 6.0; outdir=res_dir)