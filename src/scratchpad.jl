include("audio.jl")
include("aspod.jl")
using VideoIO
include("readImages.jl")
include("utils.jl")
include("synchronization.jl")

aufname = "/Users/abel/Documents/data/calf/coop/20231128/20231128_15.16.48_log.flac"
vidfname = "/Users/abel/Documents/data/calf/coop/20231128/20231128_15.16.48_log.mkv"
res_dir = "/Users/abel/Documents/data_res/calf/coop/2023-11-28_3"
res_dir = "/Users/abel/Documents/data_res/calf/coop/2023-11-28_6-v0.0.5"
res_dir = "/Users/abel/Documents/data_res/calf/coop/2023-11-28_7-v0.0.5b-whistle"
res_dir = "/Users/abel/Documents/data_res/calf/coop/2023-11-28_7-v0.0.5c-whistle-alpha"

d = load(filter(x -> (y=splitext(x); y[2]==".jld2" && occursin(splitext(basename(aufname))[1], basename(y[1])) ), readdir(res_dir; join=true)))

data, fs, _,_, timestamp = readAudio(aufname)

dets = process_detections(aufname, vidfname; res_dir=res_dir)

d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t116.9741055601481_d15000__cps0.375.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t116.9741055601481_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t12_d15000.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t12_d15000.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t116.15144729819004_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t116.15144729819004_d15000__cps0.375.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t12_d15000.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t116.14395330535885_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t116.14395330535885_d15000__cps0.375.jld2"])
plot_ang(d["res_impulsetrain"], d["ang_impulsive"][1][1]; label=["azimuth" "inclination"], type="Power")

d = load(["/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev_old/20231128_15.16.48_log_t4_d800.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev_old/20231128_15.16.48_log_t46.36293944139811_d800__cps60.0_angles.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev_old/20231128_15.16.48_log_t46.36293944139811_d800__cps60.0.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/coop/2023-11-28/20231128_15.16.48_log_t28.855573029963317_d800__cps60.0.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/20231128_15.16.48_log_t28.855573029963317_d800__cps60.0_angles.jld2"])

data_filt = filter_simple(data[:,get_relevant_channels(rx_vect)], impulsive_band_pass; fs=fs)
data_filt = filter_simple(data[:,get_relevant_channels(rx_vect)], [5_000 180_000], [[30_000,31_000],[61_000,62_000]]; fs=fs)
data_filt = filter_simple(data[:,get_relevant_channels(rx_vect)], [70_000 130_000]; fs=fs, butterworth_size=32)

data_filt = mapslices2(x->denoise(x.|>Float64, TI=true), data_filt)

i = 310#400 #easier
i = 1160
i = 11029
i = 10952
i = 7556
i = 10536
i = 5410 # easier
i = 1685
i = 3409
plot_one_event(res_dir, vidfname, d["res_impulsetrain"].pind_good_inS, data_filt, d["res_impulsetrain"].pind_good, tdoas, window_impulsive, bearings; i=i, plotfunc=gr)


plot_one_event(res_dir, vidfname, d["res_impulsetrain"].pind_good_inS, data_filt, d["res_impulsetrain"].pind_good, d["ang_impulsive"][2], window_impulsive, d["ang_impulsive"][1][1]; i=i, plotfunc=gr)
img = readImage(vidfname,d["res_impulsetrain"].pind_good_inS[i])
snip = windowing(data_filt, d["res_impulsetrain"].pind_good[i], window_impulsive);
norm_snip = norm_max(snip)
# snip = mapslices2(x->denoise(x.|>Float64, TI=true), snip)

plotlyjs(); plot_norm(snip);  vline2(tdoas[i,:], ans, 1); plot!()
tdoa = findTrigger.(eachcol(snip), 1)'; vline2(tdoa, ans, 1)#maximum(snip)) #tdoa |> vline!

tdoa = findTrigger.(eachcol(abs.(hilbert(snip))), 1)'; vline2(tdoa, ans, 1)#maximum(snip)) #tdoa |> vline!

tdoas = get_tdoa_findTrigger(data_filt, d["res_impulsetrain"].pind_good; window=window_impulsive)
tdoas = get_tdoa_raw_MaxEnergyRefChannel_resample(data_filt, d["res_impulsetrain"].pind_good; window=window_impulsive)

@time tdoas = get_tdoa_findTrigger_waveletdenoise(data_filt, d["res_impulsetrain"].pind_good; window=window_impulsive)
bearings = tdoa2dir(tdoas, rx_vect, fs)
plot_one_event(res_dir, vidfname, d["res_impulsetrain"].pind_good_inS, data_filt, d["res_impulsetrain"].pind_good, tdoas, window_impulsive, bearings; i=i, plotfunc=gr)

[(plot(snip[:,ch]); plot!(abs.(hilbert(snip[:,ch]))); hline!([median(abs.(hilbert(snip[:,ch]))) * .75]); title!(string(ch))|>display) for ch in 1:3]
tdoas[i,:] = tdoa #[259 246 252] #[245 238 239]
plot(snip |> norm_max); vline2(tdoas[i,:], ans, 1); plot!()#; hline!([median(abs.(hilbert(snip|>norm_max)), dims=1) .*.75])
plot( mapslices(x->denoise(x, TI=true), snip, dims=1) |> norm_max); vline2(tdoas[i,:], ans, 1); hline!([median(abs.(hilbert(snip|>norm_max)), dims=1) .*.75])

@time tdoas, ref_signals = get_tdoa_raw_MaxPeakRefChannel(data_filt, d["res_impulsetrain"].pind_good; window=window_impulsive)

# vline!([tdoa[i] for i in 1:length(tdoa)], color = [palette(:default)[i] for i in 1:length(tdoa)], linewidth = 3)
# vline!(tdoa;  color=palette(:default)[1:3] , width=3)
finddelay
snip = data_filt[window_impulsive .+ d["res_impulsetrain"].pind_good[i],get_relevant_channels(rx_vect)]
# plot_time_fft(snip, fs)
plotlyjs();
plot(snip)
vline!(d["ang_impulsive"][2][i,:])
plot(abs.(hilbert(snip)))
chs = get_relevant_channels(rx_vect)

# tdoas = []
# tdoas = Array{NamedTuple}(undef, length(chs))
# tdoas = Dict{Int, Any}()
tdoas_full = Vector{Any}(undef, length(chs));
Threads.@threads for ch in chs
    tdoas_full[ch] = findsignal.( Ref(signal(snip[:,ch],fs)), eachcol(signal(snip,fs)); prominence=0.1, finetune=2, mfo=true)# finetune=20)
end
tdoas = map( x-> map( y-> isempty(y.time) ? NaN : y.time[1], x), tdoas_full)
map( x-> (x.-x[1]) .* fs, tdoas)
# tdoas
nonan = map(x->isnan.(x) |> sum, tdoas)
goodindex = findall(==(0), nonan)
tdoas = tdoas[goodindex]



tdoas = findsignal.( Ref(signal(snip[:,ch],fs)), eachcol(signal(snip,fs)); prominence=0.25)# ; mfo=true)# finetune=20)
map(x->x.mfo .|> abs, tdoas_full) |> plot; title!("Channel "*string(ch))
plots = plot.( map( tdoa-> map(x->x.mfo .|> abs, tdoa), tdoas_full))

for (i, p) in enumerate(plots)
    plot(p; title="Channel "*string(i)) |> display
end

tdoas = findsignal2.( Ref(snip[:,ch]), eachcol(snip))# ; mfo=true)# finetune=20)

get_tdoa_raw.( Ref(data_filt), Ref(d["res_impulsetrain"].pind_good[i]), chs ; window=window_impulsive)
get_tdoa_envelope.( Ref(data_filt), Ref(d["res_impulsetrain"].pind_good[i]), chs ; window=window_impulsive)
get_tdoa_max( data_filt, d["res_impulsetrain"].pind_good[i], 1 ; window=window_impulsive)
get_tdoa_min( data_filt, d["res_impulsetrain"].pind_good[i], 1 ; window=window_impulsive)
get_tdoa_raw_MaxEnergyRefChannel(data_filt, d["res_impulsetrain"].pind_good[i]; window=window_impulsive)
tdoa = get_tdoa_minmax(data_filt, d["res_impulsetrain"].pind_good[i]; window=window_impulsive)
tdoa .- tdoa[1]
extrema_and_indices(snip)
# get_tdoa_raw(data_filt, d["res_impulsetrain"].pind_good[i] ; window=window_impulsive, ref_channel=ch)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vid = VideoIO.openvideo(vidfname)
vid_resolution = raw_frame_size(vid)

xx = LinRange(-fov_angle[1]÷2,fov_angle[1]÷2, fov_angle[1]+1); yy = LinRange(-fov_angle[2]÷2,fov_angle[2]÷2, fov_angle[2]+1);

xx = LinRange(-fov_angle[1]÷2,fov_angle[1]÷2, vid_resolution[1]); yy = LinRange(-fov_angle[2]÷2,fov_angle[2]÷2, vid_resolution[2]);
θ = deg2rad.(reduce(vcat, hcat.(xx', yy)))
sd = steering2(rx_vect, c, θ)

bfo = beamform_output( signal(snip, fs), sd; output_func=x->-reduce(-, extrema(x)) )
reshape(bfo, length(xx), length(yy))' |> x -> heatmap(x; aspect_ratio=:equal)

bfo = beamform( signal(snip, fs), sd)
# pow_output = p2p_db(bfo|>collect)
# reshape(pow_output, length(xx), length(yy))' |> x -> heatmap(x; aspect_ratio=:equal)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include("dsp.jl")
xt = filter_simple(data,[1000, 160_000]; fs=fs)
xt = mapslices2(x->denoise(x.|>Float64, TI=true), xt)

sample_dir = "/Users/abel/Documents/data/calf/sample"
aufnames = filter(endswith(r".flac|.wav"), readdir(sample_dir; join=true) )
aufname = aufnames[9]
data, fs = readAudio(aufname)

iseven(size(data,2)) || (data = data[begin:end-1,:])
xt1 = filter_simple(data,[1000, 160_000]; fs=fs)
xt = mapslices2(x->denoise(x.|>Float64, TI=true), xt1)

plot([xt1 xt])
plot(signal([xt1 xt], fs))

plot(norm_snip[:,1])
plot!([zeros(70,); norm_snip[:,2]])
plot!([zeros(14,); norm_snip[:,3]])

resample_ratio = 10
snip_r = resample(snip, resample_ratio; dims=1)
snip_r = snip_r |> norm_max
plot(snip_r[:,1])
plot!([zeros(689,); snip_r[:,2]])
plot!([zeros(140,); snip_r[:,3]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

extrema_in_file.(["/Users/abel/Documents/data/calf/coop/20230927", "/Users/abel/Documents/data/calf/coop/20231120", "/Users/abel/Documents/data/calf/coop/20231128"]; res_dir="/Users/abel/Documents/data_res/calf/test3")

extrema_in_file("/Users/abel/Documents/data/calf/coop/20231120", "/Users/abel/Documents/data_res/calf/test")


process_dir(["/Users/abel/Documents/data/calf/coop/20230927", "/Users/abel/Documents/data/calf/coop/20231120", "/Users/abel/Documents/data/calf/coop/20231128"]; func=extrema_in_file, arg="/Users/abel/Documents/data_res/calf/test")

process_dir("/Users/abel/Documents/data/calf/Clicktest"; func=process_folder, arg="/Users/abel/Documents/data_res/calf/click_test/new_res_20231219_2")



#~ amplification
foldername = "/Users/abel/Documents/data_res/calf/amplification/Eszter-2_res/Gn-Ao"
using JSON, DataFrames, Glob

# Get a list of all JSON files in the directory
files = glob("*.json", foldername)

# Initialize an empty DataFrame
df = DataFrame(aufname = String[], energy1 = Float64[], energy2 = Float64[], energy3 = Float64[], energy4 = Float64[], 
               extrema1 = Float64[], extrema2 = Float64[], extrema3 = Int64[], extrema4 = Int64[])

# Loop over the files
for file in files
    # Read the JSON file
    data = JSON.parsefile(file)

    # Append the data to the DataFrame
    push!(df, (data["aufname"], data["energy"][1], data["energy"][2], data["energy"][3], data["energy"][4], 
               data["extrema_and_indices"][1][1], data["extrema_and_indices"][1][2], data["extrema_and_indices"][1][3], data["extrema_and_indices"][1][4]))
end

##########
foldername = "/Users/abel/Documents/data_res/calf/amplification/Eszter-2_res/Gn-Ao"
foldername = "/Users/abel/Documents/data_res/calf/amplification/Eszter-2_res/To-Do-An"
using JSON, DataFrames, Glob

# Get a list of all JSON files in the directory
files = glob("*.json", foldername)

# Initialize an empty DataFrame
df = DataFrame(aufname = String[], energy1 = Float64[], energy2 = Float64[], energy3 = Float64[], energy4 = Float64[], 
               extrema1_1 = Float64[], extrema1_2 = Float64[], extrema1_3 = Int64[], extrema1_4 = Int64[],
               extrema2_1 = Float64[], extrema2_2 = Float64[], extrema2_3 = Int64[], extrema2_4 = Int64[],
               extrema3_1 = Float64[], extrema3_2 = Float64[], extrema3_3 = Int64[], extrema3_4 = Int64[],
               extrema4_1 = Float64[], extrema4_2 = Float64[], extrema4_3 = Int64[], extrema4_4 = Int64[])

# Loop over the files
for file in files
    # Read the JSON file
    data = JSON.parsefile(file)
    if occursin("Tap_test", data["aufname"])
        @info "skip: " * data["aufname"]
        continue
    end
    # Append the data to the DataFrame
    push!(df, (data["aufname"], data["energy"][1], data["energy"][2], data["energy"][3], data["energy"][4], 
               data["extrema_and_indices"][1][1], data["extrema_and_indices"][1][2], data["extrema_and_indices"][1][3], data["extrema_and_indices"][1][4],
               data["extrema_and_indices"][2][1], data["extrema_and_indices"][2][2], data["extrema_and_indices"][2][3], data["extrema_and_indices"][2][4],
               data["extrema_and_indices"][3][1], data["extrema_and_indices"][3][2], data["extrema_and_indices"][3][3], data["extrema_and_indices"][3][4],
               data["extrema_and_indices"][4][1], data["extrema_and_indices"][4][2], data["extrema_and_indices"][4][3], data["extrema_and_indices"][4][4]))
end

plotlyjs()
Plots.plot(df[:,[6,7]] |> Matrix; color="red");
Plots.plot!(df[:,[10,11]] |> Matrix; color="green");
Plots.plot!(df[:,[14,15]] |> Matrix; color="blue");
a=Plots.plot!(df[:,[18,19]] |> Matrix; color="yellow")

b=Plots.plot(df[:, 2:4] |> Matrix)

Plots.plot(a,b, layout = Plots.@layout[a;b] )


sha = @async run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Shakeela/2024/03";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Shakeela/03"); func=stack_audio_videos)
ella = @async run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Ella/2024/03";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Ella/03"); func=stack_audio_video, skipdone=true)
run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Shiye/2024/01";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Shiye/01"); func=stack_audio_videos, skipdone=true); run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Shiye/2024/02";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Shiye/02"); func=stack_audio_videos, skipdone=true)

run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/WeiShi/2024/04";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/WeiShi/04"); func=stack_audio_videos)

Base.throwto(sha, InterruptException())

ella
sha



using Images

# Define your function here
function my_function(x, y)
    r = x / 1000  # Red channel
    g = y / 1000  # Green channel
    b = (x + y) / 2000  # Blue channel
    return RGB(r, g, b)
end

# Create the image
image = [my_function(x, y) for y in 1:1000, x in 1:1000]

# Save the image
save("image.png", image)