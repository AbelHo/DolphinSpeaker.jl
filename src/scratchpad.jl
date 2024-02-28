include("audio.jl")
include("aspod.jl")
using VideoIO
include("readImages.jl")
include("utils.jl")

aufname = "/Users/abel/Documents/data/calf/coop/20231128/20231128_15.16.48_log.flac"
vidfname = "/Users/abel/Documents/data/calf/coop/20231128/20231128_15.16.48_log.mkv"
res_dir = "/Users/abel/Documents/data_res/calf/coop/2023-11-28/fixed_minmax"

data, fs, _,_, timestamp = readAudio(aufname)

dets = process_detections(aufname, vidfname; res_dir=res_dir)

d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t116.9741055601481_d15000__cps0.375.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t116.9741055601481_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t12_d15000.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t12_d15000.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t116.15144729819004_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t116.15144729819004_d15000__cps0.375.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t12_d15000.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t116.14395330535885_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t116.14395330535885_d15000__cps0.375.jld2"])
plot_ang(d["res_impulsetrain"], d["ang_impulsive"][1][1]; label=["azimuth" "inclination"], type="Power")

d = load(["/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev_old/20231128_15.16.48_log_t4_d800.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev_old/20231128_15.16.48_log_t46.36293944139811_d800__cps60.0_angles.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev_old/20231128_15.16.48_log_t46.36293944139811_d800__cps60.0.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/coop/2023-11-28/20231128_15.16.48_log_t28.855573029963317_d800__cps60.0.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/20231128_15.16.48_log_t28.855573029963317_d800__cps60.0_angles.jld2"])

data_filt = filter_simple(data[:,1:size(rx_vect,2)], impulsive_band_pass; fs=fs)
data_filt = filter_simple(data[:,1:size(rx_vect,2)], [5_000 180_000], [[30_000,31_000],[61_000,62_000]]; fs=fs)

i = 310#400
i = 1160
plot_one_event(res_dir, vidfname, d["res_impulsetrain"].pind_good_inS, data_filt, d["res_impulsetrain"].pind_good, d["ang_impulsive"][2], window_impulsive, d["ang_impulsive"][1][1]; i=i, plotfunc=gr)


snip = data_filt[window_impulsive .+ d["res_impulsetrain"].pind_good[i],1:size(rx_vect,2)]
# plot_time_fft(snip, fs)
plotlyjs();
plot(snip)
vline!(d["ang_impulsive"][2][i,:])
plot(abs.(hilbert(snip)))
chs = 1:size(rx_vect,2)

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
