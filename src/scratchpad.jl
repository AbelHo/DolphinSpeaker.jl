include("audio.jl")
include("aspod.jl")
using VideoIO
include("readImages.jl")
include("utils.jl")

aufname = "/Users/abel/Documents/data/calf/coop/20231128/20231128_15.16.48_log.flac"
vidfname = "/Users/abel/Documents/data/calf/coop/20231128/20231128_15.16.48_log.mkv"
res_dir = "/Users/abel/Documents/data_res/calf/coop/2023-11-28"

data, fs, _,_, timestamp = readAudio(aufname)

dets = process_detections(aufname, vidfname; res_dir=res_dir)

d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t116.9741055601481_d15000__cps0.375.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t116.9741055601481_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.12.51_log_t12_d15000.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t12_d15000.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t116.15144729819004_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.18.52_log_t116.15144729819004_d15000__cps0.375.jld2"])
d = load(["/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t12_d15000.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t116.14395330535885_d15000__cps0.375_angles.jld2", "/Users/abel/Documents/data_res/calf/click_test/20211214__get_tdoa_raw_MaxEnergyRefChannel/20211214_12.24.53_log_t116.14395330535885_d15000__cps0.375.jld2"])
plot_ang(d["res_impulsetrain"], d["ang_impulsive"][1][1]; label=["azimuth" "inclination"], type="Power")

d = load(["/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev/20231128_15.16.48_log_t4_d800.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev/20231128_15.16.48_log_t46.36293944139811_d800__cps60.0_angles.jld2", "/Users/abel/Documents/data_res/calf/coop/2023-11-28/more_dev/20231128_15.16.48_log_t46.36293944139811_d800__cps60.0.jld2"])

data_filt = filter_simple(data[:,1:size(rx_vect,2)], impulsive_band_pass; fs=fs)
data_filt = filter_simple(data[:,1:size(rx_vect,2)], [5_000 180_000], []; fs=fs)

i = 301#400
i = 1164
plot_one_event(res_dir, vidfname, d["res_impulsetrain"].pind_good_inS, data_filt, d["res_impulsetrain"].pind_good, d["ang_impulsive"][2], window_impulsive, d["ang_impulsive"][1][1]; i=i, plotfunc=gr)


snip = data_filt[window_impulsive .+ d["res_impulsetrain"].pind_good[i],1:size(rx_vect,2)]
plot_time_fft(snip, fs)
plotlyjs();
plot(snip)
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
extrema_and_indices(snip)
# get_tdoa_raw(data_filt, d["res_impulsetrain"].pind_good[i] ; window=window_impulsive, ref_channel=ch)

