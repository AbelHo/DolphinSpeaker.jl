include("detector.jl")
include("plotting.jl"); plotlyjs()
include("tabulate_data.jl")

set_device__soundtrap()

# tabulate data
summary_filepath = "/media/spin/anas2/data/marecet/datai/summary.csv"
tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024"; output=summary_filepath, filetype="wav", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024"; output=summary_filepath, filetype="flac", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
df = CSV.read(summary_filepath, DataFrame)
df = sort(df, :datetime)
CSV.write(summary_filepath, df)

include("src/detector.jl")
include("src/plotting.jl"); plotlyjs()
include("src/tabulate_data.jl")
set_device__soundtrap()
# impulsive_autothreshold_median_ratio = 4
# threshold_impulsive = nothing
aufname = "/media/spin/anas2/data/marecet/pinger_survey/17012025_Sc2_H9/20250117_102741.WAV"
aufname = "/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3/8338.240215004543.wav"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/test_2025-06-05/h4"
clip_len = 40
min_freq = 80_000

aufname="/media/spin/anas2/data/marecet/Acoustics-Abel_TestDataBenchmark/SoundTrap - PAM & towed/Irrawaddy - towed/8745.240814100510.wav"
aufname = "/media/spin/anas2/data/marecet/Acoustics-Abel/SoundTrap - PAM & towed/Indo-Pacific humpback dolphins/14052024 - Sc6/7003.240514151859.wav"
res_dir="/media/spin/anas2/data_res/dolphin/marecet/temp/test3_ST_SC"
threshold_tonal = 5
#others
min_freq = 10_000

res = detect_impulseNtonal(aufname, res_dir);

data, fs = readAudio(aufname)
data = filter_simple(data, [1000, Inf]; fs=fs)
clips = extract_clips(data, res.res_impulse.pind_good, clip_len)

ffts = map(x-> abs2.(x) |> x->x[3:end], rfft.(clips, 1))
freqss = fs÷clip_len*2:fs÷clip_len:fs÷2
peak_freqs = argmax.(ffts) .|> x-> freqss[x]
# mean_freqs = map(x-> sum(x .* freqss)/sum(x), ffts)

# pows = map( x-> pow2db.(x), ffts)
# bandwidth_db = 0535


nbhf = findall( >(min_freq), peak_freqs)

res.res_impulse.pind_good[nbhf]
audacity_label(res.res_impulse.pind_good[nbhf] ./ fs, joinpath(res_dir, res.res_impulse.outfname *"__cps"*string((click_train_minlen+1)/click_train_check_interval)*  "_nbhf_$min_freq.txt" |> basename))
raven_label(res.res_impulse.pind[nbhf] ./ fs, joinpath(res_dir, "raven_" * res.res_impulse.outfname *"__cps"*string((click_train_minlen+1)/click_train_check_interval)*  "_nbhf_$min_freq.txt" |> basename))
            

res3 = detect_impulsetrain(res.res_impulse.pind_good[nbhf], res.res_impulse.pind_good[nbhf] ./ fs, res.res_impulse.ppeak[nbhf], res.res_impulse.len_data, res.res_impulse.outfname, res_dir)

res4 = detect_impulsetrain(res.pind_good, res.pind_good ./ fs, res.ppeak, res.res_impulse.len_data, res.res_impulse.outfname, res_dir)


function get_bandwidth_db(pow, freqss, db=3)
    # pows = map( x-> pow2db.(x), ffts)
    # maxi = map( x-> maximum(x), pows) |> argmax
    maxi = argmax(pow)
    count(findall( >(maxi-db), pow) ) * freqss[2]
end

# function get_min_baseline()

############################################################################################################
#~ temp single analysis, NP click correlation factor
aufname = "/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard4/8338.240119124455.wav"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250723"
set_device__soundtrap()
data, fs = readAudio(aufname; fname2timestamp_func=fname2dt_soundtrap)
data_filt = filter_simple(data, [1000 Inf]; fs=fs)
window_extract = [-50 205]
clips = extract_clips(data_filt, window_extract .+ res.res_impulse.pind_good, NaN; flag_matrix=true)
# clips = extract_clips(data, res.res_impulse.pind_good, length(window_impulsive); flag_matrix=true)
clips = hcat(res.res_impulse.clips[res.res_impulse.nbhf]...); ffts = hcat(res.res_impulse.ffts[res.res_impulse.nbhf]...); peak_freqs = hcat(res.res_impulse.peak_freqs[res.res_impulse.nbhf]...)
selections = 1:size(clips,2) #29800:29950 #29867:29880 #35926:36140 #35636:35651 #35359:35400
plot(clips[:,selections], labels=reshape(string.(selections),1,length(selections)) )
plot(res.res_impulse.freqss, ffts[:,selections], labels=reshape(string.(selections),1,length(selections)) )
# selections = 35388:35399
# c = clips[:, [35388, 35393, 35399]]# |> plot
c = clips[:, selections]

# correl = map(x-> mfilter(c[:,i], x),  eachcol(c)) .|> energy
correls = map(c|>eachcol) do ref
    # map(x-> mfilter(ref, x) |> hilbert .|> abs,  eachcol(c)) .|> sum |> sum
    map(x-> mfilter(ref, x),  eachcol(c)) .|> energy |> sum
end;
thresh = 1.3e12 #7e14#good #3e14
plot(correls); hline!([thresh], label="Threshold")
selections[findall( >(thresh), correls) ] |> show
[1:length(selections) selections] |> showall

selections = selections[findall( >(thresh), correls) ]# |> show
# energy.(r) |> plot!
res.res_impulse.pind_good[selections] |> diff |> plot

# threshold over peak value
thresh = 30
selections = selections[findall( >(thresh), res.res_impulse.ppeak[selections] )]


############################################################################################################
#~ clips frequency analysis and correlation
aufname = "/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard3/8338.240426130000.flac"
respath = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/20240418_100000/8338.240426130000_t1.7870080480686543_d200__cps40.0.jld2"
aufname = "/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard3/8338.240419050000.flac"
respath =  "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/20240418_100000/8338.240419050000_t1.9312512448952104_d200__cps40.0.jld2"
aufname = "/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard3/8338.240418120000.flac"

target = "2025-05-03T15:26" #"2024-05-18T11:00" #"2024-05-21T13:00"# "2024-01-19T15:44"
summary_fname = "/media/spin/anas2/data/marecet/datai/summary.csv"
result_directory = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648"

result = find_closest_row(summary_fname, target)
aufname = result.filepath

# readdirjoin(x) = readdir(x; join=true)
respath = readdir(result_directory; join=true) |> 
    filter(isdir) .|> readdirjoin .|> #x->readdir(x; join=true) .|> #
    filter(endswith(".jld2")) .|>
    filter(contains(splitext(basename(aufname))[1])) |>
    filter(!isempty) |> first
    
res = load(respath)
res = dict2namedtuple(res)
# get clips
data, fs,_,_,timestamp = readAudio(aufname; fname2timestamp_func=fname2dt_soundtrap)
# data_filt = filter_simple(data, nb_bandpass; fs=fs) # nb_bandpass=[1000, Inf]
data_filt = filter_simple(data, impulsive_band_pass; fs=fs) # nb_bandpass=[1000, Inf]
window_extract = [-25 102] #[-50 205]
clips_fixed = extract_clips(data_filt, window_extract .+ res.res_impulsetrain.pind_good, NaN; flag_matrix=true)
# clips = extract_clips(data_filt, window_extract .+ res.res_impulse.pind_good, NaN; flag_matrix=true)

# extract any len signals
clips = extract_clips(data_filt, [res.res_impulsetrain.train_start res.res_impulsetrain.train_end], NaN; flag_matrix=false)

# plot clips
clips_plot_dir = "temp/clips_train_" * Dates.format(timestamp, "yyyymmdd_HHMMSS")
mkpath(clips_plot_dir)
# clips = filter_simple.(clips, Ref([80_000, Inf]); fs=fs) # clips = clips .|> norm_max
train_start_ind = [res.res_impulsetrain.train_start_ind... length(res.res_impulsetrain.pind_good)+1]
nfunc(x) = sqrt(sum(abs2.(x)))
threshold_autocor = 4; threshold_n_autocor = 3; autocor_clips = Int[]
for i in eachindex(clips)
    a = plot(signal(clips[i],fs); title=string(i))#, xlims=(1, length(clips[i])))
    b = specgram(clips[i]; fs=fs, colorbar=nothing, nfft=128, crange=80)#, crange=50, downsample=nothing)
    c = specgram(clips[i]; fs=fs, colorbar=nothing, nfft=round(Int,fs*.01)|>nextfastfft )
    # selections = res.res_impulsetrain.train_start_ind[i]:res.res_impulsetrain.train_start_ind[i+1]-1
    selections = train_start_ind[i]:train_start_ind[i+1]-1
    snip = clips_fixed[:, selections]
    d = plot_time_fft(snip, fs; legend_position=:outerbottom, labels=reshape(string.(selections),1,length(selections)))#; legend_position=:outerright)
    # plot(a, b, c, d; size=(1000, 600), layout=(4,1))#, legend=nothing)# |> display

    correls = map(snip|>eachcol) do ref
        map(x-> mfilter(norm_max(ref; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy
    end
    # e = Plots.heatmap(selections, selections, hcat(correls...))

    e = Plots.bar(selections, map(x-> mfilter(norm_max(x; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy)
    Plots.hline!([threshold_autocor], label="Threshold", color=:red)
    # plot(a, b, c, e, d; size=(1000, 1000), layout=(5,1))# |> display#, legend=nothing)# |> display
    
    savefig(joinpath(clips_plot_dir, "clip_$(i).html"))

    if count(map(x-> mfilter(norm_max(x; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy .> threshold_autocor) > threshold_n_autocor
        push!(autocor_clips, i)
    end
end
autocor_clips |> show
make_clip_index_html(clips_plot_dir; outname=basename(clips_plot_dir)*".html")

PlotlyJS.heatmap(z=spectrogram(clips[i][:], 256, 128; nfft=256, fs=fs).power) |> PlotlyJS.plot
##################
for i in 1:20
    winlen = 128 + 4*i      # window length varies from 132 to 208
    hop = 32 + 2*i          # hop size varies from 34 to 72
    fsval = fs * (0.95 + 0.01*i)  # fs varies slightly
    clipidx = (i <= length(clips)) ? i : 1  # fallback to 1 if not enough clips
    PlotlyJS.heatmap(z=spectrogram(clips[clipidx][:], winlen, hop; nfft=256, fs=fsval).power) |> PlotlyJS.plot |> display
end
##################

fftval, freqss = audiofft(clips, fs)
plot_time_fft(clips, fs; legend_position=:outerright)#  size=(1000, 400), legend=:outerright)
waterfall_3d(clips)
waterfall_3d(fftval, freqss)
plot(clips, legend=:outerright, legendfontsize=10, legend_title="Clip Index", legend_position=:right, legend_margin=30)

for i = 1:length(res.res_impulsetrain.train_start_ind)-1
    snip = clips[:, res.res_impulsetrain.train_start_ind[i]:res.res_impulsetrain.train_start_ind[i+1]-1]
    # snip_fft, freqss = audiofft(snip, fs)
    plot_time_fft(snip, fs; legend_position=:outerright, title="Train $i") |> display
    # waterfall_3d(snip, freqss)
    # waterfall_3d(snip_fft, freqss)
end
snip = clips[:, res.res_impulsetrain.train_start_ind[end]:end]
plot_time_fft(snip, fs; legend_position=:outerright, title="Train $(length(res.res_impulsetrain.train_start_ind))") |> display

i=32; snip = clips_fixed[:, res.res_impulsetrain.train_start_ind[i]:res.res_impulsetrain.train_start_ind[i+1]-1]
nfunc(x) = sqrt(sum(abs2.(x)))
nfunc(x) = energy(x)
nfunc(x) = 1
nfunc(x) = maximum(abs.(x); dims=1)
correls = map(snip|>eachcol) do ref
    map(x-> mfilter(norm_max(ref; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy
end
hcat(correls...) |> plot
hcat(correls...) |> Plots.heatmap
hcat(correls...) .|> log10 |> Plots.heatmap
plot_time_fft(snip, fs; legend_position=:outerright)

diff.(res.res_impulsetrain.train_impulse_list)

####
aufname = "/media/spin/anas2/data/marecet/labels/Porpoise clicks-20250701T064734Z-1-001/Porpoise clicks/03052025/Np6/8746.250503152639.wav"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_$(Dates.format(now(),"yyyymmdd"))/labeled/$(Dates.format(now(),"HHMMSS"))"

## no animals
aufname = "/media/spin/anas2/data/marecet/labels/NoAnimals/2025-08-05/8746.250805123209.wav"
aufname = "/media/spin/anas2/data/marecet/labels/NoAnimals/2025-08-05/8746.250805131708.wav"
#~ NBHF
res = detect_impulseNtonal(aufname, res_dir; detect_impulse=detect_impulseNarrowBand);

# 20250724
ftype = ".wav"
aufol = "/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/res_20250724/deployment_18012024_18042024"
res2 = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand, processed_skip_flag = true)
aufol = "/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard4"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/res_20250724/deployment_18012024_18042024"

ftype = ".flac"
aufol = "/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard2"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/res_20250724/deployment_20032024_30052024"
res2 = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand)
aufol = "/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard3"
res2 = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand)



aufol = "/media/spin/Extreme SSD/data/marecet/dataidata/deployment1_18012024_18042024/SDcard3"
res_dir = "/media/spin/anas/data_res/dolphin/marecet/datai_3"

aufol = "/media/spin/anas/data/marecet/datai/deployment_20032024_30052024/sdcard2"
res_dir = "/media/spin/anas/data_res/dolphin/marecet/datai/20032024_30052024"
ftype = ".flac"

## doing
aufol = "/media/spin/anas1/data/marecet/datai/deployment_20032024_30052024/sdcard3"
res_dir = "/media/spin/anas1/data_res/dolphin/marecet/datai/2024-04-18"
ftype = ".flac"

aufol="/media/spin/anas1/data/marecet/datai/deployment_18012024_18042024/SDcard3"
res_dir = "/media/spin/anas1/data_res/dolphin/marecet/datai/2024-02-14_test"

#doing
aufol = "/media/spin/One Touch/data/marecet/dataidata/deployment_20032024_30052024/sd4"
res_dir = "/media/spin/anas1/data_res/dolphin/marecet/datai/2024-03-20"
ftype = ".flac"
########################################################################################
aufol = "/media/spin/anas/data/marecet/datai/deployment1_18012024_18042024/SDcard4"
res_dir = "/media/spin/anas/data_res/dolphin/marecet/datai/18012024_18042024"
ftype = ".wav"

##### hydromoth perlis
ftype = ".WAV"
set_device__hydromoth()
aufol = "/media/spin/anas2/data/marecet/perlis/Mardiana - Soundscape Hydromoth/October 2021/Perlis/H2 17102021"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/perlis/fyp_t2/$(basename(aufol))"
# global DEFAULT_fname2timestamp_func = fname2dt_date_time
res2 = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulse);

aufol = "/media/spin/anas2/data/marecet/perlis/Mardiana - Soundscape Hydromoth/October 2021/Perlis/H2 20102021"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/perlis/fyp/$(basename(aufol))"
res2 = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand)
aufol = "/media/spin/anas2/data/marecet/perlis/Mardiana - Soundscape Hydromoth/October 2021/Perlis/H4 20102021"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/perlis/fyp/$(basename(aufol))"
res2 = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand)

##
aufol = "/home/spin/Documents/data/marecet"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/temp/hydromoth/perlis"

res_dir = "/home/spin/Documents/data/marecet/res2"
ftype = ".WAV"
global impulsive_autothreshold_median_ratio = 3
res = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulse);


########################################################################################

res2 = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand)
dir_summary = joinpath(res_dir, "summary")
mkpath(dir_summary)
plot_detection_summary(joinpath(res_dir, "counts.csv"); res_dir=dir_summary , filetype=".html")

aufol = "/media/spin/anas/data/marecet/hydromoth_flac"
res_dir = "/media/spin/anas/data_res/dolphin/marecet/hydromoth_flac"
ftype = ".flac"
DEFAULT_fname2timestamp_func = aufname -> fname2dt_date_time(aufname; skip_front=3)

dir_summary = joinpath(res_dir, "summary")
mkpath(dir_summary)
tabulate_data(aufol, output=joinpath(dir_summary,"files.csv"), fname2timestamp_func=fname2dt_soundtrap, filetype="flac")
res = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand, processed_skip_flag=true)
detectionsfiles2plot(joinpath(res_dir, "counts.csv"); res_dir=dir_summary , plottype=PlotlyJS.bar)
# plot_detection_summary( filetype=".html")

# plot_recordings_summary("/media/spin/anas/data/megafauna/data_index/summary_CJ_N_1.csv", "/media/spin/anas/data/megafauna/data_index/plots")

# tonal segment length
segments = [(res[i].res_tonalsegment.train_end - res[i].res_tonalsegment.train_start) ./ fs for i in eachindex(res)]
segments .|> length

# analysis ici of each train
"""
    Process and visualize impulse train data from a loaded result file.

This code performs the following steps:
- Loads a result file containing impulse train data.
- Initializes arrays to store inter-click intervals (`ici`) and normalized maximum histogram weights (`max_weight`).
- Iterates over each impulse train in the data:
    - Detects impulse groups using `detect_impulsegroup`.
    - Fits a histogram to the positive values of the detected impulses, with bins ranging from 0 to the maximum value in steps of 300.
    - Sets the second histogram bin's weight to zero to ignore it.
    - Plots the histogram as a bar plot.
    - Determines the most frequent inter-click interval (`ici`) as the bin edge with the highest weight.
    - Normalizes the maximum histogram weight by the number of impulses in the train and stores it.
    - Collects histogram weights and edges for further visualization.
- Plots a scatter plot of the most frequent inter-click intervals (`ici`), with point transparency (`alpha`) proportional to the normalized maximum weight.
- Displays bar plots for each histogram of impulse intervals.

# Arguments
- None directly; the code expects a specific file structure and the presence of `detect_impulsegroup` and `Histogram` fitting.

# Requirements
- The file at the specified path must exist and contain the expected data structure.
- The functions `detect_impulsegroup`, `fit`, and the `Histogram` type must be available.
- The `Plots` package must be loaded for visualization.

# Output
- Visualizations of the inter-click intervals and their distributions for each impulse train.
"""
res = load("/media/spin/anas/data_res/dolphin/marecet/datai_4/8338.240216013039_t2.332525802366085_d200__cps40.0.jl")
ind = 245 #240
ici = Array{Int}(undef,length(res["res_impulsetrain"].train_impulse_list))
max_weight = Array{Float64}(undef,length(res["res_impulsetrain"].train_impulse_list))
h_weights = []; h_edges = [];
for ind=eachindex(res["res_impulsetrain"].train_impulse_list)
    a = res["res_impulsetrain"].train_impulse_list[ind] |> detect_impulsegroup
    h = fit(Histogram, filter(>(0), a[:]), 0:300:maximum(a))
    h.weights[2] = 0
    Plots.bar(h.edges, h.weights)
    ici[ind] = h.edges[1][argmax(h.weights)]
    max_weight[ind] = maximum(h.weights) / length(res["res_impulsetrain"].train_impulse_list[ind])

    push!(h_weights, h.weights); push!(h_edges, h.edges)
end
# num_clicks = length.(res["res_impulsetrain"].train_impulse_list)
Plots.scatter(ici; alpha=max_weight/maximum(max_weight))#, label="num_clicks")
for ind=eachindex(h_edges)
    Plots.bar(h_edges[ind], h_weights[ind]; title=(string(ind))) |> display
end



## temp
aufol = "/home/spin/Documents/data/marecet"
aufname = joinpath(aufol, "H2_20220608_190000_0535.flac")
global impulsive_autothreshold_median_ratio = 2
ftype = ".flac"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/temp/hydromoth/test2"
res = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulse);


data, fs = readAudio(aufname)
sets = windowing.(Ref(data), res.res_impulse.pind_good[27617 .+ collect(0:2:4)], Ref(-20:107))

indices = 22601:22700 #27617 .+ collect(0:2:4)
sets = windowing.(Ref(data), res.res_impulse.pind_good[indices], Ref(-20:107))
plot(indices, map(x->x .* sets[1] |> sum, sets) )

plot_time_fft.(sets, Ref(fs)) .|> display

res.res_impulsetrain.pind_good |> diff |> plot
res.res_impulse.pind_good |> diff |> plot