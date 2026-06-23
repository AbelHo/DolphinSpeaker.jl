using ProgressMeter
include("detector.jl")
include("plotting.jl"); plotlyjs()
include("tabulate_data.jl")

set_device__soundtrap()

#~ tabulate raw data
data_folders = 
["/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard4",
"/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3",
"/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard2",
"/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard3",
"/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard4",
"/media/spin/anas2/data/marecet/datai/deployment_3_27022025_090702025",
]
# data_folders = 
# ["/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3/wav"]
result_directory = "/media/spin/anas2/data_res/dolphin/marecet/datai/res_$(Dates.format(now(),"yyyymmdd"))/raw/$(Dates.format(now(),"HHMMSS"))"
# result_directory = "/media/spin/anas2/data_res/dolphin/marecet/datai/res_20260619/raw/012437"

ftype = ".flac"
# tabulate data
# summary_filepath = "/media/spin/anas2/data/marecet/datai/summary.csv"
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024"; output=summary_filepath, filetype="wav", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024"; output=summary_filepath, filetype="flac", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
# df = CSV.read(summary_filepath, DataFrame)
# df = sort(df, :datetime)
# CSV.write(summary_filepath, df)

# debug only one file
# ENV["JULIA_DEBUG"] = "detector"



# analyse
for aufol in data_folders
    @info "Processing folder: $aufol"
    files = readdir(aufol; join=true) |> filter(endswith(ftype))
    files = vcat(files, readdir(aufol; join=true) |> filter(endswith("wav")))
    res_dir = joinpath(result_directory, Dates.format(DEFAULT_fname2timestamp_func(files[1]), "yyyymmdd_HHMMSS"))
    isnothing(res_dir) || mkpath(res_dir); isfile(joinpath(res_dir, "index.html")) || cp("web/index.html", joinpath(res_dir, "index.html"))
    @showprogress for file in files
        fname, ext = splitext(basename(file))
        @info "Processing file: $(dirname(file))/\033[32m$fname\033[0m$ext"
        try
            detect_impulseNtonal(file, res_dir; detect_impulse=detect_impulseNarrowBand, processed_skip_flag=true, processed_skip_strict=false)
        catch e
            @error "Error processing file $file: $e"
        end
        GC.gc() # run garbage collector to free memory
    end
end



res_dir = result_directory
pp = detectionsfiles2plot2.(
    joinpath.( filter(isdir,readdir(res_dir; join=true)), "counts.csv");
    res_dir=res_dir,
    plottype=PlotlyJS.bar, add_daynight=true);



# detectionsfiles2plot("/media/spin/anas2/data_res/dolphin/marecet/datai/res_20250724/deployment_18012024_18042024/counts.csv"; res_dir="/media/spin/anas2/data_res/dolphin/marecet/datai/res_20250724/deployment_18012024_18042024", plottype=PlotlyJS.bar)


#~ try just one
using Images, Printf
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_$(Dates.format(now(),"yyyymmdd_HHMMSS"))"
aufname = "/media/spin/anas2/data/marecet/datai/2025/03052025_NP6/8746.250503152639.wav"
res = detect_impulseNtonal(aufname, res_dir; detect_impulse=detect_impulseNarrowBand)

data, fs = readAudio(aufname)
data_filt = filter_simple(@view(data[:,get_relevant_channels(rx_vect)]), impulsive_band_pass; fs=fs)
data_filt = data

# res.res_impulse.pind_good_inS|>diff|>Plots.plot
# Plots.hline!([0.5 1.0 2.0]; labels=["0.5s" "1.0s" "2.0s"])
# cluster into groups when there is 1s gap
threshond_clusterseperation = 1.0  # seconds
# click_train_minlen = 3

t = res.res_impulse.pind_good_inS; t_inS=t
dt = res.res_impulse.pind_good_inS |> diff
sep = findall(dt .> threshond_clusterseperation)
clus_starts = [1; sep.+1]
clus_ends = [sep ; length(res.res_impulse.pind_good_inS)]
windows_ct = [t[clus_starts] t[clus_ends]] #.|> (x -> round.(Int, res.res_impulse.pind_good_inS[x] .* fs) .+ 1)
windows = res.res_impulse.pind_good .+ [window_impulsive[1] window_impulsive[end]]
c_in_ct = Array{Array{Int}}(undef, size(windows_ct, 1))
Threads.@threads for i in 1:size(windows_ct, 1)
    c_in_ct[i] = findall(x -> x .>= windows_ct[i, 1] .&& x .<= windows_ct[i, 2], t)
end
# remove less than click_train_minlen
c_in_ct = filter(x -> length(x) >= click_train_minlen, c_in_ct)

# extract clicks
data_filt = filter_simple(data, impulsive_band_pass; fs=fs); data=nothing;
clipsm = extract_clips_matrix(data_filt, windows)

# clipsm
if size(data,2) > 1
    clips_max = mapblocks(clipsm; dims=3) do x
        # @info size(x)
        chan = mapslices(energy, x; dims=1) |> vec |> argmax
        x[:, chan]
    end
else
    clips_max = clipsm[:, 1, :]
end
ffts, freqss = compute_rfft(clips_max, fs)#; type=:log)
ffts

#~ clustering
freqss=freqss[5:end]; ffts=ffts[5:end, :]
include("hdbscan.jl")
this_window=nothing; ct_ind=nothing
clusttering_method = :HDBSCAN # :OPTICS
rgbs_alpha_offset = 0.2 # for color clicks plot

outdir = "temp/$(clusttering_method)_ct2_minClus3_FFTnTime_colorclicks_datai-fullfreq3-highp"
mkpath(outdir)
open( joinpath(outdir, "$(clusttering_method)_clusters.csv"), "w") do io
    writedlm(io, ["c_ind" "time_inS" "ct_ind" "cluster_label"], ',')
end
# variable don't exist in this scope, define here
if !@isdefined rgb_bands
    freq_split = fs/2/3
    rgb_bands=[[0,freq_split], [freq_split,freq_split*2], [freq_split*2, fs/2]]
end

for (ct_idx,now_window) in enumerate(c_in_ct)
    ct_ind = ct_idx
    this_window = now_window
    include("cluster_plot.jl")
    CSV.write(joinpath(outdir, "$(clusttering_method)_clusters.csv")),
        DataFrame(c_ind=this_window, time_inS=t[this_window], ct_ind=ct_ind, cluster_label=labels;
        append=true)
    open( joinpath(outdir, "$(clusttering_method)_clusters.csv"), "a") do io
        writedlm(io, [this_window t[this_window] (ones(size(this_window)) .* ct_idx) labels], ',')
    end
    # ct_ind==2 && break
end