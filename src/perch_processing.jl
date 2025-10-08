using CSV, DataFrames, Dates
using PlotlyJS
include("tabulate_data.jl")
include("audio.jl"); include("plotting.jl")
include("audacity.jl"); include("raven.jl")

# fname = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs/run_20250929_test4/ct.csv"
fname = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs_fake32kHz/run_20250929_test2/ct.csv"
# fname = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs_fake32kHz/run_20250929_test2/240214201545_ct2.csv"
# fname = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs_fake32kHz/run_20250929_test2/ct3_2.csv"
df = CSV.read(fname, DataFrame)
label_name = splitext(basename(fname))[1]
DEFAULT_fname2timestamp_func = fname2dt_soundtrap
df.timestamp = DEFAULT_fname2timestamp_func.(df.source_id)
real_fs=384_000; model_fs=32_000;
perch2_windowlen = 5; #seconds
real_windowlen = perch2_windowlen * model_fs / real_fs
df.offset_real = df.offset .* (model_fs/real_fs)

df_summary = combine(groupby(df, :source_id),
    :offset => mean,
    :offset => median,
    nrow => :count,
)
df_summary.timestamp = DEFAULT_fname2timestamp_func.(df_summary.source_id)
df_summary.offset = [x.offset_real for x in groupby(df, :timestamp)]
sort!(df_summary, :timestamp)

# check offset values uniqueness, shouldnt have overlapping detections
df_summary.count_unique = df_summary.offset .|> unique .|> length 
sum(df_summary.count_unique) != size(df,1) && 
    @warn "Some offset values are not unique, there might be overlapping detections! unique:counts $(sum(df_summary.count_unique)):$(size(df,1))"

# res_dir = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs/run_20250929_test4/labels"
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs_fake32kHz/run_20250929_test2/labels"
# res_dir = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs_fake32kHz/run_20250929_test2/labels2"
# res_dir = "/media/spin/anas2/data_res/dolphin/marecet/clustering/perch2/raw/datai/full_fs_fake32kHz/run_20250929_test2/labels3_2"

detectionsfiles2plot2(df_summary; res_dir=res_dir,
    latitude=6.3675, longitude=99.79774, add_daynight=true,
    trace_list = [:count => "counts", :count_unique => "counts_unique"])


mkpath(res_dir)
for row in eachrow(df_summary)
    @info "Timestamp: $(row.timestamp), Offset array: $(row.count)"
    # @info joinpath(res_dir, "perch-audacity_$(row.source_id).txt")
    audacity_label([row.offset row.offset .+ real_windowlen], joinpath(res_dir, "perch-audacity_$(row.source_id).txt"); prefix="perch_")
    raven_label([row.offset row.offset .+ real_windowlen], joinpath(res_dir, "perch-raven_$(row.source_id).txt"); prefix="perch_")
end


## click train selections from SB
outfolder = "/media/spin/anas2/data/marecet/labels/Porpoise clicks/clicktrain/clips_8338.240214201545"
labelfname = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/verify_20240116_edit/raven_8338.240214201545_t2_375081633500817_d200_cps40_0_train_only_removed.txt"
labels = CSV.read(labelfname, DataFrame)
dts = DEFAULT_fname2timestamp_func(labelfname)
data_table_fname = "/media/spin/anas2/data/marecet/datai/summary.csv"
df_data = CSV.read(data_table_fname, DataFrame)
aufname = findfirst( ==(dts), df_data.datetime) |> (i -> df_data.filepath[i])

extract_audio_segment(aufname, labels[:,["Begin Time (s)","End Time (s)"]], outfolder; win_size=1.0)

extract_segments_from_ravenlabels(labelfname, data_table_fname, outfolder)


## extract audio segments based on detections
res_dir_clips = joinpath(res_dir, "clips")
data_table_fname = "/media/spin/anas2/data/marecet/datai/summary.csv"
df_data = CSV.read(data_table_fname, DataFrame)
for row in eachrow(df_summary)
    @info "Timestamp: $(row.timestamp), Offset array: $(row.count)"
    # labels = DataFrame("Begin Time (s)" => row.offset, "End Time (s)" => row.offset .+ real_windowlen)
    dts = row.timestamp
    aufname = findfirst( ==(dts), df_data.datetime) |> (i -> df_data.filepath[i])
    outfolder = joinpath(res_dir_clips, Dates.format(dts, "yyyymmdd_HHMMSS"))
    mkpath(outfolder)
    extract_audio_segment(aufname, [row.offset row.offset .+ (model_fs * perch2_windowlen / real_fs)], outfolder; win_size=1.0)
end

df_summary