# using SignalAlignment
include("audio.jl")
include("dsp.jl")
using PlotlyJS
include("plotting.jl")

aufname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/Ball_3way_S2/250519_002_0001.WAV"
vidfname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/Ball_3way_S2/GH019986.MP4"

aufname = "/media/spin/anas2/data/calf/upload/Double_ball_S1/250403_002_0001.WAV"
vidfname = "/media/spin/anas2/data/calf/upload/Double_ball_S1/GH019888.MP4"


#~ get sync time between video and audio file #######################################
data, fs = readAudio(aufname)
data_v, fs_v = get_videos_audiodata(vidfname)

# @info "duration of audio file: $(size(data,1)/fs)s"
# @info "duration of video file: $(size(data_v,1)/fs_v)s"

data_down = resample(@view(data[1:fs*120,1]), fs_v/fs)
template_start = fs_v*60
template_end = fs_v*2*60
template_sig = data_down[template_start:template_end]
# max_time = fs2*60
delays, delay_conf = finddelay2(@view(data_v[:,1]), template_sig)
# @info "Delay between video and audio signal: $((template_start - delays)/fs2)s, confidence: $delay_conf"
@info "Delay of audio against video signal: $((delays - template_start)/fs_v)s, confidence: $delay_conf"
###################################################

#~ Organized split video and audio files
data_dir = "/media/spin/anas2/data/calf/upload/Double_ball_S1"
data_dir = "/media/spin/anas2/data/calf/new_2025_freeplay_report/Ball_3way_S2"
res_dir = "/media/spin/anas2/data_res/dolphin/calf/temp/delete"
occursin.( Ref(Regex(join(vidtypes, '|'))), readdir(data_dir))

vidlist = filter( x -> occursin(Regex(join(vidtypes, '|')), x|>lowercase), readdir(data_dir; join=true))
audlist = filter( x -> occursin(Regex(join(autypes, '|')), x|>lowercase), readdir(data_dir; join=true))

# combine all video files into one file
mkpath(res_dir)
temp_filelist = joinpath(res_dir, "temp_filelist.txt")
write(temp_filelist, join(["file '$v'" for v in vidlist], "\n"))

cmd = `ffmpeg -f concat -safe 0 -i $temp_filelist -c copy $res_dir/combined__$(join(basename.(vidlist), '_')).mp4`
print(cmd)
try
    @ffmpeg_env run(cmd)
    rm(temp_filelist)
catch e
    @error "FFmpeg command failed: $e"
end

find_vid_vs_audio_syncdiff_timesegment(vidlist[1], audlist[1]; flag_verbose=true)

#~ only need to run this
include("run_example.jl")

res_dir = "/media/spin/anas2/data_res/dolphin/calf/temp/delete"
res_dir = ""
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Double_ball_S1"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Double_ball_S1"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_4males_19_09_2025"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_4males_19_09_2025"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Calibration_19_09_2025"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Calibration_19_09_2025"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_3males_19_09_2025"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_3males_19_09_2025"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Double_ball_S1"; res_dir=res_dir)

####
#~ run analysis on split video and audio files
using NaturalSort
data_dir = "/media/spin/anas2/data/calf/Single_ball"
result_directory = "/media/spin/anas2/data_res/dolphin/calf/Single_ball/res_20250910"

data_folders = sort(readdir(data_dir; join=true); lt=natural) #sort(readdir(data_dir; join=true) |> filter(isdir) |> filter(x -> occursin("Single", x)); lt=natural)
ftype = ".flac"
mkpath(result_directory)

# tabulate data
summary_filepath = joinpath(result_directory, "summary.csv")
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024"; output=summary_filepath, filetype="wav", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
tabulate_data(data_dir; output=summary_filepath, filetype="flac", fname2timestamp_func=DEFAULT_fname2timestamp_func, flag_filepath=true)
# df = CSV.read(summary_filepath, DataFrame)
# df = sort(df, :datetime)
# CSV.write(summary_filepath, df)

for aufol in data_folders
    @info "--- Processing folder: $aufol"
    files = readdir(aufol; join=true) |> filter(endswith(ftype))
    res_dir = joinpath(result_directory, Dates.format(DEFAULT_fname2timestamp_func(files[1]), "yyyymmdd_HHMMSS"))
    isnothing(res_dir) || mkpath(res_dir); cp("web/index.html", joinpath(res_dir, "index.html"))
    for file in files
        fname, ext = splitext(basename(file))
        @info "Processing file: $(dirname(file))/\033[32m$fname\033[0m$ext"
        try
            detect_impulseNtonal(file, res_dir)
        catch e
            @error "Error processing file $file: $e"
        end
        GC.gc() # run garbage collector to free memory
    end
end














asdf
#~ old code below ##########################################################
template_window = (-fs2÷2+1:fs2÷2) #(-fs2÷50+1:fs2÷50)
max_time = argmax(d_z)


@info "Audio sync time: $(max_time/fs2)s"
win2 = template_window .+ max_time
template_sig = d_z[win2]
template_sig |> PlotlyJS.plot
# m = mfilter(template_sig, data2[:,1])


m = mfilter(template_sig, @view(data2[:,1]))
video_sync_time = argmax(abs.(m))
correl_sign = m[video_sync_time] |> sign

[data2[video_sync_time .+ (1:length(template_window))] template_sig*correl_sign] |> plot_norm
@info "Video sync time: $(video_sync_time/fs2)s"
@info "Audio sync time: $(max_time/fs2)s"
vid_aud_sync_diff = video_sync_time - max_time
@info "Video-Audio sync time difference: $(vid_aud_sync_diff/fs2)s"


aud_vid_sync_diff = audio_sync_time - max_time
@info "Audio-Video sync time difference: $(aud_vid_sync_diff/fs2)s"
@info "Audio-Video sync time difference(samples): $(aud_vid_sync_diff)"
win_test = (1:2000) .+ max_time
plot_norm([d_z[win_test .+ aud_vid_sync_diff] data2[win_test,1]])#; title="Test plot around video max correlation point")

########################################
# Further analysis
########################################
data2[]
d_z[argmax(m) .+ (-fs2÷25+1:fs2÷25)] |> PlotlyJS.plot
#  |> PlotlyJS.plot

delays = align_signals([d_z[1:fs2*5], data2[1:fs2*5,1]], Delay(delay_method=DTWDelay()))

plot(d_z[win2[1]+delays[1][1]:win2[1]+delays[1][1]+1000,1])
plot(data2[win2[1]:win2[1]+1000,1])

win = (23*fs + 87843):(23*fs + 91223)

# sum(data[win,3:5]; dims=2) |> plot
a=resample(data[win,1], fs2/fs)
win2 = 19*fs2+5000+58:19*fs2+5845+58
[a -data2[19*fs2+5000+58:19*fs2+5845+58,1]] |> plot_norm
d_z=resample(data[:,1], fs2/fs)

# sanity check
[sum(data[win,3:5]; dims=2) data[win,1]] |> plot_norm

(data2[:,1] - data2[:,2]) .|> abs |> mean # both channels of video audio are identical
( norm_max(sum(data[:,3:5]; dims=2)) - norm_max(data[:,1])) .|> abs |> mean # sum of last 3 channels of audio are identical to the first two channels