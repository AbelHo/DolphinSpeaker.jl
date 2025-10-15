include("detector_impulsive.jl")
include("detector_tonal.jl")
include("detector_band.jl")

using DelimitedFiles
function write_csv(fname, input_arr, delim=','; header=["datetime" "num_impulse" "num_impulseINtrain" "num_impulsetrain" "num_tonal" "num_tonalsegment" "num_noise" "filepath"])
    if !isfile(fname)
        writedlm(fname, header, delim)
    end

    if !Sys.islinux()
        open(fname, "a") do io
            writedlm(io, input_arr, delim)
        end
    else
        text = read(fname, String)
        @debug text
        open(fname, "w") do io
            write(io, text)
            writedlm(io, input_arr, delim)
        end
    end
end

# conv_onesided(u,v) = conv(u, v)[length(v):end-(length(v)-1)]
# conv_onesided_centreMID(u,v) = conv(u, v)[ (length(v)+1:end-(length(v)-1)) .- length(v)÷2]

#~ click detection
# aufname = "/Users/abel/Documents/data/Hawaii_2022-09/punnet_red/2022-09-27/2022-09-27_08.05.34__2-1664301934595.wav"
# res_dir = "/Users/abel/Documents/data_res/aspod/Hawaii_2022-09/punnet_red"
# threshold_impulsive = 0.01 #0.21 #0.215
# click_train_minlen = 7 #[n+1]clicks within (click_train_check_interval)seconds
# click_train_check_interval = 1 #seconds has at least (click_train_minlen) clicks
# rx = 0.5/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) )
# rx_vect = [real(rx); imag(rx); zeros(1,3)]

function detect_impulseNtonal_dir(aufname, res_dir; processed_skip_flag=false)
    if isdir(aufname)
        _,files = split_vid_au(aufname)
        return process_one_file.(files, Ref(res_dir); magnetic_declination_deg=9, overlay_flag=false, processed_skip_flag=processed_skip_flag)
    else
        return detect_impulseNtonal(aufname, res_dir; processed_skip_flag=processed_skip_flag)
    end
end

function detect_impulseNtonal(aufname::String, res_dir; kwargs...)
    data, fs, _, opt, timestamp = readAudio(aufname)
    detect_impulseNtonal((aufname, data, fs, timestamp), res_dir; opt=opt, kwargs...)
end
function detect_impulseNtonal(aufname_data_fs_timestamp::Tuple, res_dir;
    threshold_tonal = threshold_tonal, freq_maxbandwidth = freq_maxbandwidth, freq_width_db=freq_width_db,
    percent_quiet = percent_quiet, tonal_band_pass=tonal_band_pass,
    impulsive_band_pass=impulsive_band_pass,
    processed_skip_flag = false, opt=nothing,
    rx_vect=rx_vect, ref_channel=ref_channel,
    detect_impulse=detect_impulse, detect_tonal=detect_tonal, detect_boat=detect_boat, detect_impulsetrain=detect_impulsetrain,
    kwargs...
    )
    aufname, data, fs, timestamp = aufname_data_fs_timestamp
    @info "-------"*aufname
    fname = splitext(aufname)[1]*"_t"*string(threshold_impulsive)*"_d"*string(dist_impulsive) *"__cps"*string((click_train_minlen+1)/click_train_check_interval)*  ".jld2" |> basename
    @debug fname
    @debug "res_dir: ", res_dir
    isnothing(res_dir) || mkpath(res_dir)
    processed_skip_flag && isfile(joinpath(res_dir, fname)) && (@info("___skipped!..."); return nothing)
    
    # data, fs, _, opt, timestamp = readAudio(aufname)
    isempty(data) && (@warn("___skipped!^^^^^^^Empty Data..........."); return nothing)
    # res_impulse = detect_impulse((aufname, data, fs), res_dir; ref_channel=ref_channel, threshold=0.21)

    # band_pass = impulsive_band_pass #[60, fs/2*.98]
    data_filt = data;
    # if !iszero(band_pass[1]) || !isinf(band_pass[2])
    #     filter_type = nothing
    #     if iszero(band_pass[1]) 
    #         filter_type = Lowpass(band_pass[2]; fs=fs)
    #     elseif isinf(band_pass[2]) 
    #         filter_type = Highpass(band_pass[1]; fs=fs)
    #     else
    #         filter_type = Bandpass(band_pass[1], band_pass[2]; fs=fs)
    #     end

    #     filter_weight = digitalfilter(filter_type, Butterworth(butterworth_size))
    #     data_filt = mapslices( x -> filtfilt( filter_weight, x), data, dims=1)
    # end

    res_impulse = nothing
    res_impulsetrain = nothing
    task_impulsive = @async begin
        res_filt = detect_impulse((aufname, data_filt, fs), res_dir; ref_channel=ref_channel, kwargs...)
        res_impulse=res_filt
        if isnothing(threshold_impulsive)
            fname = splitext(aufname)[1]*"_t"*string(res_impulse.threshold)*"_d"*string(dist_impulsive) *"__cps"*string((click_train_minlen+1)/click_train_check_interval)*  ".jld2" |> basename
        end
        # tdoas = get_tdoa_raw(data_filt, res.pind_good ; window=window_impulsive, ref_channel=ref_channel)
        res_impulsetrain = detect_impulsetrain(res_filt, res_dir)
    end


    ###################################################################
    #~ whistles/tonal_indices
    # band_pass = [5000, 14001]
    # threshold_tonal = -95#without minusing noise 
    # band_pass = [5000 64000] #[5000 24000] #[5000 28000]
    # threshold_tonal = 24#26# -95#without minusing noise -93# -6# -68 #-90 #-95 #
    # freq_maxbandwidth = 10000 #300
    # freq_width_db=3
    # percent_quiet = 0.001

    res_tonal = nothing
    res_tonalsegment = nothing
    task_tonal = @async begin
        res_tonal = detect_tonal((aufname, @view(data[:,get_relevant_channels(rx_vect)]), fs), res_dir; 
                    ref_channel=ref_channel, thresh_tonal=threshold_tonal,
                    band_pass=tonal_band_pass,
                    freq_maxbandwidth=freq_maxbandwidth, freq_width_db=freq_width_db,
                    percent_quiet=percent_quiet)

        res_tonalsegment = combine_detections_conv(data, res_tonal; res_dir=res_dir)
    end

    res_boat = detect_boat((aufname,data,fs), res_dir)

    wait(task_impulsive)
    wait(task_tonal)
    # ang_impulsetrain, tdoas = detection2angle(data_filt, res_impulsetrain.pind_good, rx_vect;
    #             fs=fs, window=window_impulsive)#, getTDOA_func=get_tdoa_max)

    if !isnothing(res_dir)
        try #if !Sys.islinux()
            savejld(joinpath(res_dir, fname); aufname, timestamp,
                res_impulse= hasproperty(res_impulse, :data_filt) ? Base.structdiff(res_impulse, NamedTuple{(:data_filt,)}) : res_impulse,
                res_impulsetrain, res_tonal, res_tonalsegment, res_boat,
                threshold_impulsive, dist_impulsive, click_train_check_interval, fs
            )   #,ang_impulsetrain)
        catch err #else
            @info "..Cant directly save jld file"
            jldsave(joinpath("",fname); aufname, timestamp,
            res_impulse= hasproperty(res_impulse, :data_filt) ? Base.structdiff(res_impulse, NamedTuple{(:data_filt,)}) : res_impulse,
            res_impulsetrain, res_tonal, res_tonalsegment, res_boat,
            threshold_impulsive, dist_impulsive, click_train_check_interval, fs
            ) #,ang_impulsetrain)
            if isfile(joinpath(res_dir, fname)) 
                println("File already exist, replace? [y/n]")
                s = readline()
                if lowercase(s[1])=='y'
                mv( joinpath("", fname), joinpath(res_dir, fname), force=true)
                else
                    @warn "skip saving"
                end
            else
                mv( joinpath("", fname), joinpath(res_dir, fname))
            end
        end
    end
    
    if isnothing(timestamp)
        try
            timestamp = split(basename(aufname),"__")[1]
        catch err
            timestamp = basename(aufname)
        end
    end
    isnothing(res_dir) || write_csv(joinpath(res_dir,"counts.csv"), 
        [timestamp length(res_impulsetrain.pind) res_impulsetrain.num_click_in_trains res_impulsetrain.num_detection res_tonal.num_detection res_tonalsegment.num_detection res_boat.num_detection aufname])

    return (;res_impulse, res_impulsetrain, res_tonal, res_tonalsegment, click_train_minlen, res_boat, fname, fs)
end

"""
    detect_all(aufol::String, res_dir::String; ftype=".flac")

Detects impulses and tonal sounds in audio files within a specified directory and generates a summary of the detections.

# Arguments
- `aufol::String`: The directory containing the audio files to be processed.
- `res_dir::String`: The directory where the results and summary will be saved.
- `ftype::String`: The file type of the audio files to be processed (default is ".flac").

# Description
This function performs the following steps:
1. Creates a summary directory within the specified results directory.
2. Tabulates data from the audio files and saves it as a CSV file in the summary directory.
3. Detects impulses and tonal sounds in each audio file and saves the results in the specified results directory.
4. Generates a summary plot of the detections and saves it as an HTML file in the summary directory.

# Example
```julia
detect_all("/path/to/audio/files", "/path/to/results", ftype=".wav")
```
"""
function detect_all(aufol, res_dir; ftype=".flac")
    dir_summary = joinpath(res_dir, "summary")
    @info dir_summary
    mkpath(dir_summary)
    tabulate_data(aufol, output=joinpath(dir_summary,"files.csv"), fname2timestamp_func=fname2dt_soundtrap, filetype="flac")
    res = detect_impulseNtonal.(readdir(aufol; join=true) |> filter(endswith(ftype)), Ref(res_dir); detect_impulse=detect_impulseNarrowBand, processed_skip_flag=true)
    detectionsfiles2plot(joinpath(res_dir, "counts.csv"); res_dir=dir_summary , plottype=PlotlyJS.bar)
end