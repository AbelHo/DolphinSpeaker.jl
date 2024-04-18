using Peaks, SignalAnalysis, SignalAnalysis.Units, JLD2, DSP
using Statistics
include("audacity.jl")
include("config.jl") #impulsive_autothreshold_median_ratio
include("audio.jl")
include("dsp.jl")
""" requires the following default parameters to be set in config.jl
ref_channel
res_dir
dist_impulsive
window_impulsive
threshold_impulsive
"""

"""
    findPings(fname::String; save_fn=res_dir*fname[end-22:end]*"_label.jld2", ref_channel=1::Int, dist=nothing, channels=5, fs=400000)
    findPings(fname::AbstractArray; ref_channel=1::Int, dist=nothing)
Find peaks in the the signal by specifying a reference channel
# Examples
# For 1 second interval pinger & chatterbox XL with 5 channels 
```julia-repl
julia> peaks, peak_times, peak_vals = findPings("/somepath/filename.bin", save_fn=nothing, dist=Int(400000*0.9), channels=5)
1
```
"""
function findPings(fname::String; save_fn=res_dir*fname[end-22:end]*"_label.jld2", ref_channel=1::Int, dist=nothing, channels=5, fs=400000)
    if !isnothing(save_fn) && isfile(save_fn)
        peaks=readAudio(save_fn) #peaks=load(save_fn)
        peaks = peaks["peaks"]
    else   
        data, fs = load(fname) #FIXME to be more general with filetype #data=loadDataBin(fname, channels=channels)
        peaks = findPings(data, ref_channel=ref_channel, dist=dist)
        if !isnothing(save_fn)
            save(save_fn, "peaks", peaks, "peak_times", peaks[1]./fs, "peak_vals", peaks[2]);
        end
    end
    println(string(length(peaks[1]))*" pings")
    return peaks, peaks[1]./fs, peaks[2]
end


function findPings(data::AbstractArray; ref_channel=1::Int, dist=nothing, flagabs=DEFAULT_FLAGABS_findPings)
    if isnothing(dist)
        peaks = findmaxima(data[:,ref_channel])
    else 
        @debug size(data)
        @debug ref_channel
        @debug dist
        d = data[:,ref_channel]
        flagabs && (d = abs.(d))
        peaks = findmaxima(d, dist)
    end
    return peaks
end

function findNauda(fname::String, outFname=res_dir*fname[end-22:end]*"_label.txt"::String, fs=400000 ) #
    peaks,_,_ = findPings(fname; dist=Int(fs*.9))
    if !isempty(outFname)
        audacity_label(peaks[1]./fs, outFname)
        println(outFname)
    end
    return peaks
end

function detect_impulse(aufname::String, res_dir=nothing;  band_pass=impulsive_band_pass, kwargs...)
    @info basename(aufname)
    data, fs, _, opt, timestamp = readAudio(aufname)
    
    # if band_pass[2] == Inf
    #     filter_weight = digitalfilter(Highpass(band_pass[1], fs=fs), Butterworth(butterworth_size))
    # else
    #     filter_weight = digitalfilter(Bandpass(band_pass[1], band_pass[2], fs=fs), Butterworth(butterworth_size))
    # end
    # data_filt = filtfilt( filter_weight, data );

    detect_impulse((aufname, data, fs, timestamp), res_dir; kwargs...)#, data_filt, fs
end

function detectImpulseNangle(aufname::String, res_dir=nothing; 
    rx_vect=rx_vect)
    res_impulse, data_filt, fs = detect_impulse(aufname, res_dir)
    angs,tdoas = detection2angle(data_filt, res_impulse.pind_good, rx_vect; fs=fs)
     jldsave(joinpath(res_dir, basename(aufname) *"_.jld"); res_impulse, angs, tdoas)
    return (;res_impulse, angs, tdoas)
end


function detect_impulse(aufname_data_fs::Tuple, res_dir=nothing; band_pass=impulsive_band_pass,
    ref_channel=ref_channel, dist=dist_impulsive, window=window_impulsive, threshold=threshold_impulsive,
    return_datafilt=false)

    aufname, data, fs = aufname_data_fs
   
    @info ("Duration: " * string(size(data,1)/fs) *"seconds")
    @debug "filtering......."
    data_filt = filter_simple(data[:,1:size(rx_vect,2)], band_pass; fs=fs)
    @debug "hilberting...."
    # data_hil = data_filt|>hilbert.|>abs
    data_hil = data_filt[:,ref_channel]|>hilbert.|>abs
    @debug "finding peaks...."
    # pind, ppeak_all = findPings(data_hil; ref_channel=ref_channel, dist=dist)
    pind, ppeak_all = findPings(data_hil; ref_channel=1, dist=dist)

    if isnothing(threshold)
        # threshold = median(data_hil)*impulsive_autothreshold_median_ratio
        threshold = quantile(data_hil, .75) * impulsive_autothreshold_median_ratio
        @info "__auto thresholding: " * string(threshold)
    end
    threshold_indices = findall(>(threshold), ppeak_all)
    pind_good = pind[threshold_indices]
    pind_good_inS = (pind_good.-1) ./fs
    num_detection = length(pind_good_inS)
    ppeak = ppeak_all[threshold_indices]
    @info "Num of Clicks Detected: " * string(num_detection)

    outfname = ""
    if res_dir isa String
        isdir(res_dir) || mkpath(res_dir)
        outfname = splitext(aufname)[1]*"___t"*string(threshold)*"_d"*string(dist) |> basename
        audacity_label(pind_good./fs, joinpath(res_dir,outfname*".txt"))
    else
        @warn "detection not saved to file, no result folder(res_dir) given"
    end

    if return_datafilt
        return (;pind, ppeak, ppeak_all, pind_good, pind_good_inS, threshold_indices, threshold, dist, aufname, num_detection, len_data=size(data,1), outfname, fs, data_filt)
    else
        return (;pind, ppeak, ppeak_all, pind_good, pind_good_inS, threshold_indices, threshold, dist, aufname, num_detection, len_data=size(data,1), outfname, fs)
    end
end

conv_onesided(u,v) = conv(u, v)[length(v):end-(length(v)-1)]

function detect_impulsetrain(res, res_dir=res_dir; 
    click_train_minlen=click_train_minlen, click_train_check_interval=click_train_check_interval)
    
    if length(res.pind_good_inS) > click_train_minlen
        cl_minlen = conv_onesided(diff(res.pind_good_inS), ones(click_train_minlen,)) #conv(diff(res.pind_good_inS), ones(click_train_minlen,))[click_train_minlen:end-(click_train_minlen-1)]
        cltrain_start_inds = cl_minlen .< click_train_check_interval
        click_inds = conv(cltrain_start_inds, ones(click_train_minlen+1,)) .> .5
    else
        click_inds = []
    end

    click_accum = Int[]
    click_in_trains = Int[]
    train_switch = false
    train_count = 0
    train_start = Int[]
    train_end = Int[]
    train_start_ind = Int[]
    for i in 1:length(click_inds)
        if click_inds[i]
            push!(click_accum, i)
            push!(click_in_trains, res.pind_good[i])
            if !train_switch
                train_count += 1
                push!(train_start, res.pind_good[i])
                push!(train_start_ind, length(click_in_trains))
            end
            train_switch = true
        else
            if train_switch
                push!(train_end, res.pind_good[i])
            end
            train_switch = false
        end
    end
    if length(train_start) ==  length(train_end)+1
        push!(train_end, res.len_data)
    end

    res_new = (;pind_good_inS = res.pind_good_inS[click_accum], 
        pind_good = click_in_trains,
        ppeak = res.ppeak[click_accum],
        pind = res.pind_good,
        train_start, train_end,
        num_detection=train_count, num_click_in_trains=length(click_in_trains),
        train_start_ind, click_train_check_interval, click_train_minlen)

        
    if res_dir isa String
        !isdir(res_dir) || mkpath(res_dir)

        audacity_label(res_new.pind_good_inS, joinpath(res_dir, res.outfname *"__cps"*string((click_train_minlen+1)/click_train_check_interval)*  ".txt" |> basename))
        # label train
        if isempty(click_accum)
            @info "no click trains detected"
        else
            fs = round(Int, (res_new.pind_good[1] - 1) / res_new.pind_good_inS[1])
            audacity_label([train_start train_end] ./ fs, joinpath(res_dir, res.outfname *"__cps"*string((click_train_minlen+1)/click_train_check_interval)*  "_train-only.txt" |> basename))
            # audacity_label([train_start train_end] ./ fs, joinpath(res_dir, splitext(res.aufname)[1]*"_t"*string(res.threshold)*"_d"*string(res.dist) *"__cps"*string((click_train_minlen+1)/click_train_check_interval)*  "_train-only.txt" |> basename))
        end
    else
        @warn "detection not saved to file, no result folder(res_dir) given"
    end
    
    # @info "Number of clicktrains: " * string(sum(cltrain_start_inds ))
    @info "Number of clicktrains: " * string(train_count)
    @info "Number of impulse_in_train/impulse_detected: " * string(length(click_in_trains)) *"/"* string(length(res.pind_good))

    return res_new
end

#~ teager kaiser energy operator
function tkeo(data; type=Float64)
    data_n = zeros(size(data))
    data_n[1] = data[1]^2
    for i in 2:length(data)-1
        data_n[i] = data[i]^2 - (data[i-1]* type(data[i+1]))
    end
    data_n[end] = data[end]^2
    return data_n
end

using LinearAlgebra, StatsBase, Plots
"""
get the histogram of the inter pulse interval of the detected impulses
"""
function detect_impulsetrain2(aufname; 
    res_dir=nothing, ref_channel=1, dist_impulsive=80,
    bin_interval=100, time_interval = 0.1,
    fullplot=false, plot_everytimeintervalhistogram=false, display=x->x)
    # aufname = "/Users/abel/Documents/data_res/megafauana/clicks - all/click_HantuW_20210210_1.wav"
    # res_dir = "/Users/abel/Documents/data_res/megafauana/clicks_res"
    # ref_channel = 1
    # dist_impulsive=200

    data, fs, nbits, opt, timestamp = readAudio(aufname);
    res_impulse = detect_impulse((aufname,data,fs), res_dir; ref_channel=ref_channel, band_pass=[500 Inf], dist=dist_impulsive, threshold=nothing)#.01)
    times = res_impulse.pind_good
    timediff = zeros(Int,length(times),length(times)) 
    # using LinearAlgebra, StatsBase, Plots
    timediff = UpperTriangular(timediff)
    Threads.@threads for i in 1:length(times) 
        for j in i:length(times)
            timediff[i,j] = times[j] - times[i]
        end
    end
    # win=findall( x -> x>(283) && x<(284), times)
    # h = fit(Histogram, filter( x -> x > 0, timediff[win, win .+ i][:]); nbins=1000); println((i/fs, h.edges[1][argmax(h.weights)], maximum(h.weights)) ); plot(h; title=string(i))|>display

    max_t = ceil(res_impulse.pind_good_inS[end]) #300
    #~ good(time_interval = 0.1, nbins = 1000)
    # time_interval = 0.1; nbins = 1000;
    # time_interval = 0.05 #0.1   #0.2    #0.5
    # nbins = 200 #1000
    cum_res = zeros(max_t/time_interval |> Int, 3)

    numbers=[];
    ipi_direct=[];weights_direct=[];hs2=[];
    tt=[];ipi=[];weights=[];hs=[];
    for t in 0:time_interval:max_t-time_interval
        win=findall( x -> x>(t) && x<(t+time_interval), res_impulse.pind_good_inS)

        # h = fit(Histogram, filter( x -> x > 0, timediff[win, win][:]) ); 
        # h = fit(Histogram, filter( x -> x > 0, timediff[win, win][:]); nbins=1000 ); 
        # h = fit(Histogram, filter( x -> x > 0, timediff[win, win][:]), 0:100:10000 ); #best
        h = fit(Histogram, filter( x -> x > 0, timediff[win, win][:]), 0:bin_interval:10000 ); #best
        @debug h
        isempty(h.weights) && continue
        @debug ((t, h.edges[1][argmax(h.weights)], maximum(h.weights)) );
        push!(tt, t); 
        push!(ipi, h.edges[1][argmax(h.weights)]); 
        push!(weights, maximum(h.weights))
        push!(hs, h)

        h2 = fit(Histogram, filter( x -> x > 0, res_impulse.pind_good[win] |> diff), 0:bin_interval:10000 ); 
        push!(ipi_direct, h2.edges[1][argmax(h2.weights)]); 
        push!(weights_direct, maximum(h2.weights))
        push!(hs2, h2)

        push!(numbers, length(win))

        if plot_everytimeintervalhistogram
            p_ipi = plot(h; title=string(t))# |>display; 
            # savefig(joinpath(res_dir, basename(aufname)*"_click-segmentipi_all_" *string(t)* "s.png"))
            #### scatter!( h.edges[1][argmax(h.weights)], maximum(h.weights))|>display

            p_ipi_direct = plot(h2; title="direct-ipi "*string(t))# |>display; 
            # savefig(joinpath(res_dir, basename(aufname)*"_click-segmentipi_direct_" *string(t)* "s.png"))
            plot(p_ipi, p_ipi_direct; layout=@layout[a b], xlabel="ipi(sample)", ylable="weights(counts)");
            savefig(joinpath(res_dir, basename(aufname)*"_click-segmentipi_h" *string(t)* "s.png"))
            plot(p_ipi, p_ipi_direct; layout=@layout[a;b], xlabel="ipi(sample)", ylable="weights(counts)");
            savefig(joinpath(res_dir, basename(aufname)*"_click-segmentipi_v" *string(t)* "s.png"))
        end

    end
    arg_max = argmax(weights)
    plot(hs[arg_max]; title="max-weight_"* string(tt[arg_max]) *"s", xlabel="Inter Pulse Interval(samples)", ylabel="Counts") |> display
    if !isnothing(res_dir)
        savefig(joinpath(res_dir, basename(aufname)*"_click-max-weight.html"))
        savefig(joinpath(res_dir, basename(aufname)*"_click-max-weight.png"))
    end

    @info "max weighted ipi(weight): " * string(ipi[argmax(weights)]) *"("*string(maximum(weights))*") @ " * string(tt[argmax(weights)]) *"s"
    a= plot(tt, weights; label="counts", title=basename(aufname), ylabel="counts");
    b= plot(tt, ipi./res_impulse.fs; color=:red, label="ICI", ylabel="ICI(s)");
    c= plot(tt, ipi; color=:red, label="ICI_sample", ylabel="ICI(sample)");
    plot(a,b,c; layout=@layout[a;b;c], xlabel="Time(s)") |> display
    if !isnothing(res_dir)
        savefig(joinpath(res_dir, basename(aufname)*"_autoclick.html"))
        savefig(joinpath(res_dir, basename(aufname)*"_autoclick.png"))

        if fullplot
            gr()
            a= plot(tt, weights; label="weight", title=basename(aufname), ylabel="weight"); xlims!(0,size(data,1)/fs); xlabel!("Time(s)")
            b= plot(tt, ipi./res_impulse.fs; color=:red, label="ICI", ylabel="ICI(s)"); xlims!(0,size(data,1)/fs); xlabel!("Time(s)")
            c= plot(tt, ipi; color=:red, label="ICI_sample", ylabel="ICI(sample)"); xlims!(0,size(data,1)/fs); xlabel!("Time(s)")
            d= plot(signal(data,fs)); xlims!(0,size(data,1)/fs); xlims!(0,size(data,1)/fs); xlabel!("Time(s)")
            e= specgram(data; fs=fs,colorbar=nothing); xlabel!("Time(s)")
            f= heatmap( tt, hs[1].edges[1][2:end], hcat(map(x->x.weights, hs)...); colorbar=nothing, ylabel="ipi_all\n(sample)"); ylims!(0,3000); xlabel!("Time(s)")
            g= heatmap( tt, hs2[1].edges[1][2:end], hcat(map(x->x.weights, hs2)...); colorbar=nothing, ylabel="ipi\n(sample)"); ylims!(0,3000); xlabel!("Time(s)")
            counting = plot(tt, numbers; ylabel=("pulse\ncounts")); xlims!(0,size(data,1)/fs); xlabel!("Time(s)")
            histo_all = histogram(res_impulse.pind_good |> diff; bins=0:10:600, title=basename(aufname)*" IPI raw histogram", xlabel="Inter Pulse Interval(samples)", ylabel="Counts") #histogram(res_impulse.pind_good |> diff; bins=400, title=basename(aufname)*" IPI raw histogram", xlabel="Inter Pulse Interval(samples)", ylabel="Counts")
            weights_direct_plot= plot(tt, weights_direct; label="weight_direct", title=basename(aufname), ylabel="weight_direct"); xlims!(0,size(data,1)/fs); xlims!(0,size(data,1)/fs); xlabel!("Time(s)")
            ipi_direct_plot = c= plot(tt, ipi_direct; color=:red, label="ICI_sample", ylabel="ICI_direct\n(sample)"); xlims!(0,size(data,1)/fs); xlims!(0,size(data,1)/fs); ; ylims!(0,3000); xlabel!("Time(s)")

            plot(d,e,a,f,g,weights_direct_plot,counting,ipi_direct_plot,histo_all; layout=@layout[a;b;c;d;e;f;g;h;i;j], size=(1080,1080), legend=nothing) |> display
            # plot(a,d,e,b,f,g,counting,histo_all; layout=@layout[a;b;c;d;e;f;g;h;i], xlabel="Time(s)", size=(1080,920), legend=nothing) |> display

            # f= plot(tt, weights_direct; label="counts", ylabel="counts_direct");
            # plot(a,f,d,e; layout=@layout[a;b;c;d;e], xlabel="Time(s)", size=(1080,720)) |> display

            savefig(joinpath(res_dir, basename(aufname)*"_autoclick_full.png"))
            plotlyjs()
        end
    end

    p1=histogram(res_impulse.pind_good |> diff; bins=400, title=basename(aufname)*" IPI raw histogram", xlabel="Inter Pulse Interval(samples)", ylabel="Counts") 
    p2=histogram(res_impulse.pind_good |> diff; bins=0:10:3000, xlabel="Inter Pulse Interval(samples)", ylabel="Counts")
    p3=histogram(res_impulse.pind_good |> diff; bins=0:25:3000, xlabel="Inter Pulse Interval(samples)", ylabel="Counts") 
    plot(p1,p2,p3; layout=@layout[a b c]) |> display
    if !isnothing(res_dir)
        savefig(joinpath(res_dir, basename(aufname)*"_autoclick_raw-histogram.html"))
        savefig(joinpath(res_dir, basename(aufname)*"_autoclick_raw-histogram.png"))
    end

    return [ipi[argmax(weights)] ipi[argmax(weights)]/res_impulse.fs maximum(weights) tt[argmax(weights)] ], hs, res_impulse, timediff, times, ipi_direct, weights_direct
end

using DataFrames, CSV
function detect_impulsetrain2_folder(folname; res_dir=nothing, bin_interval=100, time_interval=0.1)
    # bin_interval=100; time_interval=0.1; #bin_interval=100; time_interval=0.1;
    folname = "/Users/abel/Documents/data_res/megafauana/clicks - all"
    res_dir = "/Users/abel/Documents/data_res/megafauana/clicks_res5_$bin_interval-bin_$time_interval-s__rawhisto-goodonly_new-direct_withcountNall3_shortheatmap_rawhistox2"
    rr = detect_impulsetrain2.(readdir(folname; join=true); res_dir = res_dir, fullplot=true, bin_interval=bin_interval, time_interval=time_interval)


    histo = map( x->x[2], rr)
    res_impulses = map( x->x[3], rr)
    timediff = map( x->x[4], rr)
    times = map( x->x[5], rr)
    rr = map( x->x[1], rr)

    timediff_single = diff.(times)
    # tds = fit.(Ref(Histogram), timediff_single; nbins=400)
    # timediff_single_stats = map( td -> [td.edges[1][ argmax(td.weights) ] maximum(td.weights)] , tds)
    tds = fit.(Ref(Histogram), timediff_single, Ref(0:5:3000))
    plot2(x,y; kwargs...) = plot(x; title=y, kwargs...)
    plot2.(tds, fnames; xlims=(0,2000)) .|> display
    timediff_single_stats = map( td -> [td.edges[1][ argmax(td.weights) ] maximum(td.weights)] , tds)

    histogram(vcat(timediff_single_stats...)[:,1]; bins=0:10:600, xlabel="Inter Pulse Interval(samples)", ylabel="Counts")
    savefig(joinpath(res_dir, "raw-histogram_IPI-mode.html"))
    savefig(joinpath(res_dir, "raw-histogram_IPI-mode.png"))

    fnames = readdir(folname);
    # jldsave(joinpath(res_dir,"result.jld2"); histo,timediff,times,rr,fnames, timediff_single)

    # using DelimitedFiles
    # writedlm( joinpath(res_dir,"table_clicks.csv"),  vcat(rr...), ',')

    # rr2 = vcat(rr...);
    # writedlm( joinpath(res_dir,"table_clicks.csv"),
    #     vcat(["filename" "ICI_sample" "ICI_seconds" "weight" "time_s"],
    #         hcat(readdir(folname), rr2)),
    #     ','
    # )

    # df = DataFrame(hcat(readdir(folname), rr2), ["filename", "ICI_sample", "ICI_seconds", "weight","time_s"])

    # hcat(rr2, vcat(timediff_single_stats...))
    # using DataFrames
    df = DataFrame(hcat(readdir(folname), vcat(rr...), vcat(timediff_single_stats...)),
        ["filename", "ICI_sample", "ICI_seconds", "weight","time_s","ici_direct_mode","ici_direct_weights"])
    # using CSV
    CSV.write(joinpath(res_dir,"table_clicks.csv"), df)

    jldsave(joinpath(res_dir,"result.jld2"); histo,timediff,times,rr,fnames, timediff_single, df)
end



# for i in 1:length(d["times"]); plot(d["times"][i][1:end-1]./fs, d["timediff_single"][i]; title=d["fnames"][i]) |> display; end




    #     cum_res[t/time_interval |> Int,:] = [t h.edges[1][argmax(h.weights)] maximum(h.weights)]

    # for t in 0:time_interval:max_t
    #     win=findall( x -> x>(t) && x<(t+1), res_impulse.pind_good_inS)
    #     # isempty(win) && continue
    #     length(win) < 2 && continue
    #     h = fit(Histogram, filter( x -> x > 0, timediff[win, win .+ t][:]); nbins=1000 ); println((t, h.edges[1][argmax(h.weights)], maximum(h.weights)) ); plot(h; title=string(t))|>display
    #     cum_res[t/time_interval |> Int,:] = [t h.edges[1][argmax(h.weights)] maximum(h.weights)]
    #     println(t)
    # end

    # scatter(cum_res[:,1], cum_res[:,2]./fs)

    # for i in 0:fs:fs*300
    #     h = fit(Histogram, filter( x -> x > 0, timediff[win, win .+ i][:]), ); println((i/fs, h.edges[1][argmax(h.weights)], maximum(h.weights)) ); plot(h; title=string(i))|>display
    # end
# end

# function detect_impulsetrain2(res_impulse)

 