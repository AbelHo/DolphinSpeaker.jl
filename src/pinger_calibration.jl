include("PAM.jl")

using Peaks, SignalAnalysis, JLD2

res_dir= homedir()*"/Desktop/results/"#"/Desktop/data/aspod/results/"
if !isdir(res_dir)
    mkpath(res_dir)
end

# fs=400000
# channels=4
# ref_channel=1

"""
    findPings(fname::String; save_fn=res_dir*fname[end-22:end]*"_label.jld2", ref_channel=1::Int, dist=nothing, channels=5, fs=400000)
    findPings(fname::AbstractArray; ref_channel=1::Int, dist=nothing)
Find peaks in the the signal by specifying a reference channel



# Examples
## For 1 second interval pinger & chatterbox XL with 5 channels 
```julia-repl
julia> peaks, peak_times, peak_vals = findPings("/somepath/filename.bin", save_fn=nothing, dist=Int(400000*0.9), channels=5)
1
```
"""
function findPings(fname::String; save_fn=res_dir*fname[end-22:end]*"_label.jld2", ref_channel=1::Int, dist=nothing, channels=5, fs=400000)
    if !isnothing(save_fn) && isfile(save_fn)
        peaks=load(save_fn)
        peaks = peaks["peaks"]
    else   
        data=loadDataBin(fname, channels=channels)
        peaks = findPings(data, ref_channel=ref_channel, dist=dist)
        if !isnothing(save_fn)
            save(save_fn, "peaks", peaks, "peak_times", peaks[1]./fs, "peak_vals", peaks[2]);
        end
    end
    println(string(length(peaks[1]))*" pings")
    return peaks, peaks[1]./fs, peaks[2]
end


function findPings(data::AbstractArray; ref_channel=1::Int, dist=nothing)
    if isnothing(dist)
        peaks = findmaxima(data[:,ref_channel])
    else 
        peaks = findmaxima(data[:,ref_channel], dist)
    end
    return peaks
end



function audacity_label(event_time, io=stdout)
    # write audacity label file
    # labelfile = res_dir*fname[end-22:end]*"_label.txt"
    # f=open(io, "w")
    for i in 1:length(event_time)
        write(io, string(event_time[i]) *"\t"*string(event_time[i]) *"\t"* string(i) *"\n")
    end
    # close(io)
    # labelfile
end

function audacity_label(event_time, fname::String)
    # labelfile = res_dir*fname[end-22:end]*"_label.txt"
    f=open(fname, "w")
    audacity_label(event_time, f)
    close(f)
end

function findNauda(fname::String, outFname=res_dir*fname[end-22:end]*"_label.txt"::String, fs=400000 ) #
    peaks,_,_ = findPings(fname; dist=Int(fs*.9))
    if !isempty(outFname)
        audacity_label(peaks[1]./fs, outFname)
        println(outFname)
    end
    return peaks
end
