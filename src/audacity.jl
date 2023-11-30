export audacity_label

function audacity_label(event_time, io=stdout; prefix="")
    # write audacity label file
    # labelfile = res_dir*fname[end-22:end]*"_label.txt"
    # f=open(io, "w")
    if size(event_time,2)==1
        for i in 1:length(event_time)
            write(io, string(event_time[i]) *
            "\t"*string(event_time[i]) *
            "\t"* prefix * string(i) *"\n")
        end
    elseif size(event_time,2)==2
        for i in 1:size(event_time,1)
            write(io, string(event_time[i,1]) *
            "\t"*string(event_time[i,2]) *
            "\t"* prefix * string(i) *"\n")
        end
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