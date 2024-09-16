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


# using Dates
# using Printf

# function audacity_label2snippets(input_wav::String, labels_file::String, output_dir::String)
#     # Create the output directory if it doesn't exist
#     if !isdir(output_dir)
#         mkpath(output_dir)
#     end

#     # Open the labels file and process each line
#     open(labels_file, "r") do file
#         for line in eachline(file)
#             # Split the line by tab character
#             parts = split(line, '\t')
#             if length(parts) == 3
#                 start, last, label = parts
#                 output_file = joinpath(output_dir, "$label.wav")
                
#                 # Construct the ffmpeg command
#                 cmd = `ffmpeg -i $input_wav -ss $start -to $end -c copy $output_file`
                
#                 # Run the ffmpeg command
#                 run(cmd)
#             end
#         end
#     end
# end

# # Example usage
# input_wav = "input.wav"
# labels_file = "labels.txt"
# output_dir = "output_dir"

# audacity_label2snippets(input_wav, labels_file, output_dir)