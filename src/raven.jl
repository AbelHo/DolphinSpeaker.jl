export raven_label

header = "Selection\tView\tChannel\tBegin Time (s)\tEnd Time (s)\tNotes\n"

function raven_label(event_time, io=stdout; delim="\t", prefix="", channel=1, view_type="Spectrogram 1", label_list=nothing)
    # write raven label file
    # labelfile = res_dir*fname[end-22:end]*"_label.txt"
    # f=open(io, "w")
    write(io, header)
    if size(event_time,2)==1
        for i in 1:length(event_time)
            label = isnothing(label_list) ? prefix * string(i) : label_list[i]
            
            write(io, string(i) *
            delim * view_type *
            delim * string(channel) *
            delim*string(event_time[i]) *
            delim*string(event_time[i]) *
            delim* label *"\n")
        end
    elseif size(event_time,2)==2
        for i in 1:size(event_time,1)
            label = isnothing(label_list) ? prefix * string(i) : label_list[i]

            write(io, string(i) *
            delim * view_type *
            delim * string(channel) *
            delim*string(event_time[i,1]) *
            delim*string(event_time[i,2]) *
            delim* label *"\n")
        end
    end
    # close(io)
    # labelfile
end

function raven_label(event_time, fname::String; kwargs...)
    # labelfile = res_dir*fname[end-22:end]*"_label.txt"
    f=open(fname, "w")
    raven_label(event_time, f; kwargs...)
    close(f)
end


# using Dates
# using Printf

# function raven_label2snippets(input_wav::String, labels_file::String, output_dir::String)
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
#                 start, last, label audacitypeg command
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

# raven_label2snippets(input_wav, labels_file, output_dir)


using CSV
using DataFrames

"""
    convert_raven_to_audacity(raven_file::String, audacity_file::String)

Converts a Raven selection file into an Audacity label file format.

# Arguments
- `raven_file`: Path to the Raven selection file (tab-delimited).
- `audacity_file`: Path to save the Audacity label file.

# Notes
The Raven file should have columns `Begin Time`, `End Time`, and optionally `Label`.
"""
function convert_raven_to_audacity(raven_file::String, new_folder::String=raven_file|>dirname; suffix="_audacity")
    if isdir(raven_file)
        for rfile in readdir(raven_file; join=true)
            @debug "Processing file: $rfile"
            # Convert each Raven file to Audacity format
            return convert_raven_to_audacity_one(rfile, new_folder; suffix=suffix)
        end
    end
    # Read the Raven selection file into a DataFrame
    df = CSV.read(raven_file, DataFrame, delim='\t')
    if "Begin File" in names(df)
        sets = groupby(df, Symbol("Begin File"))
    else
        sets = [df]
    end
    # if length(sets) > 1
    df.Duration = df[:,"End Time (s)"] - df[:,"Begin Time (s)"]
    "Offset_endtime" in names(df) || (df.Offset_endtime = df[:,"File Offset (s)"] + df.Duration)
    # end

    for (i, set) in enumerate(sets)
        fname = "Begin File" in names(df) ? set[begin,:"Begin File"] : basename(raven_file)
        # Construct the output file name
        # if audacity_file == ""
        #     audacity_file = replace(raven_file, r"\.txt$" => string(suffix, ".txt"))
        # end
        # audacity_file = replace(fname, r"\.txt$" => string("_", i, suffix, ".txt"))
        # audacity_file = splitext(raven_file)
        audacity_file = joinpath(new_folder, "$fname$suffix.txt")
        @debug fname, joinpath(new_folder,audacity_file)

        # Convert the current set to Audacity format
        convert_raven_to_audacity_one(set, joinpath(new_folder,audacity_file))
    end
end

function convert_raven_to_audacity_one(df::Union{DataFrame, SubDataFrame}, audacity_file::String)
    # Read the Raven selection file into a DataFrame
    # df = CSV.read(raven_file, DataFrame, delim='\t')

    # Open the Audacity label file for writing
    open(audacity_file, "w") do io
        i = 1
        for row in eachrow(df)
            begin_time = "File Offset (s)" in names(df) ? row["File Offset (s)"] : row["Begin Time (s)"]
            end_time = "Offset_endtime" in names(df) ? row["Offset_endtime"] : row["End Time (s)"]
            label = string(i) *" " * string(haskey(row, "Notes") ? row["Notes"] : "") * string(haskey(row, "Note") ? row["Note"] : "")

            # Write the Audacity label format: start_time, end_time, label
            println(io, "$begin_time\t$end_time\t$label")

            if "Low Freq (Hz)" in names(df) && "High Freq (Hz)" in names(df)
                low_freq = row["Low Freq (Hz)"]
                high_freq = row["High Freq (Hz)"]
                if !(ismissing(low_freq) || ismissing(high_freq) || isempty(low_freq) || isempty(high_freq))
                    println(io, "\\\t$low_freq\t$high_freq")
                end
            end
            i += 1
        end
    end
end
function convert_raven_to_audacity_one(raven_file::String, audacity_file::String; delim='\t')
    # Ensure required columns exist
    # if !all(["Begin Time", "End Time"] .∈ names(df))
    #     error("Raven file must contain 'Begin Time' and 'End Time' columns.")
    # end

    # Read the Raven selection file into a DataFrame
    df = CSV.read(raven_file, DataFrame; delim=delim)
    convert_raven_to_audacity_one(df, audacity_file)
end

# Example usage:
# convert_raven_to_audacity("raven_selection.txt", "audacity_labels.txt")



# function convert_raven_to_audacity(raven_file::String, audacity_file::String)
#     # Read the Raven selection file into a DataFrame
#     df = CSV.read(raven_file, DataFrame, delim='\t')

#     # Ensure required columns exist
#     required_columns = ["Begin Time (s)", "End Time (s)"]
#     if !all(col -> col ∈ names(df), required_columns)
#         error("Raven file must contain 'Begin Time' and 'End Time' columns.")
#     end

#     # Open the Audacity label file for writing
#     open(audacity_file, "w") do io
#         for row in eachrow(df)
#             begin_time = row["Begin Time"]
#             end_time = row["End Time"]
#             label = haskey(row, "Label") ? row["Label"] : ""

#             # Write the Audacity label format: start_time, end_time, label
#             println(io, "$begin_time\t$end_time\t$label")
#         end
#     end
# end