using Pkg
Pkg.activate(".")

using ProgressMeter
include("detector.jl")
include("plotting.jl"); plotlyjs()
include("tabulate_data.jl")

set_device__hydromoth()


data_folders = 
[
# "/media/spin/anas2/data/marecet/Acoustics-Abel/Hydromoth - PAM/Feb 2023 - LGK/H1",
"/media/spin/anas2/data/marecet/Acoustics-Abel/Hydromoth - PAM/Oct 2021 - LGK/H5 - Teluk Apau/23102021 - 26102021"
]
ftype = ".WAV"
result_directory = "/media/spin/anas2/data_res/dolphin/marecet/langkawi_hydromoth/20261020_1900_nbhf"

# simple check of file name to configure NAME2DT_DATE__TIME_SKIP_FRONT automatically
files = readdir(data_folders[1]; join=true) |> filter(endswith(ftype))
files[1]
for ind=-1:100
    try
        DEFAULT_fname2timestamp_func(files[1]; skip_front=FNAME2DT_DATE__TIME_SKIP_FRONT)
    catch e
        @warn "Failed to get timestamp from filename: $e; FNAME2DT_DATE__TIME_SKIP_FRONT=$FNAME2DT_DATE__TIME_SKIP_FRONT"
        FNAME2DT_DATE__TIME_SKIP_FRONT = ind
    end
end
if FNAME2DT_DATE__TIME_SKIP_FRONT == 100
    @error "Failed to determine correct FNAME2DT_DATE__TIME_SKIP_FRONT"
else
    @info "Final FNAME2DT_DATE__TIME_SKIP_FRONT: $FNAME2DT_DATE__TIME_SKIP_FRONT\n"
end



# # num_proc_per_run = 2
# # proc_free = true(num_proc_per_run)

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
            # detect_impulseNtonal(file, res_dir; processed_skip_flag=true, processed_skip_strict=false)
        catch e
            @error "Error processing file $file: $e"
        end
        GC.gc() # run garbage collector to free memory
    end
end

sleep(10)
print("Exiting... wait 10 seconds")
exit(0)

#~ multiple processes
# num_proc_per_run = 2
# for aufol in data_folders
#     @info "Processing folder: $aufol"
#     files = readdir(aufol; join=true) |> filter(endswith(ftype))
#     files = vcat(files, readdir(aufol; join=true) |> filter(endswith("wav")))
#     res_dir = joinpath(result_directory, Dates.format(DEFAULT_fname2timestamp_func(files[1]), "yyyymmdd_HHMMSS"))
#     isnothing(res_dir) || mkpath(res_dir)
#     isfile(joinpath(res_dir, "index.html")) || cp("web/index.html", joinpath(res_dir, "index.html"))

#     # Load all files into the channel, then close so workers stop when empty
#     ch = Channel{String}(length(files))
#     foreach(f -> put!(ch, f), files)
#     close(ch)

#     workers = map(1:num_proc_per_run) do _  # exactly 2 worker tasks
#         Threads.@spawn begin
#             for file in ch           # each worker pulls next available file
#                 fname, ext = splitext(basename(file))
#                 @info "Processing file: $(dirname(file))/\033[32m$fname\033[0m$ext"
#                 try
#                     # pause(10) # for testing, to simulate long processing time
#                     detect_impulseNtonal(file, res_dir;
#                         detect_impulse        = detect_impulseNarrowBand,
#                         processed_skip_flag   = true,
#                         processed_skip_strict = false)
#                 catch e
#                     @error "Error processing file $file: $e"
#                 end
#                 GC.gc()
#             end
#         end
#     end

#     foreach(wait, workers)   # block until both workers have drained the channel
#     df = CSV.read(joinpath(res_dir, "counts.csv"), DataFrame)
#     df_new = sort(df, :datetime)
#     # check if there's any repeated rows, if so, remove them
#     # df_new!=unique(df_new) && @warn "Removing repeated rows in counts.csv"# && CSV.write(joinpath(res_dir, "counts.csv"), unique(df_new))
#     df!=df_new && CSV.write(joinpath(res_dir, "counts.csv"), df_new)
# end