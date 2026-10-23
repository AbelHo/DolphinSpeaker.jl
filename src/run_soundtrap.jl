using ProgressMeter
include("detector.jl")
include("plotting.jl"); plotlyjs()
include("tabulate_data.jl")

set_device__soundtrap()

#~ tabulate raw data
data_folders = 
["/media/spin/anas2/data/marecet/examples",
]
# data_folders = 
# ["/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3/wav"]
result_directory = "/media/spin/anas2/data_res/dolphin/marecet/temp/soundtrap/test_only/res_$(Dates.format(now(),"yyyymmdd"))/raw/$(Dates.format(now(),"HHMMSS"))"
# result_directory = "/media/spin/anas2/data_res/dolphin/marecet/datai/res_20260619/raw/012437"

ftype = ".wav"
# tabulate data
# summary_filepath = "/media/spin/anas2/data/marecet/datai/summary.csv"
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024"; output=summary_filepath, filetype="wav", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024"; output=summary_filepath, filetype="flac", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
# df = CSV.read(summary_filepath, DataFrame)
# df = sort(df, :datetime)
# CSV.write(summary_filepath, df)

# debug only one file
# ENV["JULIA_DEBUG"] = "detector"

#~ temp delete
data_folders = 
[
# "/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard4",
"/media/spin/anas2/data/marecet/datai/deployment_3_27022025_090702025",
]
result_directory = "/media/spin/anas2/data_res/dolphin/marecet/temp/soundtrap/test_only"

# num_proc_per_run = 2
# proc_free = true(num_proc_per_run)

# analyse
for aufol in data_folders
    @info "Processing folder: $aufol"
    files = readdir(aufol; join=true) |> filter(endswith(ftype))
    files = vcat(files, readdir(aufol; join=true) |> filter(endswith("wav")))
    
    res_dir = joinpath(result_directory, "no_date")
    try
        res_dir = joinpath(result_directory, Dates.format(DEFAULT_fname2timestamp_func(files[1]), "yyyymmdd_HHMMSS"))
    catch e
        # res_dir = joinpath(result_directory, "no_date")
        # continue
    end
    isnothing(res_dir) || mkpath(res_dir); isfile(joinpath(res_dir, "index.html")) || cp("web/index.html", joinpath(res_dir, "index.html"))
    @info files
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