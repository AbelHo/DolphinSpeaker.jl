include("detector.jl")
include("plotting.jl"); plotlyjs()
include("tabulate_data.jl")

set_device__soundtrap()

#~ tabulate raw data
data_folders = 
["/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard4",
"/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3",
"/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard2",
"/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024/sdcard3"
]
data_folders = 
["/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3/wav"]
result_directory = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_$(Dates.format(now(),"yyyymmdd"))/raw/$(Dates.format(now(),"HHMMSS"))"

ftype = ".flac"
# tabulate data
# summary_filepath = "/media/spin/anas2/data/marecet/datai/summary.csv"
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024"; output=summary_filepath, filetype="wav", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_20032024_30052024"; output=summary_filepath, filetype="flac", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
# df = CSV.read(summary_filepath, DataFrame)
# df = sort(df, :datetime)
# CSV.write(summary_filepath, df)


# analyse
for aufol in data_folders
    @info "Processing folder: $aufol"
    files = readdir(aufol; join=true) |> filter(endswith(ftype))
    res_dir = joinpath(result_directory, Dates.format(DEFAULT_fname2timestamp_func(files[1]), "yyyymmdd_HHMMSS"))
    isnothing(res_dir) || mkpath(res_dir); cp("web/index.html", joinpath(res_dir, "index.html"))
    for file in files
        fname, ext = splitext(basename(file))
        @info "Processing file: $(dirname(file))/\033[32m$fname\033[0m$ext"
        try
            detect_impulseNtonal(file, res_dir; detect_impulse=detect_impulseNarrowBand)
        catch e
            @error "Error processing file $file: $e"
        end
        GC.gc() # run garbage collector to free memory
    end
end




pp = detectionsfiles2plot2.(
    joinpath.( filter(isdir,readdir(res_dir; join=true)), "counts.csv");
    res_dir=res_dir,
    plottype=PlotlyJS.bar, add_daynight=true);



# detectionsfiles2plot("/media/spin/anas2/data_res/dolphin/marecet/datai/res_20250724/deployment_18012024_18042024/counts.csv"; res_dir="/media/spin/anas2/data_res/dolphin/marecet/datai/res_20250724/deployment_18012024_18042024", plottype=PlotlyJS.bar)
