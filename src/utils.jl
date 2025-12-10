# using SignalAnalysis, SignalAnalysis.Units
# using WAV
# using FFMPEG
# using Dates
# import TimeZones
# include("audio.jl")
import FileIO, FileIO.load

showall(x) = show(stdout, "text/plain", x)

function skiphiddenfiles(list)
    filter(!startswith(r"^[.@]") ∘  basename, list)
end

function savejld(savefname; kwargs...)
    # savefname = splitext(basename(aufname))[1] *"_t"*string(impulsive_autothreshold_median_ratio)*"_d"*string(dist_impulsive)*".jld2"
    if !Sys.islinux()
        jldsave(savefname; kwargs...)
        # jldsave(joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".jld2"); pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect,  vidfname, aufname, res_dir, detector_set, pixel_related)
    #@error(err)
    #@error("Cant save jld2 file, saving locally and copying instead")
    #rm(joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".jld2"))
    else
        jldsave(joinpath("", savefname|>basename); kwargs...)
        mv(joinpath("", savefname|>basename), savefname; force=true)
    end
    @info "Saved to: " * savefname
end

function FileIO.load(filenames::Array{String,1}; kwargs...)
    merge( load.(filenames; kwargs...)...)
end

"""
    process_files(folname; func=(a,b)->x, arg=nothing, no_overwrite_func=nothing)

Process files within a directory tree starting from `folname`.

# Arguments
- `folname::String`: The root directory from which to start processing files.
- `func`: A function to be applied to each file. It should accept two arguments: the source path and the destination path. It defaults to a function that does nothing.
- `arg::Union{Nothing, String}`: An optional argument. If it is a string, it is treated as a directory path, and `mkpath` is called to ensure this directory exists. This path is used as the base for the destination path in `func`.
- `no_overwrite_func`: An optional function that determines whether to overwrite existing files. It should accept the same arguments as `func`.

# Behavior
- Iterates over every file in the directory tree rooted at `folname`.
- For each file, it constructs a source path (from `folname`) and a destination path (from `arg`).
- It then attempts to apply `func` to these paths. If `func` fails, it tries to call `func` without the `no_overwrite_func` argument. If this also fails, it logs the error.

# Example
```julia
process_files("path/to/source", func=(src, dst) -> copy(src, dst), arg="path/to/destination")
"""
function process_files(folname; func=(a,b)->x, arg=nothing, no_overwrite_func=nothing)
    arg isa String && mkpath(arg)
    for (root, dirs, files) in walkdir(folname)
        # println("Directories in $root")
        # for dir in dirs
        #     println(joinpath(root, dir)) # path to directories
        #     try
        #         func(joinpath(root, dir), joinpath(arg, dir); no_overwrite_func=no_overwrite_func)
        #     catch err
        #         try 
        #             func(joinpath(root, dir), joinpath(arg, dir))
        #         catch err
        #             @error exception=(err, catch_backtrace())
        #             @error (joinpath(root, dir), joinpath(arg, dir))
        #         end
        #     end
        # end
        # println("Files in $root")
        for file in files
            println(joinpath(root, file)) # path to files
            try
                # @info joinpath(root, file), joinpath(arg, file)
                func(joinpath(root, file), joinpath(arg, file); no_overwrite_func=no_overwrite_func)
            catch err
                try 
                    func(joinpath(root, file), joinpath(arg, file))
                catch err
                    @error exception=(err, catch_backtrace())
                    @error (joinpath(root, file), joinpath(arg, file))
                end
            end
        end
    end
end

function replace_suffix(src::AbstractString, suffix, replacement=""; preview=true, kwargs...)
    if isdir(src) 
        replace_suffix.(readdir(src; join=true), Ref(suffix), Ref(replacement); preview=preview, kwargs...)
        return
    end

    splitted = splitext(src)
    if endswith(splitted[1], suffix)
        @info splitted, suffix
        indices_of_suffix = findlast(suffix, splitted[1])
        isempty(indices_of_suffix) && return
        
        dst = splitted[1][first:[1]-1] * replacement * splitted[2]
        if preview
            println(basename(src) *"\t"* basename(dst))
            return
        end
        mv(src, dst; kwargs...)
    end
end

# Function to search for files with specific string and suffix, then run a function on threshold_impulsive
function func_on(a, b; needle="", filetype="", prefix="", func=(x,y)->nothing, verbose=false)
    fname = basename(a)
    verbose && print("\t",a,b, "\n")
    if occursin(needle, fname) && endswith(fname, filetype) && startswith(fname, prefix)
        func(a, joinpath(b, fname))
    end
end

"""
# Transfer some files to a new folder with the same folder structure as the original folder.
# example:
infol = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648"
outfol = "/home/spin/Documents/data_res/marecet/temp/Datai_20250724/detection_files2"
# func_new(x, y) = func_on(x, y; needle="cps40.0", filetype=".txt", func=(a,b)->print("\t",a,"  ",b,"\n"), verbose=false)
# transfer raven files
func_new(x, y) = func_on(x, y; prefix="raven_", needle="cps40.0", filetype=".txt", func=cp, verbose=false)
process_files_hierarchy_with_func(infol, outfol; func=func_new)#verbose=true)
# transfer counts.csv files
func_new(x, y) = func_on(x, y; prefix="counts", needle="", filetype=".csv", func=cp, verbose=false)
process_files_hierarchy_with_func(infol, outfol; func=func_new)#verbose=true

count_files = glob("*/counts.csv", infol) |> sort
df = CSV.read.(count_files, DataFrame)
select!(df[2], Not(:filepath)); df = vcat(df...)
CSV.write(joinpath(outfol, "counts_all.csv"), df)
# df = vcat(CSV.read.(count_files, DataFrame)...)
"""
function process_files_hierarchy_with_func(folname, outfolder, args...; verbose=false, func=(x,y)->nothing, flag_skiphiddenfiles=true, kwargs...)
    flist = readdir(folname; join=true) 
    flag_skiphiddenfiles && (flisth = flist |> skiphiddenfiles)
    for file in flist
        if isdir(file)
            @info "directory: $file"
            mkpath(joinpath(outfolder, basename(file)))
            process_files_hierarchy_with_func( file, joinpath(outfolder, basename(file)), args...; 
                func=func, verbose=verbose, flag_skiphiddenfiles=flag_skiphiddenfiles, kwargs...)
        else
            verbose && print("file: $file, outfolder: $outfolder\n")
            func(file, outfolder, args...; kwargs...)
        end
    end
end


# using Glob
# """
    # make_clip_index_html(folder)# flag_sort_number=false)#; outname="index.html", prefix="clip_")
# """
function make_clip_index_html(folder::AbstractString; outname="index2.html", prefix="clip_", flag_sort_number=true, output_type="html")
    # Extract number from filename for sorting
    function clipnum(f)
        m = match(Regex("^" * prefix * "(\\d+)\\.$output_type\$"), basename(f))
        isnothing(m) ? typemax(Int) : parse(Int, m.captures[1])
    end

    files = filter(f -> occursin(Regex("^" * prefix * "\\d+\\.$output_type\$"), basename(f)), readdir(folder; join=true))
    if flag_sort_number 
        files = sort(files, by=clipnum)
    else
        files = sort(files)
    end
    outpath = joinpath(dirname(folder), outname)

    open(outpath, "w") do io
        println(io, "<!DOCTYPE html>\n<html>\n<head><title>Clips Index</title></head>\n<body>")
        for f in files
            fname = basename(f)
            relpath = joinpath(basename(folder), fname)
            println(io, """<h3>$fname</h3>""")
            if output_type in ["png", "jpg", "jpeg", "gif", "webp", "svg", "svg", "bmp", "tiff"]
                println(io, """<img src="$relpath" alt="$fname" style="max-width:100%; height:auto;">""")
            else
                println(io, """<iframe src="$relpath" width="100%" height="400" style="border:none;margin-bottom:1em;"></iframe>""")
            end
        end
        println(io, "</body>\n</html>")
    end
    return outpath
end



"""
    run_func_fileauto(dname, outfolder; sensor_names=["acoustic", "topview", "uw1"], sensor_filetypes=[".ogg", ".mkv", ".mkv"], func=x->x, prefix_filter="", kwargs...)

Process files in a directory `dname` based on specified `sensor_names` and `sensor_filetypes`, applying a function `func` to each file group and outputting the results to `outfolder`.

# Arguments
- `dname::String`: The directory name where the source files are located.
- `outfolder::String`: The output directory where the processed files will be saved.
- `sensor_names::Array{String}`: An array of sensor names, used to identify groups of files. Defaults to `["acoustic", "topview", "uw1"]`.
- `sensor_filetypes::Array{String}`: An array of file extensions corresponding to each sensor name. Defaults to `[".ogg", ".mkv", ".mkv"]`.
- `func`: A function to be applied to each group of files. It defaults to an identity function (`x->x`). The function should accept file paths as arguments and any number of keyword arguments (`kwargs`).
- `prefix_filter::String`: A string filter to apply to filenames, selecting only those that start with the specified prefix. Defaults to an empty string, which selects all files.
- `kwargs...`: Additional keyword arguments to be passed to `func`.

# Behavior
- Creates the `outfolder` if it does not exist.
- Iterates over files in the first sensor's directory (`sensor_names[1]`), filtering by `prefix_filter`.
- Constructs a group of file paths for each sensor based on the common prefix and specified file types.
- Applies `func` to each group of file paths, passing `outfolder` and any `kwargs` as arguments.
- Logs an error if `func` fails to process a group of files.

# Example
template function:
    ```
    stack_audio_videos(aufname, v1, v2, res_dir; kwargs...)
    ```
```julia
run_func_fileauto("data", "processed", func=(files..., outfolder; kwargs...) -> println("Processing: ", files, " into ", outfolder), prefix_filter="2021_")
"""
function run_func_fileauto(dname, outfolder,
    sensor_names = ["acoustic", "topview", "uw1"], sensor_filetypes = [".ogg", ".mkv", ".mkv"];
    func=x->x, prefix_filter="", kwargs...)

    mkpath(outfolder)
    # dname = dirname(firstfolder)
    for fname in readdir(joinpath(dname, sensor_names[1]))|>skiphiddenfiles |> x->filter(startswith(prefix_filter), x) #FIXME will fail in year 2100 onwards
        fname_split = splitext(fname)[1]
        fname_split = fname_split[1:findlast('_', fname_split)-1]

        # thisfiletype = fname_split[2]
        # fname_split = fname_split[1]
        # if thisfiletype != filetype
        #     continue
        # end
        @info fname

        try
            # @info joinpath.(Ref(dname),sensor_names,fname_split .* "_" .* sensor_names .* sensor_filetypes)
            func(joinpath.(Ref(dname),sensor_names,fname_split .* "_" .* sensor_names .* sensor_filetypes)..., outfolder; kwargs...)
            # func(joinpath(dname,"cam_topview",fname_split*"_topview.mkv"), joinpath(dname,"cam_uw1",fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
        catch err
            @error "ERROR run_func_fileauto"
            @error exception=(err, catch_backtrace())
            # func(joinpath(aufolder,fname_split*"_topview.mkv"), joinpath(aufolder,fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
        end
    end
end

"""
    print out the argument of a this function
"""
print_args(args...) = println(args)

"""
# Example
process_files("/Volumes/data/Concretecho/data/temp/try8"; func=(a,b)->template_func(a,b; template=["Ambient_", ".png"], func=cp), arg="/Volumes/data/Concretecho/data/temp/Ambient_pic") 
"""
function template_func(fp, res_file=""; template=[""], func=x->x)
    fn = basename(fp)
    if startswith(fn, template[begin]) && length(template)>1 && endswith(fn, template[end])
        @info (fp, res_file)
        return func(fp, res_file)
    end
end

# macro threaded(expr)
#     if expr.head == :comprehension
#         @info expr[1]
#         loop_vars = expr.args[2]
#         body = expr.args[1]
#         return :(
#             let
#                 items = collect($(esc(loop_vars)))
#                 Threads.@threads for i in eachindex(items)
#                     items[i] = $(esc(body))
#                 end
#                 items
#             end
#         )
#     else
#         throw(ArgumentError("The @parallel_comprehension macro expects a comprehension expression."))
#     end
# end
# Convert Dict to NamedTuple
dict2namedtuple(d) = NamedTuple(Symbol(k) => v for (k, v) in d)

function find_filse_with_string(folname, str)
    # folname = "/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3/wav"
    files = readdir(folname; join=true) |> filter(contains(str))
    if isempty(files)
        @warn "No files found containing $str in $folname"
    end
    return files
end
function find_files_with_suffix(folname, suffix)
    # folname = "/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024/SDcard3/wav"
    files = readdir(folname; join=true) |> filter(endswith(suffix))
    if isempty(files)
        @warn "No files found with suffix $suffix in $folname"
    end
    return files
end

readdirjoin(x) = readdir(x; join=true)

function concat_files(files::Vector{String}, out::String)
    open(out, "w") do o
        for f in files
            write(o, read(f))   # read returns Vector{UInt8} for binary
        end
    end
end

"""
function: returns a Vector{String} of all files (full paths) under `dir`
"""
function readdir_all(dir::AbstractString; follow_symlinks::Bool=true, kwargs...)
    paths = String[]
    for (root, _, files) in walkdir(dir, follow_symlinks=follow_symlinks)
        for f in files
            push!(paths, joinpath(root, f))
        end
    end
    return paths
end

"""
function: returns a lazy iterator of full file paths (use collect(...) to get a Vector)
"""
readdir_all_iter(dir::AbstractString; follow_symlinks::Bool=true, kwargs...) = 
    (joinpath(root, f) for (root, _, files) in walkdir(dir, follow_symlinks=follow_symlinks) for f in files)

@info "end utils.jl"

# function findTrigger(data, fs; threshold_percentMAX=0.75, plot_window_inS=nothing, ref_channel=size(data,2)) # plot_window_inS=[-.1 .1]
#     # @debug ref_channel
#     # @debug abs.(data[:,ref_channel])
#     trigger = findfirst(abs.(data[:,ref_channel]) .> maximum(data[:,ref_channel])*threshold_percentMAX )
#     # @debug trigger
#     isnothing(plot_window_inS) && return trigger
#     # @debug fs
#     plot(data[:, ref_channel])
#     @debug "after plot"
#     scatter!([trigger], [data[trigger,ref_channel]])
#     xlims!(trigger+fs*plot_window_inS[1], trigger+fs*plot_window_inS[2])
#     # return p
# end

# function findVidAudioBlip(fname; threshold_percentMAX=0.60, plot_window_inS=[-.1 .1])
#     try
#         data, fs = get_videos_audiodata(fname)
#         @debug (basename(fname), extrema(data))

#         trigger = findTrigger(data, fs; threshold_percentMAX=threshold_percentMAX, plot_window_inS=plot_window_inS)
        
#         isnothing(plot_window_inS) && return (trigger-1)/fs
#         # @debug trigger
#         plot!(trigger; title=basename(fname))
#     catch err
#         missing
#     end
# end

# function findAudioBlip(fname; threshold_percentMAX=0.81, plot_window_inS=[-.1 .1])
#     try
#         if occursin(".wav", fname)
#             data, fs = wavread(fname; subrange=0)
#             try 
#                 data, fs = wavread(fname; subrange=3*fs)
#             catch err
#                 @warn "file too short"
#                 data, fs = wavread(fname)
#             end
#         else
#             data, fs, nbits, opt, timestamp = readAudio(fname);
#             data = data[1:fs*3,:]
#         end

#         trigger = findTrigger(data, fs; threshold_percentMAX=threshold_percentMAX, plot_window_inS=plot_window_inS)

#         isnothing(plot_window_inS) && return (trigger-1)/fs
#         plot!(; title=basename(fname))

#         # plot(signal(data,fs)[:,end][1.5s:2.2s]; title=basename(fname))
#         # plot!(signal(data,fs)[:,end][1.5s:2.2s] .|> abs; title=basename(fname))
#         # ylims!(-0.001, .101)
#         # vline!([500])
#     catch err
#         missing
#     end
# end

# function findBlip_bothVidAudio(folname; filename_only=false, vidtype=r".mkv|.MP4|.avi|.mp4", autype=r".wav|.mat|.flac.mp3")
#     # folname = "/Volumes/My Passport/Bahamas_2022/2022.06.25/0002"
#     flist = readdir(folname; join=true) |> skiphiddenfiles
    
#     vids = filter( x -> occursin(vidtype, x), flist)
#     aus = filter( x -> occursin(autype, x), flist)

#     if filename_only
#         trigger_times = [basename.(vids) findVidAudioBlip.(vids; plot_window_inS=nothing) get_duration_smart.(vids) basename.(aus) findAudioBlip.(aus; plot_window_inS=nothing) get_duration_smart.(aus)]
#     else
#         trigger_times = [findVidAudioBlip.(vids; plot_window_inS=nothing) findAudioBlip.(aus; plot_window_inS=nothing)]
#     end
#     return trigger_times
# end


# function TimeZones.ZonedDateTime(dt::AbstractString)
#     plus_index = findlast('+', dt)
#     tz_ind = isnothing(plus_index) ? findlast('-',dt) : plus_index
#     return ZonedDateTime( DateTime(dt[1:tz_ind-1]), TimeZone(dt[tz_ind:end]))
# end

# function TimeZones.ZonedDateTime(dt::String)
# # function zdt(dt::String)
#     if '.' in dt
#         if length(dt) > 27
#             plus_index = findlast('+', dt)
#             tz_ind = isnothing(plus_index) ? findlast('-',dt) : plus_index

#             ts =  split(dt, ['-','T',':','.','+'])
#             return ZonedDateTime( parse.(Int, ts[1:7])..., TimeZone(dt[tz_ind:end]))
#         end
#     end

#     # ZonedDateTime( (parse.(Int,split(ts[1:23], ['-','T',':','.'])))..., TimeZone("-0400"))
# end


# findBlip_bothVidAudio(folname; filename_only=true) |> showall



#~ #FUTURE
# skipsomething(itr) = SkipSomething(itr)

# struct SkipSomething{T}
#     x::T
# end
# IteratorSize(::Type{<:SkipSomething}) = SizeUnknown()
# IteratorEltype(::Type{SkipSomething{T}}) where {T} = IteratorEltype(T)
# eltype(::Type{SkipSomething{T}}) where {T} = nonsomethingtype(eltype(T))

# function iterate(itr::SkipSomething, state...)
#     y = iterate(itr.x, state...)
#     y === nothing && return nothing
#     item, state = y
#     while f(item)
#         y = iterate(itr.x, state)
#         y === nothing && return nothing
#         item, state = y
#     end
#     item, state
# end

# IndexStyle(::Type{<:SkipSomething{T}}) where {T} = IndexStyle(T)
# eachindex(itr::SkipSomething) =
#     Iterators.filter(i -> @inbounds(itr.x[i]) |> f, eachindex(itr.x))
# keys(itr::SkipSomething) =
#     Iterators.filter(i -> @inbounds(itr.x[i]) |> f, keys(itr.x))
# @propagate_inbounds function getindex(itr::SkipSomething, I...)
#     v = itr.x[I...]
#     f(v) && throw(SomethingException("the value at index $I is something"))
#     v
# end