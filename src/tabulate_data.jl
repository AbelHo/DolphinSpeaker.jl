using Dates
using VideoIO
include("media_info.jl")
include("utils.jl")

title = "filename,duration,datetime"#,filesize,filepath"

# output = stdout # open(joinpath(homedir(), "Documents/temp/test.csv"), "w") #stdout
# if typeof(output) == IOStream
#   println(output, "filename,duration,datetime")
# end

# fol = "/Volumes/data-megafauna/data/SD_WN_1/acoustic/2022/02/2022-02-23" # "/Volumes/data/Aspod/Hawaii_2022-09/bumblebee/2022-09-24"
# filetype = "flac" #"MP4"

# fmt_str = "yyyymmdd_HHMMSS" # "yyyy-mm-dd_HH.MM.SS"
# fmt = DateFormat(fmt_str)
# separator = "," # "\t"
# postfix = "" #"_0001"

## fmt_str="epochtime-10h"
"""
Example usage:
tabulate_data("/Users/aa/Documents/data/calf"; filetype=:media, flag_recursive=true, flag_filepath=true, flag_filesize=true, output="/Users/aa/Documents/data/calf/summary_media.csv")
"""
function tabulate_data(fol::String;

  output = stdout, filetype::Union{String,Symbol} = "flac",
  fmt_str::String = "yyyymmdd_HHMMSS",
  fname_fmt = "end",
  separator = ",", postfix::String = "",
  flag_filesize=false,
  flag_normalmedia=true,
  flag_filepath=false,
  flag_lastmodified=true,
  flag_overwrite = false,
  title="filename,duration,datetime",
  fname2timestamp_func = nothing,
  flag_recursive::Bool = false
  )

  title_final = title * (flag_filesize ? ",filesize" : "") * (flag_filepath ? ",filepath" : "") * (flag_lastmodified ? ",last_modified" : "")
  # Prepare output stream
  output_stream = output
  output_is_stream = false
  if typeof(output) == IOStream
    println(output, title_final)
    output_stream = output
    output_is_stream = true
  elseif typeof(output) == String
    if flag_overwrite
      output_stream = open(output,"w")
      println(output_stream, title_final)
    else
      if !isfile(output)
        output_stream = open(output,"w")
        @info("Creating new file: $output")
        println(output_stream, title_final)
      else
        # @warn("File $output already exists. Use flag_overwrite=true to overwrite.")
        output_stream = open(output,"a")
      end
    end
  end

  if length(fmt_str)>12 && fmt_str[1:9]=="epochtime"
    timediff = 3600 * parse(Int, fmt_str[10:12])
    epoch = true
  else
    fmt = DateFormat(fmt_str)
    epoch = false
  end

  function process_dir(dir)
    for (root, dirs, files) in walkdir(dir)
      # If not recursive, only process the top-level directory
      if !flag_recursive && root != dir
        break
      end
      println(stdout, "Files in $root")
      for file in files
        file_ext = splitext(file)[2]
        if filetype == :all
            (filetype==:media && occursin(Regex(join( vcat(autypes,vidtypes), '|')), file_ext |> lowercase)) || 
            (filetype==:audio && occursin(Regex(join( autypes, '|')), file_ext |> lowercase)) || 
            (filetype==:video && occursin(Regex(join( vidtypes, '|')), file_ext |> lowercase)) ||
            file_ext==filetype

          print(output_stream, "$file$separator")
          try
            dur = get_duration(joinpath(root,file))
            dur = ismissing(dur) ? "" : dur
            flag_normalmedia && print(output_stream, dur )
          catch
            @error("cant open this file: \t$file")
            print(output_stream, "-99999999" )
          end
          print(output_stream, "$separator")
          minlen = (length(file_ext)+length(postfix))
          @debug basename(file), minlen
          if length(basename(file)) < minlen + length(fmt_str)
            fname_time = file
          elseif fname_fmt == "end"
            fname_time = file[end-minlen-length(fmt_str)+1:end-minlen]
          else
            fname_time = file[1:length(fmt_str)]
          end
          @debug fname_time

          try
            if epoch
              dt = unix2datetime(parse(Int, fname_time)/1000 + timediff)
            else
              if isnothing(fname2timestamp_func)
                dt = DateTime(fname_time, fmt)
              else
                dt = fname2timestamp_func(file)
              end
            end
            print(output_stream, dt)
          catch
            @error("cant convert this time: \t$fname_time")
            dt = missing
          end
          flag_filesize && print(output_stream, ","*string(filesize(joinpath(root,file))))
          flag_filepath && print(output_stream, ",\""*joinpath(root,file)*"\"")
          flag_lastmodified && print(output_stream, ",\""* (mtime(joinpath(root,file))|>unix2datetime|>string) *"\"")
          print(output_stream, "\n")
        end
      end
      # If not recursive, do not process subdirectories
      if !flag_recursive
        break
      end
    end
  end

  process_dir(fol)

  if !output_is_stream
    close(output_stream)
  end
  # if typeof(output) == IOStream
  #   close(output)
  # end
end

# join_dataframe2(df1, df2) = outerjoin(df1,df2, on=intersect(df1|>names,df2|>names), matchmissing=:equal)

join_dataframe(args::Array{String}) = join_dataframe(map(x->CSV.File(x)|>DataFrame, args)) 
function join_dataframe(args)
  # Base case: if there's only one dataframe, return it as is
  if length(args) == 1
      return args[1]
  end

  # Recursive case: join the first two dataframes, then recursively join the rest
  joined_df = outerjoin(args[1], args[2], on=intersect(names(args[1]), names(args[2])), matchmissing=:equal)
  for df in args[3:end]
      joined_df = outerjoin(joined_df, df, on=intersect(names(joined_df), names(df)), matchmissing=:equal)
  end

  return joined_df
end

using CSV, DataFrames, Dates

# Example usage:
# df = CSV.read("/media/spin/anas2/data/marecet/datai/summary.csv", DataFrame)
# target = DateTime("2024-01-18T23:30:00")
# result = find_closest_row(df, target)
# result.delta
# println(result)
# println("Full file path: ", result.filepath)
function find_closest_row(df::DataFrame, target_dt::DateTime)
    # Parse the datetime column if not already DateTime
    if !(eltype(df.datetime) <: DateTime)
        df.datetime = DateTime.(df.datetime)
    end
    # Find the index of the closest datetime
    idx = findmin(abs.(df.datetime .- target_dt))[2]
    row = df[idx, :]
    delta = row.datetime - target_dt
    # Return as NamedTuple with delta and full path
    @info (row.datetime, delta)
    return merge(NamedTuple(row), (delta=delta|>Second,))
end
find_closest_row(df::DataFrame, target_dt::String) = find_closest_row(df, DateTime(target_dt))
find_closest_row(summary_fname::String, target_dt) = find_closest_row(CSV.read(summary_fname, DataFrame), target_dt)

# df = CSV.read("/media/spin/anas2/data/marecet/datai/summary.csv", DataFrame)
# summary_fname = "/media/spin/anas2/data/marecet/datai/summary.csv"
# target = "2024-05-21T13:00"#"2024-05-18T11:00" #"2024-01-19T15:44"
# result = find_closest_row(summary_fname, target)
# result.delta
# println(result)
# println("Full file path: ", result.filepath)

count_arraysize(x; delim=';') = ismissing(x) ? 0 : 1 + count(==(delim), x)
array_tostring(x; delim=';') = join(x, delim)





# include("tabulate_data.jl");
# infol = "/volume1/data-megafauna/data/SD_WN_1"
# outfol = "/volume1/data4-S2S/Megafauna/results/data_index"
# tabulate_data(infol; output=joinpath(outfol,"summary_$(basename(infol))_new.csv"), flag_filesize=true, flag_filepath=true, flag_normalmedia=true, filetype="flac")
# infol = "/volume1/data-megafauna/data/SD_WN_2"
# tabulate_data(infol; output=joinpath(outfol,"summary_$(basename(infol))_new.csv"), flag_filesize=true, flag_filepath=true, flag_normalmedia=true, filetype="flac")
# infol = "/volume1/data-megafauna/data/SD_WN_3"
# tabulate_data(infol; output=joinpath(outfol,"summary_$(basename(infol))_new.csv"), flag_filesize=true, flag_filepath=true, flag_normalmedia=true, filetype="flac")
# infol = "/volume1/data-megafauna/data/SD_WS_1"
# tabulate_data(infol; output=joinpath(outfol,"summary_$(basename(infol))_new.csv"), flag_filesize=true, flag_filepath=true, flag_normalmedia=true, filetype="flac")
# infol = "/volume1/data-megafauna/temp"
# tabulate_data(infol; output=joinpath(outfol,"summary_tempbin.csv"), flag_filesize=true, flag_filepath=true, flag_normalmedia=false, filetype="bin")
# tabulate_data(infol; output=joinpath(outfol,"summary_$(basename(infol)).csv"), flag_filesize=true, flag_filepath=true, flag_normalmedia=true, filetype="flac")

