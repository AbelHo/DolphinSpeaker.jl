include("detector.jl")
include("plotting.jl"); plotlyjs()
include("tabulate_data.jl")

set_device__soundtrap()

using Glob, CSV, DataFrames
detection_folder = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/20240116_184938" #"/home/spin/Documents/data_res/marecet/temp/Datai_20250724/detection_files_20250724/20240116_184938"
verification_folder = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/verify_20240116"

vfol = glob("*train-only.txt", verification_folder)
dfol = glob("raven_*train-only.txt", detection_folder)
jldfol = glob("*.jld2", detection_folder)

jldfol_dt = DEFAULT_fname2timestamp_func.(jldfol)
dfol_dt = DEFAULT_fname2timestamp_func.(dfol)

# d_summary = CSV.read("/media/spin/anas2/data/marecet/datai/summary.csv", DataFrame)

autocor_correlated = []
for (i, df) in enumerate(CSV.read.(vfol, DataFrame))
    @debug "Processing file: $(dfol[i])"
    @info dfol[i]|>basename, size(df, 1)

    df = CSV.read(vfol[i], DataFrame)
    dt = DEFAULT_fname2timestamp_func(vfol[i])
    jldfname = jldfol[findfirst( ==(jldfol_dt[findfirst(==(dt), jldfol_dt)]), jldfol_dt)]
    # dfname = dfol[findfirst( ==(dfol_dt[findfirst(==(dt), dfol_dt)]), dfol_dt)]

    df = filter_check(df)
    # df_ori = CSV.read(dfname, DataFrame)
    # d = load(jldfname)
    # aufname = d_summary[findfirst(==(dt), d_summary.datetime),:filepath]

    autocor_clips, plot_dir = analyze_and_plot_clips(
        string(dt),
        "/media/spin/anas2/data/marecet/datai/summary.csv",
        "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648",
        impulsive_band_pass,
        [-25 102];
        plot_dir_prefix = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/autcor_clicktrains/" ,
        output_types = ["png", "html"], flag_extra_plot = true
    )
    push!(autocor_correlated, autocor_clips)
    
    # threshold_autocor, threshold_n_autocor

    # df = df[!, Not([:File, :Begin File, :End File])]
    # df = rename(df, :Begin => :Begin_Time, :End => :End_Time)
    # df = select(df, [:Begin_Time, :End_Time, :Label])
    # df.Label = replace.(df.Label, r"^$" => "unknown")
    # CSV.write(joinpath(verification_folder, "train-only-$i.txt"), df)
end
d = load(jldfol[7])

dd = filter_check.(vfol) .|> x -> findall( ==(true), x.verification)
autocor_correlated |> showall
dd |> showall

corrects = DataFrame( :datetime=>vfol .|> DEFAULT_fname2timestamp_func, :ct_yes=>true)
counts_summary = CSV.read("/home/spin/Documents/data_res/marecet/temp/Datai_20250724/detection_files_20250724/counts_all.csv", DataFrame)
dff = sort(outerjoin(counts_summary, corrects, on=:datetime), :datetime)
d_merge = filter(x-> x.num_impulsetrain > 0, dff)

function filter_check(df::DataFrame)
    # rename checked by SB to check
    if "checked by SB" in names(df)
        rename!(df, Symbol("checked by SB") => :check)
    else
        rename!(df, Symbol("Checked by SB") => :check)
    end

    df.verification = map(df[!, :check]) do x
        ismissing(x) ? false : occursin("correct", x) || occursin("buzz", x)
    end
    return df
end
filter_check(dfname::String) = CSV.read(dfname, DataFrame) |> filter_check



filter(x-> x.verification, df).Selection |> showall



file1="/home/spin/Downloads/Telegram Desktop/raven_8338_240214201545_t2_375081633500817_d200_cps40_0_train_only.txt"
file2="/home/spin/Downloads/Telegram Desktop/raven_8338_prev_convention.txt"


#~ analyse autocorrelation
d_summary = CSV.read("/media/spin/anas2/data/marecet/datai/summary.csv", DataFrame)

autocor_correlated = Array{Array}(undef, size(d_summary, 1)); dts_list = d_summary.datetime; #Array{DateTime}(undef, size(d_summary, 1))
# for dts in d_summary.datetime
res_dir = "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/autcor_clicktrains/"
mkpath(res_dir)
csv_out = joinpath(res_dir, "autocorrelated.csv")
open(csv_out, "w") do io
    println(io, "datetime,impulsetrain_correlated");
    # for (i, dts) in enumerate(d_summary.datetime)
    for i in 1:size(d_summary, 1)
        dts = d_summary.datetime[i]
        autocor_clips, clips_plot_dir, threshold_autocor, threshold_n_autocor = analyze_and_plot_clips(
            dts,
            "/media/spin/anas2/data/marecet/datai/summary.csv",
            "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648",
            impulsive_band_pass,
            [-25 102];
            plot_dir_prefix = res_dir,
            output_types = ["png", "html"], flag_extra_plot = true
        )
        autocor_correlated[i] = autocor_clips
        # push!(autocor_correlated, autocor_clips)
        # push!(dts_list, dts)
        println(io, string(dts) * "," * join(autocor_clips, ";")); flush(io);
    end
end


CSV.write(
    joinpath(res_dir, "autocorrelated2.csv"),
    DataFrame(datetime=dts_list, 
    num_impulsetrain_correlated=length.(autocor_correlated),
    impulsetrain_correlated=array_tostring.(autocor_correlated))
)

# show data available datetime
PlotlyJS.scatter(
    x = string.(dts_list[1:end-1]),
    y = ones(length(dts_list)-1),
    mode = "markers",
    marker = attr(size=8, color="blue", opacity=0.25)
) |> PlotlyJS.plot

df_new = DataFrame(datetime=dts_list, 
    num_impulsetrain_correlated=length.(autocor_correlated),
    impulsetrain_correlated=array_tostring.(autocor_correlated))
df = CSV.read("/home/spin/Documents/data_res/marecet/temp/Datai_20250724/detection_files_20250724/counts_all.csv", DataFrame)
df2 = outerjoin(df, df_new; on=:datetime)
df2 = df2[1:end-1,:]
select!(df2, Not([:num_tonal]))
rename!(df2, :num_impulsetrain_correlated => :num_tonal)
combine_autocor_csv = joinpath(res_dir,"combine_autocor-tonal.csv");
CSV.write(combine_autocor_csv, df2)
detectionsfiles2plot2(combine_autocor_csv; res_dir=res_dir, plottype=PlotlyJS.bar,
    latitude=6.3675, longitude=99.79774, add_daynight=true)

# plot each separately #clickT_correl
df_filt = filter(x-> x.datetime < Date("2024-02-01"), df2)
detectionsfiles2plot2(df_filt; res_dir=res_dir, plottype=PlotlyJS.bar,
    latitude=6.3675, longitude=99.79774, add_daynight=true)
df_filt = filter(x-> Date("2024-04-04") > x.datetime > Date("2024-02-01"), df2)
detectionsfiles2plot2(df_filt; res_dir=res_dir, plottype=PlotlyJS.bar,
    latitude=6.3675, longitude=99.79774, add_daynight=true)
df_filt = filter(x-> x.datetime > Date("2024-04-01"), df2)
detectionsfiles2plot2(df_filt; res_dir=res_dir, plottype=PlotlyJS.bar,
    latitude=6.3675, longitude=99.79774, add_daynight=true)


function analyze_and_plot_clips(
    target,#::AbstractString,
    summary_fname::AbstractString,
    result_directory::AbstractString,
    impulsive_band_pass::AbstractVector,
    window_extract,
    threshold_autocor::Real=4,
    threshold_n_autocor::Int=3
    ;
    plot_dir_prefix::AbstractString="temp/clips_train_",
    output_types = ["html"],
    flag_extra_plot = false
)
    # Find closest row and get audio file
    result = find_closest_row(summary_fname, target)
    aufname = result.filepath
    @info "Processing audio file: $aufname"

    # Find result path
    readdirjoin(x) = readdir(x; join=true)
    respath = readdir(result_directory; join=true) |>
        filter(isdir) .|> readdirjoin .|>
        filter(endswith(".jld2")) .|>
        filter(contains(splitext(basename(aufname))[1])) |>
        filter(!isempty) |> first

    @info "loading result from: $respath"
    res = load(respath)
    res = dict2namedtuple(res)

    if res.res_impulsetrain.train_start |> isempty
        @warn "No impulsive train detected in the result."
        return Int[], "", threshold_autocor, threshold_n_autocor
    end

    # Read and filter audio
    data, fs, _, _, timestamp = readAudio(aufname; fname2timestamp_func=fname2dt_soundtrap)
    data_filt = filter_simple(data, impulsive_band_pass; fs=fs)

    # Extract clips
    clips_fixed = extract_clips(data_filt, window_extract .+ res.res_impulsetrain.pind_good, NaN; flag_matrix=true)
    clips = extract_clips(data_filt, [res.res_impulsetrain.train_start res.res_impulsetrain.train_end], NaN; flag_matrix=false)

    # Prepare output directory
    clips_plot_dir = plot_dir_prefix * Dates.format(timestamp, "yyyymmdd_HHMMSS")
    mkpath(clips_plot_dir)

    train_start_ind = [res.res_impulsetrain.train_start_ind... length(res.res_impulsetrain.pind_good)+1]
    nfunc(x) = sqrt(sum(abs2.(x)))
    autocor_clips = Int[]

    for i in eachindex(clips)
        a = plot(signal(clips[i],fs); title=string(i))
        b = specgram(clips[i]; fs=fs, colorbar=nothing, nfft=128, crange=80)
        c = specgram(clips[i]; fs=fs, colorbar=nothing, nfft=round(Int,fs*.01)|>nextfastfft )
        selections = train_start_ind[i]:train_start_ind[i+1]-1
        snip = clips_fixed[:, selections]
        d = plot_time_fft(snip, fs; legend_position=:outerbottom, labels=reshape(string.(selections),1,length(selections)))

        correls = map(snip|>eachcol) do ref
            map(x-> mfilter(norm_max(ref; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy
        end

        e = Plots.bar(selections, map(x-> mfilter(norm_max(x; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy)
        Plots.hline!([threshold_autocor], label="Threshold", color=:red)

        for ftype = output_types
             savefig(joinpath(clips_plot_dir, "clip_$(i).$ftype"))
             if flag_extra_plot
                Plots.plot(d, e; layout=(2,1), legend_position=:outerbottom)#, size=(1200,800), title=string(i))
                savefig(joinpath(clips_plot_dir, "all_$(i).$ftype"))
             end
        end
        # savefig(joinpath(clips_plot_dir, "clip_$(i).html"))

        if count(map(x-> mfilter(norm_max(x; norm_func=nfunc), norm_max(x; norm_func=nfunc)),  eachcol(snip)) .|> energy .> threshold_autocor) > threshold_n_autocor
            push!(autocor_clips, i)
        end
    end
    autocor_clips |> show
    for ftype = output_types
            make_clip_index_html(clips_plot_dir; outname=basename(clips_plot_dir)*"_$ftype.html", output_type=ftype)
    end
    # make_clip_index_html(clips_plot_dir; outname=basename(clips_plot_dir)*".html", output_types)
    return autocor_clips, clips_plot_dir, threshold_autocor, threshold_n_autocor
end

target_dt = "2024-02-14T20:15"
autocor_clips, plot_dir = analyze_and_plot_clips(
    target_dt,
    "/media/spin/anas2/data/marecet/datai/summary.csv",
    "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648",
    impulsive_band_pass,
    [-25 102];
    output_types = ["png", "html"]
)

resname = "/media/spin/anas2/data/marecet/labels/raven_8338.240214201545___t2.375081633500817_d200__cps40.0_nbhf_100000.txt" # "/media/spin/anas2/data_res/dolphin/marecet/datai/new_20250724/raw/173648/verify_20240116/raven_8338_240214201545_t2_375081633500817_d200_cps40_0_train_only.txt"
df = CSV.read(resname, DataFrame, delim='\t')
filter(row -> !ismissing(row["checked by SB"]), df)
df_skipped = filter(row -> !ismissing(row["checked by SB"]), df)
convert_raven_to_audacity_one(rename(rename(df_skipped, Dict(:Notes => :tmp, Symbol("checked by SB") => :Notes)), Dict(:tmp => :old_notes)),
    splitext(replace(resname, "raven"=>"audacity"))[1] * "_removed.txt")

audacity_label(df[autocor_clips, ["Begin Time (s)", "End Time (s)"]], splitext(replace(resname, "raven"=>"audacity"))[1] * "_autocor-n3t4.txt")

# CSV.write(splitext(resname)[1] * "_removed.txt", df_skipped)