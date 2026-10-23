# Audacity label files for the multi-arrival methods, on the same combined timeline as the pipeline's
# combined_impulse__audacity.txt (files concatenated in name order, no A/V sync delay). One point label per click:
#   combined_Impulse9sicwt__audacity.txt                track rule (round 2):  "az el direct|echo|untracked"
#   combined_Impulse10sicw_longesttrack-inherit__audacity.txt  longest track + followers inherit (round 3): "az el [follower]"
#   combined_Impulse11sicw_sources_n<N>__audacity.txt   dolphin separation (round 4): "dolphin k|unassigned role az el"
# Each is written only if its input table exists (run the corresponding step of run_on_recording.sh first).
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/export_audacity.jl <res_folder> [n]
using CSV, DataFrames, JLD2, Printf
res_folder = ARGS[1]
n_given = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
folder = basename(rstrip(res_folder, '/'))
res = joinpath(@__DIR__, "..", "results", "data")
# cumulative start time of each audio file (in the pipeline's order: sorted by name)
jl = sort(filter(f -> occursin("__cps", f) && endswith(f, ".jld2") && !endswith(f, "_angles.jld2"), readdir(res_folder)))
offset = Dict{String,Float64}(); cum = 0.0
for f in jl
    j = jldopen(joinpath(res_folder, f)); ri = j["res_impulse"]; fs = j["fs"]; close(j)
    offset[basename(ri.aufname)] = cum
    global cum += ri.len_data / fs
end
function write_labels(name, t, labels)
    fname = joinpath(res_folder, "combined_$(name)__audacity.txt")
    open(fname, "w") do io
        for (ti, l) in zip(t, labels); @printf(io, "%.6f\t%.6f\t%s\n", ti, ti, l); end
    end
    println("wrote ", fname, " (", length(t), " labels)")
end
tt(df) = df.t_s .+ [offset[f] for f in df.file]

p = joinpath(res, "click_features_sicw_$folder.csv")
if isfile(p)
    f = CSV.read(p, DataFrame)
    az = [f[i, "az_$(f.selected_track[i])"] for i in 1:nrow(f)]; el = [f[i, "el_$(f.selected_track[i])"] for i in 1:nrow(f)]
    write_labels("Impulse9sicwt", tt(f), [@sprintf("%.1f %.1f %s", az[i], el[i], f.track_status[i]) for i in 1:nrow(f)])
end
p = joinpath(res, "real_rules_sicw_default_longesttrack-inherit_$folder.csv")
if isfile(p)
    f = CSV.read(p, DataFrame)
    write_labels("Impulse10sicw_longesttrack-inherit", tt(f), [@sprintf("%.1f %.1f%s", f.az[i], f.el[i], f.triangle[i] ? " follower" : "") for i in 1:nrow(f)])
end
p = joinpath(res, "real_sources_sicw_tuned_$folder.csv")
if isfile(p)
    f = CSV.read(p, DataFrame)
    lab(l) = l == 0 ? "unassigned" : "dolphin $l"
    write_labels("Impulse11sicw_sources_n$(n_given)", tt(f), [@sprintf("%s %s %.1f %.1f", lab(f.label_n[i]), f.role_n[i], f.az_out[i], f.el_out[i]) for i in 1:nrow(f)])
end
