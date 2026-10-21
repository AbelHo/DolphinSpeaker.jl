# Checks the output of run_contiguous_folders_multitdoa (src/run_multitdoa.jl) in one result folder.
# usage: julia --project=. improvements/tdoa_localization_2026-10-21/scripts/check_multitdoa_outputs.jl <res_dir> [reference_webui_csv]
#   res_dir: the folder that holds detection_pixels__color255_webui__Impulse*.csv
#   reference_webui_csv: detection_pixels__color255_webui.csv of an earlier run_contiguous_folders run on the same recording
using CSV, DataFrames, Printf

res_dir = ARGS[1]
ref_csv = length(ARGS) >= 2 ? ARGS[2] : nothing

std = CSV.read(joinpath(res_dir, "detection_pixels__color255_webui.csv"), DataFrame)
files = sort(filter(f -> startswith(f, "detection_pixels__color255_webui__Impulse"), readdir(res_dir)))
variants = Dict(replace(f, "detection_pixels__color255_webui__" => "", ".csv" => "") => CSV.read(joinpath(res_dir, f), DataFrame) for f in files)

# max absolute difference per numeric column and number of differing shapes between two webui tables
function compare(a, b)
    size(a) == size(b) || return "different size $(size(a)) vs $(size(b))"
    num = [maximum(abs, a[!, c] .- b[!, c]) for c in ["frame", "px", "py", "radius", "r", "g", "b", "a"]]
    return @sprintf("max |diff| frame %.3g, px %.3g, py %.3g, radius %.3g, rgb %.3g, a %.3g; shape differs in %d rows",
        num[1], num[2], num[3], num[4], maximum(num[5:7]), num[8], count(a.shape .!= b.shape))
end

println("variants found: ", join(sort(collect(keys(variants))), ", "))
n_tonal = count(!=("circle"), std.shape)   # impulse rows are circles in the standard file
n_impulse = nrow(std) - n_tonal
println("\nstandard file: $(nrow(std)) rows = $n_impulse impulse + $n_tonal tonal")
for (key, df) in sort(collect(variants); by=first)
    tonal_rows = df[df.shape .== "square", :]
    println(@sprintf("%-32s rows %7d  triangles %6d  squares(tonal) %6d  %s", key, nrow(df), count(==("triangle"), df.shape), nrow(tonal_rows),
        startswith(key, "Impulse5g") || startswith(key, "Impulse6g") ? "(gated twin)" : ""))
end

# the ungated variants have the same click set: equal impulse row counts
ungated = [k for k in keys(variants) if !occursin("g_", k)]
counts = Dict(k => nrow(variants[k]) - count(==("square"), variants[k].shape) for k in ungated)
println("\nimpulse rows per ungated variant: ", join(["$k=$(counts[k])" for k in sort(ungated)], ", "))
println("equal impulse row counts across ungated variants: ", length(unique(values(counts))) == 1)

k1 = first(filter(k -> startswith(k, "Impulse1_"), keys(variants)))
println("\nImpulse1 file vs standard file in the same folder: ", compare(variants[k1], std))
if !isnothing(ref_csv)
    ref = CSV.read(ref_csv, DataFrame)
    println("standard file vs reference $(basename(dirname(ref_csv))): ", compare(std, ref))
end

summary_csv = joinpath(res_dir, "multitdoa_summary.csv")
if isfile(summary_csv)
    println("\nmultitdoa_summary.csv:")
    show(stdout, MIME("text/plain"), CSV.read(summary_csv, DataFrame)); println()
end
