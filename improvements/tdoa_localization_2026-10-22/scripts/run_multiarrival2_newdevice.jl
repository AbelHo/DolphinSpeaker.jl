# Round 2 on the New_device multitdoa result folder: SIC candidates (min_rel 0.7) + earliest-arrival selection,
# with the pipeline window (tag "sic", -50:50) or the wider window that reaches further back (tag "sicw", -110:50).
# usage: julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/run_multiarrival2_newdevice.jl sic|sicw [res_folder]
src = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src, "multiarrival.jl"))
tag = length(ARGS) >= 1 ? ARGS[1] : "sic"
res_folder = length(ARGS) >= 2 ? ARGS[2] :
    "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device"
window = Dict("sic" => window_impulsive, "sicw" => -110:50)[tag]
@time run_multiarrival2(res_folder; out_dir=joinpath(@__DIR__, "..", "results"), candidate_func=get_srp_candidates_sic,
    window=window, min_rel=0.7, select_rule=:earliest, tag=tag)
