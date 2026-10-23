# Direction candidates + tracks for the New_device multitdoa result folder (post-processing; see src/multiarrival.jl).
# usage: julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/run_multiarrival_newdevice.jl [res_folder]
src = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src, "multiarrival.jl"))
res_folder = length(ARGS) >= 1 ? ARGS[1] :
    "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device"
@time run_multiarrival(res_folder; out_dir=joinpath(@__DIR__, "..", "results"))
