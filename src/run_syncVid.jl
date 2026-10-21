include("run_example.jl")
set_device__ophk_acoustic_D3()


folder = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
res_dir = "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_nooverlay2"
res = run_contiguous_folders(folder; res_dir = res_dir, flag_overlayvideo=false)