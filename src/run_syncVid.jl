include("run_example.jl")
set_device__ophk_acoustic_D3()
include("run_multitdoa.jl")   # run_contiguous_folders_multitdoa: one webui CSV per TDOA method (Impulse1-6), see reviews/plans/multitdoa-webui-csvs_2026-10-22.md

#~ old
# folder = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
# res_dir = "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_nooverlay2"
# res = run_contiguous_folders(folder; res_dir = res_dir, flag_overlayvideo=false)

#~ new
folder = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
res_dir = "/media/spin/anas2/data/calf/upload/results/tmp_20261022/Two-way_2males_06_10_2025_multitdoa"
# all audio files of the folder (0001-0003) are processed as one contiguous set
res = run_contiguous_folders_multitdoa(folder; res_dir = res_dir, closure_max = 1.0)
