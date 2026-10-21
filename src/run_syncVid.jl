include("run_example.jl")
set_device__ophk_acoustic_D3()


folder = ""
res_dir = ""
res = run_contiguous_folders(folder; res_dir = res_dir, flag_overlayvideo=false)