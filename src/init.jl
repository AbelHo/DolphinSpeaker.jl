version = "2025-09-04T13:00"
@info "version v$version =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-"
include("utils.jl")
export showall, skiphiddenfiles, process_files

include("audio.jl")
export mat2flac, mat2wav, mat2flac_check, bin2flac, bin2flac_check

include("synchronization.jl")
export findBlip_bothVidAudio, flac2signal, multisync

include("run_example.jl")
export process_one_set, process_folder, process_dir

# from dsp.jl
export extrema_in_file

include("beampattern.jl")
export stack_audio_videos, run_func_fileauto


DETECTION_TYPES = 1:2
impulsive_autothreshold_median_ratio = 7
OVERLAY_DEFAULT_ALPHA = 0.2
pt_config = [((1,0,0),30), ((1,1,0),20), ((0,1,0),15), ((1,0,1),10), ((1,1,1),7)]
process_folder("/media/spin/anas2/data/calf/new_2025_freeplay_report"; outfolder="/media/spin/anas2/data_res/dolphin/calf/res_20250904/rope2")
