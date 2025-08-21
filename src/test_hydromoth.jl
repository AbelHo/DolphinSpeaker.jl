res_dir = "/media/spin/anas/data_res/dolphin/marecet/temp/aqua/julia2"
aufname = "/media/spin/anas/data/marecet/hydromoth"
a=readdir(aufname; join=true)

fname2dt_date_time_hydromoth(aufname) = DateTime(basename(aufname)[4:18], dateformat"yyyymmdd_HHMMSS")
DEFAULT_fname2timestamp_func = fname2dt_date_time_hydromoth

res = detect_impulseNtonal.(readdir(aufname; join=true), res_dir; threshold_tonal=1, freq_maxbandwidth=3000, freq_width_db=10);