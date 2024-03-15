using Dates
# rx = 0.14722/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) ) #aspod 2
# rx = 0.4/sqrt(3) .* exp.(im.* deg2rad.([30 -90 150]) ) # calf_hk
# rx_vect = [real(rx); imag(rx); zeros(1,3)]

device_name = "aspod2"
device_name = "punnet"
device_name = "calf_hk"
# device_name = "calf_hk_clicker"
# device_name = "aspod2"

ref_channel = 3
c = 1540 # m/s speed of sound

#~ impulsive sound detector parameters
dist_impulsive = 800# 15000#pinger
window_impulsive = -50:200; #-5000:10000 #-50:300
impulsive_band_pass = [60, Inf] #fs/2*.98]
impulsive_autothreshold_median_ratio = 4;

#~ tonal sound detector parameters
window_tonal = 0.01
nfft_inS = 0.01
tonal_band_pass = [5000, 24000] # [5000, 25000] # 

butterworth_size = 4
percent_quiet = 0.01
freq_maxbandwidth = 10000
freq_width_db = 3
winlen=0.01s # tdoa estimation window length

#~ boat detector
threshold_boat = 8e-6
band_pass_boat = [1 500]

#~ bin file parameters
FILE_device_ID = ""
FILE_location_ID = ""
FILE_gain_setting = ""
FILE_device_location = ""

#~ audio
DEFAULT_fname2timestamp_func = nothing

#~ bin file parameters
DEFAULT_bin_channels = 4

#~ video overlay parameters
pt_config = [((1,0,0),63), ((1,1,0),58), ((0,1,0),53), ((1,0,1),48), ((1,1,1),33)]
DETECTION_TYPES = 1:2
# pt_config = [((1,1,0),25), ((1,0,0),30), ((0,1,0),20), ((1,0,1),15), ((1,1,1),10)]

#~ audio merge parameters
MERGE_VID_AU_DYNAMIC_NORM = false


calf_timestamp_func(fname) = DateTime(basename(fname)[1:17], DateFormat("yyyymmdd_HH.MM.SS"))

function set_device__hk_clicker()
    global impulsive_band_pass = [1000, Inf] #fs/2*.98]
    global threshold_impulsive = nothing# .003#calf_hk
    global impulsive_autothreshold_median_ratio = 12

    global threshold_tonal = nothing #-110#calf_hk -15#aspod2
    global click_train_minlen = 2#calf 7 #5 20
    global click_train_check_interval = 8#calf 1 #.01
    global window_impulsive = -300:700 #-50:100
    global dist_impulsive = 15000# 15000#pinger

    global rx = 0.4/sqrt(3) .* exp.(im.* deg2rad.([150 -90 30]) ) # calf_hk
    global imsize=(2160, 3840)
    global fov_angle=[54,34,0]
end

function set_device__calf_hk()
    global ref_channel = 3

    global impulsive_band_pass = [5000, 180_000] #fs/2*.98]
    global threshold_impulsive = nothing #.1# .003#calf_hk
    global impulsive_autothreshold_median_ratio = 9;#4 use for for full detection but less accurate localization;
    global dist_impulsive = 200
    # global # dist_impulsive = 400#800# 15000#pinger

    global click_train_minlen = 5#calf 7 #5 20  #[n+1]clicks within (click_train_check_interval)seconds
    global click_train_check_interval = 0.1#calf 1 #.01  #seconds has at least (click_train_minlen) clicks

    global window_impulsive = -130:130 #-299:700; #-50:200; #-5000:10000 #-50:300

    global threshold_tonal = nothing #-110#calf_hk -15#aspod2
    global tonal_band_pass = [2500, 24000]

    global rx = 0.4/sqrt(3) .* exp.(im.* deg2rad.([150 -90 30]) ) # calf_hk
    global imsize=(2160, 3840)
    global fov_angle=[54,34,0]

    global threshold_boat = Inf #1e3
    global band_pass_boat = [1 500]
    # global default_getTDOA_func = get_tdoa_raw_MaxEnergyRefChannel

    global DEFAULT_fname2timestamp_func = calf_timestamp_func

    # bin file parameters
    global FILE_device_ID = "ophk_acoustic_D1.2.0"
    global FILE_location_ID = "MMBRC_pool2"
    global FILE_gain_setting = "0011"

    #~ overlay parameters
    global DETECTION_TYPES = 1:1
end

function set_device__aspod2()
    global band_pass = [5000, 24000]; winlen=0.01s;
    global threshold_impulsive = .005*32767#aspod2 .003#calf_hk #0.01 # 0.1#pinger/clickler #0.2
    global threshold_tonal = -20#aspod2 -110#calf_hk 
    global rx = 0.14722/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) ) #aspod 2
    global imsize=(2160, 3840)
    global fov_angle=[62.61721188568244,35.793211268714096,71.6855447884958]
end

function set_device__rwsnus()
    global DEFAULT_fname2timestamp_func = fname2dt_dashdot

    # bin file parameters
    global FILE_device_ID = "rwsnus_array_1.0.0"
    global FILE_location_ID = "rws"
    global FILE_gain_setting = "0000"
end

function set_device__LS1()
    global impulsive_band_pass = [5000, Inf] #fs/2*.98]
    global ref_channel = 1
    global tonal_band_pass = [900 Inf]
    global freq_width_db = 3 #3
    global freq_maxbandwidth = 1910
    global nfft_inS = 0.02

    global DEFAULT_fname2timestamp_func = fname2dt_ls1

    global rx_vect = [0;0;0]

    # bin file parameters
    # global FILE_device_ID = "rwsnus_array_1.0.0"
    # global FILE_location_ID = "rws"
    # global FILE_gain_setting = "0000"
end

#~ default
    impulsive_band_pass = [1000, Inf] #fs/2*.98]
    threshold_impulsive = nothing #.1# .003#calf_hk
    impulsive_autothreshold_median_ratio = 4;#10;
    # dist_impulsive = 400#800# 15000#pinger

    click_train_minlen = 5#calf 7 #5 20  #[n+1]clicks within (click_train_check_interval)seconds
    click_train_check_interval = 0.1#calf 1 #.01  #seconds has at least (click_train_minlen) clicks

    window_impulsive = -300:700; #-50:200; #-5000:10000 #-50:300

    threshold_tonal = nothing #-110#calf_hk -15#aspod2
    tonal_band_pass = [2500, 24000]

    rx = 0.4/sqrt(3) .* exp.(im.* deg2rad.([150 -90 30]) ) # calf_hk
    imsize=(2160, 3840)
    fov_angle=[54,34,0]

    threshold_boat = Inf #1e3
    band_pass_boat = [1 500]



if device_name=="aspod2"
    band_pass = [5000, 24000]; winlen=0.01s;
    threshold_impulsive = .005*32767#aspod2 .003#calf_hk #0.01 # 0.1#pinger/clickler #0.2
    threshold_tonal = -20#aspod2 -110#calf_hk 
    rx = 0.14722/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) ) #aspod 2
    imsize=(2160, 3840)
    fov_angle=[62.61721188568244,35.793211268714096,71.6855447884958]
elseif device_name=="calf_hk"
    set_device__calf_hk()
    # impulsive_band_pass = [1000, Inf] #fs/2*.98]
    # threshold_impulsive = nothing #.1# .003#calf_hk
    # impulsive_autothreshold_median_ratio = 4;#10;
    # # dist_impulsive = 400#800# 15000#pinger

    # click_train_minlen = 5#calf 7 #5 20  #[n+1]clicks within (click_train_check_interval)seconds
    # click_train_check_interval = 0.1#calf 1 #.01  #seconds has at least (click_train_minlen) clicks

    # window_impulsive = -300:700; #-50:200; #-5000:10000 #-50:300

    # threshold_tonal = nothing #-110#calf_hk -15#aspod2
    # tonal_band_pass = [2500, 24000]

    # rx = 0.4/sqrt(3) .* exp.(im.* deg2rad.([150 -90 30]) ) # calf_hk
    # imsize=(2160, 3840)
    # fov_angle=[54,34,0]

    # threshold_boat = Inf #1e3
    # band_pass_boat = [1 500]

elseif device_name=="punnet"
    # threshold_impulsive = .005*32767#aspod2 .003#calf_hk #0.01 # 0.1#pinger/clickler #0.2
    threshold_tonal = -20#aspod2 -110#calf_hk 
    threshold_impulsive = 0.01 #0.21 #0.215
    click_train_minlen = 7 #[n+1]clicks within (click_train_check_interval)seconds
    click_train_check_interval = 1 #seconds has at least (click_train_minlen) clicks
    # click_train_check_interval = .2;click_train_minlen = 4;threshold_impulsive = nothing
    rx = 0.5/sqrt(3) .* exp.(im.* deg2rad.([90 -30 210]) ) #0.5/sqrt(3) .* exp.(im.* deg2rad.([-30 90 -150]) )
elseif device_name=="calf_hk_clicker"
    set_device__hk_clicker()
    # impulsive_band_pass = [1000, Inf] #fs/2*.98]
    # threshold_impulsive = .1# .003#calf_hk
    # threshold_tonal = nothing #-110#calf_hk -15#aspod2
    # click_train_minlen = 2#calf 7 #5 20
    # click_train_check_interval = 5#calf 1 #.01
    # window_impulsive = -300:700 #-50:100
    # dist_impulsive = 15000# 15000#pinger

    # rx = 0.4/sqrt(3) .* exp.(im.* deg2rad.([150 -90 30]) ) # calf_hk
    # imsize=(2160, 3840)
    # fov_angle=[54,34,0]
end

rx_vect = [real(rx); imag(rx); zeros(1,3)]

# folname = "/Users/abel/Documents/data/aspod/field/bahamas_2022"
# res_dir = "/Users/abel/Documents/data_res/aspod/real/bahamas_2022_10_test"

# if !Sys.islinux()
#     folname = "/Volumes/My Passport/Bahamas_2022" #"/Users/abel/Documents/data/aspod/field/bahamas_2022" #"/Volumes/dd/Bahamas_2022"
#     res_dir = "/Volumes/My Passport/res/bahamas_2022_4_test" #"/Users/abel/Documents/data_res/aspod/real/bahamas_2022_redo" #"/Volumes/dd/result/bahamas_2022" #"/Users/abel/Documents/data_res/aspod/real/bahamas_2022"
# else
#     res_dir = "/run/user/1000/gvfs/smb-share:server=10.246.128.21,share=others/sg_res/Aspod/bahamas_2022" #/run/user/1000/gvfs/smb-share:server=10.246.128.21,share=others/sg_res/Aspod/temp7" #"/Users/abel/Documents/data_res/aspod/real/bahamas_2022"
#     folname = "/run/user/1000/gvfs/smb-share:server=10.246.128.21,share=data/Aspod/Aspod4/Bahamas_2022" #"/media/sousa/Space/data/aspod2/Bahamas_2022" #"/media/sousa/My Passport/Bahamas_2022"
# end

set_device__calf_hk()