# include("calibration_aspod4_2.jl")
include("test_run.jl")
include("test_make_video.jl")
using WAV
using DSP
using Plots
using Images
import GLMakie

vidfname = "/Volumes/dd/Bahamas_2022/2022.06.26/0001/Vid_2022-06-26_101510.mkv"
aufname = "/Volumes/dd/Bahamas_2022/2022.06.26/0001/Aud_2022-06-26_101511.wav"
data, fs = wavread(aufname; format="native")
@info "passed"

bp = Bandpass(5000,14000, fs=fs)
bs = Bandstop(155400, 156000, fs=fs) #Bandstop(154000, 157000, fs=fs)
data_filt = mapslices( x -> filtfilt( digitalfilter(bp, Butterworth(4)), x), data, dims=1)
# data_filt = filtfilt( digitalfilter(bs, Butterworth(4)), data[:,3])

d=signal(data, fs)
d_filt=signal(data_filt, fs)
## plot
# plotlyjs()
# psd(data[1:400000,1], fs=fs)
# data_filt = filtfilt( digitalfilter(bs, Butterworth(4)), data[:,3])
# psd!(data_filt[1:400000,1], fs=fs)

#~ angle estimate
# win1 = 1.0s;

function twindow2dird(win1, winlen=0.1s)
    window = win1:(win1+winlen);
    tdoa = [finddelay(d_filt[window,1],d_filt[window,3]); finddelay(d_filt[window,2],d_filt[window,3]); finddelay(d_filt[window,3],d_filt[window,3])];
    tdoa_ori = [finddelay(d[window,1],d[window,3]); finddelay(d[window,2],d[window,3]); finddelay(d[window,3],d[window,3])];
    tdoas = [tdoa tdoa_ori]
    # tdoa2dir(tdoas', rx_vect) .|> rad2deg
    tdoa2dir(tdoa', rx_vect) .|> rad2deg
    # tdoa2dir(tdoa_ori', rx_vect) .|> rad2deg
end


function get_image(vid, time=1.0)
    seek(vid, time)
    img = read(vid)
end

# function get_image2(time, vid)
#     @info time
#     seek(vid, time)
#     img = read(vid)
# end

# using ImageDraw, VideoIO
vid = VideoIO.openvideo(vidfname)
winlen = 0.04s
time_range = 65.5s:winlen:67.2s #27.0s:0.1s:29.0s #1.0s:0.1s:3.0s
# imgstack = get_image.(vid, time_range)
# b=twindow2dird.(1.0s:0.1s:3.0s)
pxs = twindow2dird.(time_range, winlen) .|> angle2pxd
seek(vid, time_range[1].val)
img = read(vid)
display_corners!(img, pxs')

# using Images
# p_list = []
# for i in eachindex(time_range)
#     try
#         push!(p_list,display_corners!(get_image(vid, time_range[i].val), [pxs[i]']))
#     catch err
#         @error (i)
#     end
# end
# i=1;println(time_range[i]);p_list[i]
# i=i+1;println(time_range[i]);p_list[i]

@info "Done audio_analysis"