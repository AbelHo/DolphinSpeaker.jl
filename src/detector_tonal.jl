using FFTW
using DSP, SignalAnalysis
using JLD2
using WAV
using Statistics
include("config.jl")

include("audacity.jl") # if required reading from audiofile directly instead of just calling detector

# aufname = "/Volumes/dd/Bahamas_2022/2022.06.27/0000/Aud_2022-06-27_11.16.17.wav"
# vidfname = "/Volumes/dd/Bahamas_2022/2022.06.27/0000/Vid_2022-06-27_11.16.16.mkv"
#res_dir
#fs
#data, assumise last channel is synch channel and not used for processing

amp2pow(ampl) = 20 * log10.(abs.(ampl))
# ft = stft_norm(cww, nfft, 0; fs=fs)
function stft_norm(args...; kwargs...) 
    fs = kwargs[:fs]
    nfft = args[2]
 
    ft = stft(args...; kwargs...) ./nfft * 2
    # return ft
    return (;ampl=ft, 
            time=0:nfft/fs:length(args[1])/fs-nfft/fs, 
            freq=0:fs/nfft:fs/2)
end

function data2ft(data, fs; nfft_inS=nfft_inS, ref_channel=1)
    fs = Int(fs)
    nfft = nextfastfft(nfft_inS*fs)
    # stft_norm(data[:,ref_channel], nfft, 0; fs=fs)
    mapslices( x -> stft_norm(x, nfft, 0; fs=fs), data, dims=1)
end

function plot_ft(ft::NamedTuple)
    plot( ft.freq, maximum(abs.(ft.ampl); dims=2) )
end

function specgram_ft(ft::NamedTuple)
    heatmap(ft.time, ft.freq, 20 * log10.(abs.(ft.ampl)) ; colorbar=false)
end

specgram_flatten(ft) = plot(ft.freq, ft.ampl .|> abs; legend=nothing)
specgram_flatten_pow(ft) = plot(ft.freq, ft.ampl .|> abs2; legend=nothing)
specgram_flatten_log(ft) = plot(ft.freq, 10* log10.(ft.ampl .|> abs2 ); legend=nothing)

# function detect_tonal2(fname; nfft_inS=nfft_inS, band_pass=band_pass)
#     data, fs = wavread(fname)
#     nfft = nextfastfft(nfft_inS*fs)
#     freq_filt = Int(round(band_pass[1]/fs*nfft)):Int(round(band_pass[2]/fs*nfft))

#     mag_ft = abs.(ft.ampl)
#     quietests = minimum(mag_ft; dims=2)

#     fft_max = maximum(mag_ft[freq_filt,:]; dims=1)'[:]
# end

function makie_plot_dual(ft, nfft, ind, ft_all=ft; freq_filt=1:size(ft.ampl,1))
    fig = GLMakie.Figure(;resolution=(1500,1080))
    # ax = GLMakie.Axis(fig[1,2];)
    sig=[]; mask=[];
    if ft_all != ft
        mask = zeros(size(ft.ampl[:,ind])) # mask = zeros(size(ft.ampl)) #mask = similar(ft.ampl) #
        mask[freq_filt,:] .= 1
        foreach( x->push!(sig, irfft(x.ampl[:,ind] .* mask,nfft)), ft_all)
        # foreach( x -> push!(sig2,x), ft)
    elseif length(freq_filt) != size(ft.ampl,1)
        mask = zeros(size(ft.ampl[:,ind])) # mask = zeros(size(ft.ampl)) #mask = similar(ft.ampl) #
        mask[freq_filt,:] .= 1
        sig = [irfft(ft.ampl[:,ind], nfft), irfft(ft.ampl[:,ind] .* mask, nfft)]
    else
        sig = irfft(ft.ampl[:,ind], nfft)
    end
    # @info sig|>size
    # @info 0:1/ft.freq[2]:length(sig[:,1])*1/ft.freq[2]-1/ft.freq[2]
    axs = []#Array{GLMakie.Axis}
    push!(axs, GLMakie.Axis(fig[1,2:8]) )
    last(axs).title = string(ind)
    GLMakie.lines!.(Ref(last(axs)), Ref(0:ft.time[2]/nfft:ft.time[2]-ft.time[2]/nfft), sig; label=string(1:length(sig)))
    GLMakie.tightlimits!(last(axs))
    # GLMakie.lines!(ax, 0:ft.time[2]/nfft:ft.time[2]-ft.time[2]/nfft, sig[:,2])
    
    push!(axs, GLMakie.Axis(fig[1:2,1]) )
    GLMakie.lines!(last(axs), ft.freq, abs.(ft.ampl[:,ind]) )
    GLMakie.lines!(last(axs), ft.freq, (abs.(ft.ampl[:,ind])) .* mask )
    # GLMakie.Axis(fig[1,2]).title = string(ind)

    push!(axs, GLMakie.Axis(fig[2,2:8]) )
    specgram_zoom_makie(ft, ind, last(axs))

    fig
end
# slider_ind = GLMakie.Slider(fig[end,1:2], range=1:size(ft.ampl,2), startvalue=1)
# for ind=50:size(ft[1].ampl,2)
#     makie_plot_dual(ft[1], nfft, ind, ft; freq_filt=freq_filt) |> display
#     readline()
# end


function specgram_zoom(ft, ind)
    xlim_min = ind-20 < 1 ? 1 : ind-20
    heatmap(20*log10.(abs.(ft.ampl)); colorbar=false); xlims!(xlim_min, ind+20); ylims!(0,500); vline!(ind .+ [-.5 .5]; color=:black, width=3)
end

function specgram_zoom_makie(ft, ind, ax)
    # fig = GLMakie.Figure(;resolution=(1400,600))
    xlim_min = ind-20 < 1 ? 1 : ind-20
    GLMakie.heatmap!(ax, 20*log10.(abs.(ft.ampl))'; colorbar=false)
    GLMakie.xlims!(ax, xlim_min, ind+20)
    GLMakie.ylims!(ax, 0,500)
    GLMakie.vlines!(ax, ind .+ [-.5 .5]; color=:black, width=3)
    fig
end
# function plot_fft(args...; fs=1)
#     ft = stft_norm(args...; fs)
#     plot(ft)
# end

function process_audio_tonal(aufname::String, res_dir)
    filetype = splitext(aufname)[2] |> lowercase
    if ".wav" == filetype 
        data, fs = wavread(aufname, format="native")
    elseif ".mat" == filetype
        d = load(aufname)
        data = d["data"]
        haskey(d, "fs") ? fs=d["fs"] : fs=500_000
    else
        data, fs = load(aufname)
    end

    process_audio_tonal( (aufname,data,fs), res_dir)
end

#     process_audioVideo_tonal1(, vidfname, res_dir)

function process_audioVideo_tonal1(aufname_data_fs::Tuple, vidfname, res_dir;
    ref_channel = ref_channel, thresh_tonal = threshold_tonal, 
    nfft_inS = nfft_inS, band_pass = band_pass, 
    butterworth_size = butterworth_size,c = c, winlen=winlen)

    aufname, data, fs = aufname_data_fs
    # ref_channel = 3
    # thresh_tonal = -20
    # nfft_inS = 0.01
    # band_pass = [5000, 14000] # [5000, 25000] # 
    # butterworth_size = 4
    # c = 1540

    # winlen=0.01s # tdoa estimation window length

    # data, fs = wavread(aufname, format="native")
    @info ("Duration: " * string(size(data,1)/fs) *"seconds")

    nfft = nextfastfft(nfft_inS*fs)
    freq_filt = 2:Int(nfft/2+1)
    if !isnothing(band_pass)
        freq_filt = Int(round(band_pass[1]/fs*nfft)):Int(round(band_pass[2]/fs*nfft)) #50:140
    end
    filter_weight = digitalfilter(Bandpass(band_pass[1],band_pass[2], fs=fs), Butterworth(butterworth_size)) #Bandpass(5000,14000, fs=fs)

    
    specs = mapslices( x -> spectrogram(x, nfft, 0; fs=fs), data[:,1:3], dims=1)
    mag_ft = broadcast( x -> pow2db.(x.power) ,specs)

    time_index = first(specs).time#[window]
    fft_max = maximum(mag_ft[ref_channel][freq_filt,:]; dims=1)'[:]

    tonal_indices = findall( >(thresh_tonal), fft_max)
    tonal_indices = tonal_indices[1:end-1] #FIXME #tempsolution

    @info "Num of Whistles Detected: " * string(length(tonal_indices))
    detected_tonal_inS = time_index[tonal_indices]
    mkpath(res_dir)
    audacity_label(detected_tonal_inS, joinpath(res_dir, splitext(aufname)[1]*"___tonal_t"*string(thresh_tonal)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".txt" |> basename) )

    data_filt = mapslices( x -> filtfilt( filter_weight, x), data[:,1:end-1], dims=1) #data, assumise last channel is synch channel and not used for processing
    d_filt=signal(data_filt, fs)

    tdoa = twindow2tdoa.( (detected_tonal_inS)s; fs=fs, d_filt=d_filt, winlen=winlen)
    pxs = twindow2dird.( (detected_tonal_inS)s; fs=fs, d_filt=d_filt, winlen=winlen)

    ## create return variables:
    ang = vcat(pxs...) #.|> deg2rad
    p_pixels = angle2px(ang)

    pind_vidframes = round.(Int, detected_tonal_inS * get_fps(vidfname)) .+ 1    

    # writedlm( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh_tonal)*"tonal_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".csv"), ["p_pixel" "px" "py"], ',')
    open( joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh_tonal)*"tonal_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".csv"), "w") do io
        writedlm(io, ["p_pixel" "px" "py"], ',')
        writedlm(io, [pind_vidframes p_pixels], ',')
    end
    
    return pind_vidframes, p_pixels, thresh_tonal, band_pass, vcat(pxs...) .|> deg2rad, tdoa, nothing, winlen, tonal_indices, nothing, detected_tonal_inS, nothing, fft_max[tonal_indices], ref_channel, c, rx_vect, fs
end

function quietest_segment(specs)
    
end

function detect_tonal(aufname::String, args...; kwargs...)
    data, fs, nbits, opt, timestamp = readAudio(aufname)
    detect_tonal((aufname, data, fs), args...; kwargs...)
end

function detect_tonal(aufname_data_fs::Tuple, res_dir=nothing;
    ref_channel = ref_channel, thresh_tonal = threshold_tonal, 
    nfft_inS = nfft_inS, band_pass = tonal_band_pass, window=window_tonal,
    freq_maxbandwidth = freq_maxbandwidth, freq_width_db=freq_width_db,
    percent_quiet=percent_quiet,
    stats_func=maximum)

    aufname, data, fs = aufname_data_fs
    # ref_channel = 3
    # thresh_tonal = -20
    # nfft_inS = 0.01
    # band_pass = [5000, 14000] # [5000, 25000] # 
    # butterworth_size = 4
    # c = 1540

    # winlen=0.01s # tdoa estimation window length

    # data, fs = wavread(aufname, format="native")
    @info ("Duration: " * string(size(data,1)/fs) *"seconds")

    nfft = nextfastfft(nfft_inS*fs)
    freq_filt = 2:Int(nfft/2+1)
    if !isnothing(band_pass)
        freq_filt = Int(round(band_pass[1]/fs*nfft)):Int(round(band_pass[2]/fs*nfft)) #50:140
    end
    # filter_weight = digitalfilter(Bandpass(band_pass[1],band_pass[2], fs=fs), Butterworth(butterworth_size)) #Bandpass(5000,14000, fs=fs)

    
    specs = mapslices( x -> spectrogram(x, nfft, 0; fs=fs), data, dims=1)
    mag_ft = broadcast( x -> pow2db.(x.power) ,specs)

    mag = mag_ft[ref_channel][freq_filt,:]
    # magsum = zeros(size(mag))
    # for i in 1:length(mag_ft)
    #     magsum = magsum .+ mag_ft[i][freq_filt,:]
    # end
    # magsum /= length(mag_ft)

    mf = cat(mag_ft..., dims=3)
    # mag_mean = dropdims( mean(mf[freq_filt,:,:],dims=3), dims=3)
    mag_max = dropdims( maximum(mf[freq_filt,:,:],dims=3), dims=3)
    @debug mag_max|>size
    #
    quietest = nothing
    if percent_quiet > 0
        @debug percent_quiet
        # mag_max = mag_max .- mean(sort(mag_max, dims=2)[:,1:round(Int,size(mag_max,2) * percent_quiet)], dims=2)
        num_bins = round(Int,size(mag_max,2) * percent_quiet)
        if num_bins < 1
            @warn "!! Using ONE time bin only as specified quiet segment percentage is smaller than one time bin."
            num_bins = size(mag_max,2)
        end
        quietest = maximum(mag_max[:,sortperm(sum(mag_max, dims=1)[:])[1:num_bins]], dims=2)
        mag_max = mag_max .- quietest
    end
    freqss = first(specs).freq[freq_filt]
    time_index = first(specs).time#[window]
    # fft_max = maximum(mag; dims=1)'[:]
    fft_max = stats_func(mag_max, dims=1)[:]

    #~ filtering for small bandwidth timeframe only
    maxi_t = maximum(mag, dims=1)
    small_bands = (sum((maxi_t .- mag ) .< freq_width_db, dims=1) .*freqss[2] ) .< freq_maxbandwidth
    small_bands = [ele==0 ? -Inf : 0 for ele in small_bands]
    # @debug size(fft_max)
    # @debug size(small_bands)
    if isnothing(thresh_tonal)
        thresh_tonal = median(fft_max) + 8
        @info "__auto thresholding: " * string(thresh_tonal)
    end
    # @debug size(small_bands)
    tonal_indices = findall( >(thresh_tonal), fft_max .+ small_bands')
    tonal_indices = tonal_indices[1:end-1] #FIXME #tempsolution

    num_detection = length(tonal_indices)
    @info "Num of Tonals Detected: " * string(num_detection)
    pind_good_inS = time_index[tonal_indices]

    outfname = nothing
    if !isnothing(res_dir)
        mkpath(res_dir)
        # @debug "audacity label: "*joinpath(res_dir, splitext(aufname)[1]*"___tonal_t"*string(thresh_tonal)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".txt" |> basename)
        outfname = splitext(aufname)[1]*"___tonal_t"*string(thresh_tonal)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2]) *"_bw"*string(freq_maxbandwidth) *"_pquiet"*string(percent_quiet) |> basename
        audacity_label(pind_good_inS, joinpath(res_dir, outfname * ".txt") )
    end
    
    pind_good = round.(Int, pind_good_inS.*fs)
    return (;fft_max, time_index , tonal_indices, pind_good_inS, freqss, num_detection, percent_quiet, freq_maxbandwidth, threshold=thresh_tonal, ppeak=fft_max[tonal_indices], quietest, aufname, fs, outfname, pind_good, len_data=size(data,1), band_pass)#, mag_max, specs) #mag, mag_ft, mag_mean, mag_max, 
    # return (;pind=, ppeak, pind_good, pind_good_inS, threshold_indices)
end

conv_onesided_centreMID(u,v) = conv(u, v)[ (length(v)+1:end-(length(v)-1)) .- length(v)÷2]

function combine_detections_conv(data, res_new, nfft=1280, infl_len=30;
    winfunc=gaussian, σ=0.1,
    threshold=0.01, thresh_continue=0.00001,
    res_dir=nothing)
# σ=0.1; nfft=1280; infl_len=30;
new_dat = zeros(size(data,1),1)
new_dat[res_new.pind_good] .= 1

prob = conv_onesided_centreMID(new_dat, gaussian(nfft*infl_len,σ));
# prob = mfilter(gaussian(nfft*infl_len,σ), new_dat)
# prob2 = conv_onesided(new_dat, gaussian(nfft*infl_len,σ))

#~ look for sets
# click_accum = Int[]
# click_in_trains = Int[]
train_switch = false
train_count = 0
train_start = Int[]
train_end = Int[]
# train_start_ind = Int[]
for (i,val) in enumerate(prob)
    if val > threshold
        # push!(click_accum, i)
        # push!(click_in_trains, res.pind_good[i])
        if !train_switch
            train_count += 1
            push!(train_start, i)
            # push!(train_start_ind, length(click_in_trains))
        end
        train_switch = true
    elseif train_switch
        if val < thresh_continue
            push!(train_end, i)
            train_switch = false
        end
    end
end
if length(train_start) >  length(train_end)
    push!(train_end, res_new.len_data)
end

@info "Number of Tonal Segments: "* string(length(train_start))

pind_good_inS = train_start ./ res_new.fs
# @debug pind_good_inS
# @debug res_new.time_index
ppeak = map(i -> 
    res_new.fft_max[ (findlast( <(pind_good_inS[i]), res_new.time_index) |> x -> isnothing(x) ? 1 : x): (findfirst( >(train_end[i]/res_new.fs), res_new.time_index) |> x -> isnothing(x) ? length(res_new.fft_max) : x  )] |> median
, eachindex(train_start))
# @debug "ppeak: "* string(size(ppeak))

if !isnothing(res_dir)
    audacity_label([train_start train_end] ./ res_new.fs, joinpath(res_dir, res_new.outfname *"__len"* string(nfft*infl_len) *"_sigma"*string(σ)*  "_segment-only.txt" |> basename))
end

return (;pind_good_inS, 
pind_good = train_start,
ppeak,
pind = res_new.pind_good,
train_start, train_end,
num_detection=length(train_start), num_click_in_trains=length(res_new.pind_good),
)

# plot!( (1:nfft*infl_len) ./ res_new.fs, gaussian(nfft*infl_len,σ))
# plot!( ((1:nfft*infl_len) .+ nfft*i) ./ res_new.fs, gaussian(nfft*infl_len,σ))
# plot!( ((1:nfft*infl_len) .+ nfft*i) ./ fs, gaussian(round(Int,nfft_inS*fs)*30,σ))
# wavwrite(prob, fs, joinpath(res_dir,"prob_s0.1.wav"))
end

# spec = res_tonal.specs[3].power
# # function pentropy(spec; removedc_flag=true)
#     if removedc_flag
#         window = 2:size(spec,1)
#     else
#         window = 1:size(spec,1)
#     end
#     denom = sum(@view spec[window,:])

#     white_noise_norm = log2(length(spec))
#     S = Array{Float64}(undef, size(spec,2))
#     for t in axes(spec,2)
#         P = sum(@view spec[window,t]) / denom
#         S[t] = - sum( P*log2(P) ) / white_noise_norm
#     end

#     plot(S)
#     plot(res_tonal.time_index, S)


function detect_tonal2(aufname_data_fs::Tuple, res_dir=nothing;
    ref_channel = ref_channel, thresh_tonal = threshold_tonal, 
    nfft_inS = nfft_inS, band_pass = tonal_band_pass, window=window_tonal,
    freq_maxbandwidth = freq_maxbandwidth, freq_width_db=freq_width_db,
    percent_quiet=percent_quiet,
    stats_func=maximum)

    aufname, data, fs = aufname_data_fs
    # ref_channel = 3
    # thresh_tonal = -20
    # nfft_inS = 0.01
    # band_pass = [5000, 14000] # [5000, 25000] # 
    # butterworth_size = 4
    # c = 1540

    # winlen=0.01s # tdoa estimation window length

    # data, fs = wavread(aufname, format="native")
    @info ("Duration: " * string(size(data,1)/fs) *"seconds")

    nfft = nextfastfft(nfft_inS*fs)
    freq_filt = 2:Int(nfft/2+1)
    if !isnothing(band_pass)
        freq_filt = Int(round(band_pass[1]/fs*nfft)):Int(round(band_pass[2]/fs*nfft)) #50:140
    end
    # filter_weight = digitalfilter(Bandpass(band_pass[1],band_pass[2], fs=fs), Butterworth(butterworth_size)) #Bandpass(5000,14000, fs=fs)

    
    specs = mapslices( x -> spectrogram(x, nfft, 0; fs=fs), data, dims=1)
    mag_ft = broadcast( x -> sqrt.(x.power) ,specs)

    mag = mag_ft[ref_channel][freq_filt,:]
    # magsum = zeros(size(mag))
    # for i in 1:length(mag_ft)
    #     magsum = magsum .+ mag_ft[i][freq_filt,:]
    # end
    # magsum /= length(mag_ft)

    mf = cat(mag_ft..., dims=3)
    # mag_mean = dropdims( mean(mf[freq_filt,:,:],dims=3), dims=3)
    mag_max = dropdims( maximum(mf[freq_filt,:,:],dims=3), dims=3)
    @debug mag_max|>size
    #
    # quietest = nothing
    # if percent_quiet > 0
    #     @debug percent_quiet
    #     # mag_max = mag_max .- mean(sort(mag_max, dims=2)[:,1:round(Int,size(mag_max,2) * percent_quiet)], dims=2)

    #     num_bins = round(Int,size(mag_max,2) * percent_quiet)
    #     if num_bins < 1
    #         @warn "!! Using ONE time bin only as specified quiet segment percentage is smaller than one time bin."
    #         num_bins = size(mag_max,2)
    #     end
    #     quietest = maximum(mag_max[:,sortperm(sum(mag_max, dims=1)[:])[1:num_bins]], dims=2)
        
    #     mag_max = mag_max .- quietest
    # end
    quietest = minimum(mag_max; dims=2)
    mag_max = mag_max .- quietest
    mag_max = map( x -> x > 0 ? x : 0, mag_max)
    mag_max = 2 .* pow2db.(mag_max)
    mag = 2 .* pow2db.(mag)
    
    freqss = first(specs).freq[freq_filt]
    time_index = first(specs).time#[window]
    # fft_max = maximum(mag; dims=1)'[:]
    fft_max = stats_func(mag_max, dims=1)[:]

    #~ filtering for small bandwidth timeframe only
    maxi_t = maximum(mag, dims=1)
    small_bands = (sum((maxi_t .- mag ) .< freq_width_db, dims=1) .*freqss[2] ) .< freq_maxbandwidth
    small_bands = [ele==0 ? -Inf : 0 for ele in small_bands]
    # @debug size(fft_max)
    # @debug size(small_bands)
    if isnothing(thresh_tonal)
        thresh_tonal = median(fft_max) + 8
        @info "__auto thresholding: " * string(thresh_tonal)
    end
    # @debug size(small_bands)
    tonal_indices = findall( >(thresh_tonal), fft_max .+ small_bands')
    tonal_indices = tonal_indices[1:end-1] #FIXME #tempsolution

    num_detection = length(tonal_indices)
    @info "Num of Tonals Detected: " * string(num_detection)
    pind_good_inS = time_index[tonal_indices]

    outfname = nothing
    if !isnothing(res_dir)
        mkpath(res_dir)
        # @debug "audacity label: "*joinpath(res_dir, splitext(aufname)[1]*"___tonal_t"*string(thresh_tonal)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2])*".txt" |> basename)
        outfname = splitext(aufname)[1]*"___tonal_t"*string(thresh_tonal)*"_bp"*string(band_pass[1])*"_"*string(band_pass[2]) *"_bw"*string(freq_maxbandwidth) *"_pquiet"*string(percent_quiet) |> basename
        audacity_label(pind_good_inS, joinpath(res_dir, outfname * ".txt") )
    end
    
    pind_good = round.(Int, pind_good_inS.*fs)
    return (;fft_max, time_index , tonal_indices, pind_good_inS, freqss, num_detection, percent_quiet, freq_maxbandwidth, threshold=thresh_tonal, ppeak=fft_max[tonal_indices], quietest, aufname, fs, outfname, pind_good, len_data=size(data,1), small_bands, mag_max, specs) #mag, mag_ft, mag_mean, mag_max, 
    # return (;pind=, ppeak, pind_good, pind_good_inS, threshold_indices)
end