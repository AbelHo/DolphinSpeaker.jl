using Base.Threads
using Random
using DSP
using Plots
using Statistics
using Test
include("../src/audio.jl")

const FS = 8000
const DURATION = 1.5
const NFFT = 512
const HOP = 128
const SAMPLE_AUDIO_FILE = "/media/spin/anas2/data/marecet/sample/7003.240514151859_465s.flac"

"""Return the most likely ridge (one bin per frame) via a thread-safe Viterbi sweep."""
function viterbi_ridge(S::AbstractMatrix{<:Real}; lambda_smooth::Real = 10.0)
    n_f, n_t = size(S)
    V = zeros(eltype(S), n_f, n_t)
    B = zeros(Int, n_f, n_t)
    V[:, 1] .= S[:, 1]

    freq_bins = collect(1:n_f)
    for t in 2:n_t
        prev = @view V[:, t - 1]
        @threads for f in 1:n_f
            cost = prev .- lambda_smooth .* abs.(freq_bins .- f)
            _, idx = findmax(cost)
            B[f, t] = idx
            V[f, t] = S[f, t] + cost[idx]
        end
    end

    ridge = Vector{Int}(undef, n_t)
    ridge[end] = argmax(V[:, end])
    for t in (n_t):-1:2
        ridge[t - 1] = B[ridge[t], t]
    end

    return ridge
end

function viterbi_multi_ridge(S::AbstractMatrix{<:Real};
                             num_tracks::Int = 3,
                             lambda_smooth::Real = 12.0,
                             suppression_width::Int = 4)
    working = copy(S)
    ridges = Vector{Vector{Int}}()
    n_f = size(S, 1)
    for _ in 1:num_tracks
        ridge = viterbi_ridge(working; lambda_smooth=lambda_smooth)
        push!(ridges, ridge)
        for (frame, bin) in enumerate(ridge)
            lo = max(bin - suppression_width, 1)
            hi = min(bin + suppression_width, n_f)
            working[lo:hi, frame] .= 0.0
        end
    end
    return ridges
end

function find_best_active_segment(mask::AbstractVector{Bool}, min_length::Int)
    best_start = 0
    best_len = 0
    cur_start = 0
    cur_len = 0
    for (i, on) in enumerate(mask)
        if on
            cur_len += 1
            cur_start == 0 && (cur_start = i)
        elseif cur_len > 0
            if cur_len > best_len
                best_len = cur_len
                best_start = cur_start
            end
            cur_len = 0
            cur_start = 0
        end
    end
    if cur_len > best_len
        best_len = cur_len
        best_start = cur_start
    end
    if best_len >= min_length
        return best_start, best_start + best_len - 1
    elseif best_len > 0
        return best_start, best_start + best_len - 1
    end
    first_one = findfirst(mask)
    if first_one !== nothing
        return first_one, first_one
    end
    return 1, length(mask)
end

function viterbi_ridge_with_boundaries(S::AbstractMatrix{<:Real};
                                        lambda_smooth::Real = 10.0,
                                        energy_threshold_ratio::Real = 0.25,
                                        min_duration_frames::Int = 6)
    ridge = viterbi_ridge(S; lambda_smooth=lambda_smooth)
    energies = [S[ridge[t], t] for t in eachindex(ridge)]
    max_energy = maximum(energies)
    threshold = energy_threshold_ratio * max_energy
    mask = energies .>= threshold
    start_frame, stop_frame = find_best_active_segment(mask, min_duration_frames)
    return ActiveRidge(ridge, energies, start_frame, stop_frame)
end

"""Holds the STFT, time/freq axes, and ridge metadata."""
struct RidgeResult
    times::Vector{Float64}
    freq_bins::Vector{Float64}
    magnitude::Matrix{Float64}
    ridge_idx::Vector{Int}
    estimated_freq::Vector{Float64}
end

struct MultiRidgeResult
    times::Vector{Float64}
    freq_bins::Vector{Float64}
    magnitude::Matrix{Float64}
    ridges::Vector{Vector{Int}}
    estimated_freqs::Matrix{Float64}
end

struct ActiveRidge
    ridge::Vector{Int}
    energies::Vector{Float64}
    start_frame::Int
    stop_frame::Int
end

"""Compute the ridge from a raw signal and return the analysis results."""
function viterbi_ridge_signal(signal::AbstractVector{<:Real};
                              fs::Real = FS,
                              nfft::Int = NFFT,
                              hop::Int = HOP,
                              lambda_smooth::Real = 12.0,
                              window::Union{Nothing, AbstractVector{<:Real}} = nothing)
    noverlap = nfft - hop
    window = window === nothing ? DSP.hanning(nfft) : window
    S = stft(signal, nfft, noverlap; window=window, nfft=nfft, fs=fs)
    magnitude = abs.(S)
    n_frames = size(magnitude, 2)
    times = ((0:n_frames - 1) .* hop) ./ fs
    freq_bins = (0:size(magnitude, 1) .- 1) .* (fs / nfft)
    ridge_idx = viterbi_ridge(magnitude; lambda_smooth=lambda_smooth)
    estimated_freq = freq_bins[ridge_idx]
    return RidgeResult(times, freq_bins, magnitude, ridge_idx, estimated_freq)
end

function viterbi_multi_ridge_signal(signal::AbstractVector{<:Real};
                                     fs::Real = FS,
                                     nfft::Int = NFFT,
                                     hop::Int = HOP,
                                     lambda_smooth::Real = 12.0,
                                     num_tracks::Int = 3,
                                     suppression_width::Int = 4,
                                     window::Union{Nothing, AbstractVector{<:Real}} = nothing)
    noverlap = nfft - hop
    window = window === nothing ? DSP.hanning(nfft) : window
    S = stft(signal, nfft, noverlap; window=window, nfft=nfft, fs=fs)
    magnitude = abs.(S)
    n_frames = size(magnitude, 2)
    times = ((0:n_frames - 1) .* hop) ./ fs
    freq_bins = (0:size(magnitude, 1) .- 1) .* (fs / nfft)
    ridges = viterbi_multi_ridge(magnitude;
                                 num_tracks=num_tracks,
                                 lambda_smooth=lambda_smooth,
                                 suppression_width=suppression_width)
    estimated_freqs = reduce(hcat, freq_bins[ridge] for ridge in ridges)
    return MultiRidgeResult(times, freq_bins, magnitude, ridges, estimated_freqs)
end

"""Draw the STFT + ridge overlay using consistent styling."""
function plot_ridge_result(result::RidgeResult;
                           title::String = "STFT with detected ridge",
                           mask_color = :red,
                           mask_alpha = 0.35,
                           ridge_color = :white,
                           track_only::Bool = false,
                           active_segment::Union{Nothing, ActiveRidge} = nothing)
    if track_only
        plt = plot(result.times, result.estimated_freq;
               color = ridge_color,
               linewidth = 2,
               xlabel = "Time (s)",
               ylabel = "Frequency (Hz)",
               title = title,
               legend = false)
        return plt
    end
    mask = zeros(Float64, size(result.magnitude))
    peak = maximum(result.magnitude)
    for (frame, row) in enumerate(result.ridge_idx)
        mask[row, frame] = peak
    end

    plt = heatmap(
        result.times,
        result.freq_bins,
        20 .* log10.(result.magnitude .+ eps()),
        color = :viridis,
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
        title = title,
        legend = false,
        colorbar = true,
        label = ""
    )
    heatmap!(plt, result.times, result.freq_bins, mask;
             color = mask_color,
             alpha = 0.5,#mask_alpha,
             label = "ridge mask",
             colorbar = false)
    plot!(plt, result.times, result.estimated_freq;
          color = ridge_color,
          linewidth = 2,
          label = "ridge track",
          alpha = mask_alpha)
    if active_segment !== nothing
        start = result.times[clamp(active_segment.start_frame, 1, end)]
        stop = result.times[clamp(active_segment.stop_frame, 1, end)]
        vline!(plt, [start]; color = :white, linestyle = :dash, label = "start")
        vline!(plt, [stop]; color = :orange, linestyle = :dash, label = "stop")
    end
    return plt
end

function plot_multi_ridge_result(result::MultiRidgeResult;
                                  title::String = "STFT multi-track ridge",
                                  mask_color = :red,
                                  mask_alpha = 0.35,
                                  track_only::Bool = false)
    if track_only
        colors = [:white, :yellow, :cyan, :lime, :magenta]
        plt = plot(; xlabel = "Time (s)", ylabel = "Frequency (Hz)", title = title)
        for (idx, freq_track) in enumerate(eachcol(result.estimated_freqs))
            plot!(plt, result.times, freq_track;
                  color = colors[mod1(idx, length(colors))],
                  linewidth = 2,
                  linestyle = :solid,
                  alpha = 0.9,
                  label = "track $idx")
        end
        return plt
    end
    mask = zeros(Float64, size(result.magnitude))
    peak = maximum(result.magnitude)
    for ridge in result.ridges
        for (frame, row) in enumerate(ridge)
            mask[row, frame] = peak
        end
    end

    plt = heatmap(
        result.times,
        result.freq_bins,
        20 .* log10.(result.magnitude .+ eps()),
        color = :viridis,
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
        title = title,
        legend = false,
        colorbar = true,
        label = ""
    )
    heatmap!(plt, result.times, result.freq_bins, mask;
             color = mask_color,
             alpha = mask_alpha,
             label = "ridge mask",
             colorbar = false)
    colors = [:white, :yellow, :cyan, :lime, :magenta]
    for (idx, freq_track) in enumerate(eachcol(result.estimated_freqs))
        plot!(plt, result.times, freq_track;
              color = colors[mod1(idx, length(colors))],
              linewidth = 2,
              linestyle = :solid,
              alpha = 0.3,
              label = "track $idx")
    end
    return plt
end

function process_signal_multi_ridge(signal::AbstractVector{<:Real}, fs::Real;
                                    kwargs...)
    return viterbi_multi_ridge_signal(signal; fs=fs, kwargs...)
end

function process_audio_file_multi_ridge(aufname::AbstractString;
                                        channel::Int = 1,
                                        kwargs...)
    data, fs, _, _, _ = readAudio(aufname)
    signal = ndims(data) == 1 ? vec(data) : vec(mean(data, dims=2))
    return process_signal_multi_ridge(signal, fs; kwargs...)
end

function simulate_ridge_signal()
    n_samples = Int(round(FS * DURATION))
    times = range(0, step=1 / FS, length=n_samples)
    freq_fn(t) = 500 + 120 * sin(2π * 3 * t / DURATION)
    inst_freq = freq_fn.(times)
    phase = 2π .* cumsum(inst_freq) ./ FS
    signal = sin.(phase)
    Random.seed!(1234)
    signal .+= 0.08 .* randn(length(signal))
    return signal, freq_fn
end

function main(track_only::Bool=false)
    signal, freq_fn = simulate_ridge_signal()
    result = viterbi_ridge_signal(signal; lambda_smooth=12.0)
    active = viterbi_ridge_with_boundaries(result.magnitude;
                                           lambda_smooth=12.0,
                                           min_duration_frames=8)
    p = plot_ridge_result(result;
                          title = "Simulated STFT with detected ridge",
                          mask_alpha = 0.35,
                          ridge_color = :white,
                          track_only=track_only,
                          active_segment=active)
    filepath = track_only ? "results/viterbi_ridge_tracks.png" : "results/viterbi_ridge.png"
    savefig(filepath)
    println("Saved figure to $filepath")
    println("Active ridge starts at frame $(active.start_frame), stops at $(active.stop_frame)")
    target_freq = freq_fn.(result.times)
    @test mean(abs.(result.estimated_freq .- target_freq)) < 80
    println("Viterbi ridge tracks the simulated ridge within acceptable tolerance.")
    return p
end

function make_fm_signal(times, center, depth, rate)
    freq_curve = center .+ depth .* sin.(2π .* rate .* times)
    phase = 2π .* cumsum(freq_curve) ./ FS
    return sin.(phase), freq_curve
end

function match_multi_track_error(actual::Matrix{<:Real}, estimated::Matrix{<:Real})
    n_frames, n_tracks = size(actual)
    total = 0.0
    for t in 1:n_frames
        for j in 1:n_tracks
            total += minimum(abs.(estimated[t, :] .- actual[t, j]))
        end
    end
    return total / (n_frames * n_tracks)
end

function main_viterbi_ridge_signal(track_only::Bool=false)
    n_samples = Int(round(FS * DURATION))
    full_times = range(0, step=1 / FS, length=n_samples)
    component_settings = (
        (650.0,  80.0, 1.3),
        (1000.0, 140.0, 2.7),
        (1450.0, 120.0, 2.0),
        (2000.0, 100.0, 2.0),
    )
    signal = zeros(Float64, length(full_times))
    for (center, depth, rate) in component_settings
        comp_signal, _ = make_fm_signal(full_times, center, depth, rate)
        signal .+= comp_signal
    end
    Random.seed!(2025)
    signal .+= 0.04 .* randn(length(signal))

    result = viterbi_multi_ridge_signal(signal;
                                        lambda_smooth=14.0,
                                        num_tracks=length(component_settings),
                                        suppression_width=5)
    p = plot_multi_ridge_result(result;
                                title = "Composite FM STFT with detected ridges",
                                mask_alpha = 0.25,
                                track_only=track_only)

    actual_freqs = reduce(hcat, (center .+ depth .* sin.(2π .* rate .* result.times))
                               for (center, depth, rate) in component_settings)
    # for j in 1:size(actual_freqs, 2)
    #     plot!(p, result.times, actual_freqs[:, j];
    #           color = :gray,
    #           linestyle = :dash,
    #           alpha = 0.5,
    #           linewidth = 1,
    #           label = "component $j")
    # end

    filepath = track_only ? "results/viterbi_ridge_composite_tracks.png" : "results/viterbi_ridge_composite.png"
    savefig(filepath)
    println("Saved composite figure to $filepath")
    avg_error = match_multi_track_error(actual_freqs, result.estimated_freqs)
    @test avg_error < 120
    println("Multi-track ridge average error $(round(avg_error, digits=1)) Hz")
    return p
end

function main_audio_file_ridge(aufilename::AbstractString = SAMPLE_AUDIO_FILE;
                               num_tracks::Int = 4,
                               track_only::Bool = false)
    result = process_audio_file_multi_ridge(aufilename;
                                           lambda_smooth=16.0,
                                           num_tracks=num_tracks,
                                           suppression_width=8)
    p = plot_multi_ridge_result(result;
                                title = "Audio file multi-track ridge",
                                mask_alpha = 0.3,
                                track_only=track_only)
    filepath = track_only ? "results/viterbi_ridge_audio_tracks.png" : "results/viterbi_ridge_audio.png"
    savefig(filepath)
    println("Saved audio file figure to $filepath")
    return p
end

if abspath(PROGRAM_FILE) == @__FILE__
    # create both full spectrogram+overlay and track-only variations for quick inspection
    main(false)
    main(true)
    main_viterbi_ridge_signal(false)
    main_viterbi_ridge_signal(true)
    main_audio_file_ridge(SAMPLE_AUDIO_FILE, num_tracks=4, track_only=false)
    main_audio_file_ridge(SAMPLE_AUDIO_FILE, num_tracks=4, track_only=true)
end
