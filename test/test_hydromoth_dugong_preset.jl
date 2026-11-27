# Run directly: julia --project=. test/test_hydromoth_dugong_preset.jl
using Test
using DolphinSpeaker

const DS = DolphinSpeaker

@testset "HydroMoth dugong chirp/trill presets" begin
    # Invalid options must not mutate the current device configuration.
    DS.set_device__hydromoth_dugong()
    prior_band = copy(DS.tonal_band_pass)
    @test_throws ArgumentError DS.set_device__hydromoth_dugong(call_type=:whistle)
    @test_throws ArgumentError DS.set_device__hydromoth_dugong(sensitivity=:invalid)
    @test DS.tonal_band_pass == prior_band

    @test DS.set_device__hydromoth_dugong() === nothing
    @test DS.tonal_band_pass == [1_200, 10_500]
    @test DS.nfft_inS == 0.010
    @test DS.window_tonal == DS.nfft_inS
    @test DS.ref_channel == 1
    @test DS.threshold_tonal === nothing
    @test DS.tonal_threshold_median == 7
    @test DS.percent_quiet == 0.05
    @test DS.freq_width_db == 6
    @test DS.freq_maxbandwidth == 3_000
    @test DS.DEFAULT_fname2timestamp_func === DS.fname2dt_date_time

    @test DS.set_device__hydromoth_dugong(call_type=:chirp) === nothing
    @test DS.tonal_band_pass == [1_200, 10_500]
    @test DS.nfft_inS == 0.010
    @test DS.freq_maxbandwidth == 2_500

    @test DS.set_device__hydromoth_dugong(call_type=:trill) === nothing
    @test DS.tonal_band_pass == [1_200, 6_500]
    @test DS.nfft_inS == 0.020
    @test DS.freq_maxbandwidth == 2_000

    @test DS.set_device__hydromoth_dugong(call_type=:trill, sensitivity=:high_recall) === nothing
    @test DS.tonal_threshold_median == 4
    @test isinf(DS.freq_maxbandwidth)

    # Switching profiles must reset settings from earlier calls.
    DS.set_device__hydromoth_dugong()
    @test DS.tonal_band_pass == [1_200, 10_500]
    @test DS.nfft_inS == 0.010
    @test DS.freq_maxbandwidth == 3_000
    @test DS.tonal_threshold_median == 7
end

@testset "Tonal detector chirp/trill signal smoke test" begin
    # Deterministic, clean synthetic tones check the detector mechanics.
    # They do NOT establish performance on real HydroMoth dugong WAVs.
    fs = 48_000
    t = (0:(2 * fs - 1)) ./ fs
    signal = 0.0008 .* sin.(2pi .* 780 .* t) .+
             0.0003 .* sin.(2pi .* 1_090 .* t)

    # Chirp: 3 -> 5 kHz over 120 ms, beginning at 0.30 s.
    chirp_start, chirp_duration = 0.30, 0.12
    chirp_indices = (round(Int, chirp_start * fs) + 1):round(Int, (chirp_start + chirp_duration) * fs)
    chirp_t = t[chirp_indices] .- chirp_start
    signal[chirp_indices] .+= 0.2 .* sin.(2pi .* (
        3_000 .* chirp_t .+ 0.5 .* (2_000 / chirp_duration) .* chirp_t .^ 2
    ))

    # Trill: frequency-modulated 2.9--3.7 kHz over 650 ms.
    trill_start, trill_duration = 1.0, 0.65
    trill_indices = (round(Int, trill_start * fs) + 1):round(Int, (trill_start + trill_duration) * fs)
    trill_t = t[trill_indices] .- trill_start
    signal[trill_indices] .+= 0.2 .* sin.(2pi .* 3_300 .* trill_t .+
        (400 / 16) .* sin.(2pi .* 16 .* trill_t))

    data = reshape(signal, :, 1)
    DS.set_device__hydromoth_dugong()
    res = DS.detect_tonal(("synthetic_chirp_trill", data, fs))
    @test size(res.output_res) == size(res.fft_max)
    @test length(res.time_index) == length(res.output_res)
    @test all(p -> 1 <= p <= length(signal), res.pind_good)
    @test minimum(res.freqss) >= 1_200
    @test maximum(res.freqss) <= 10_500
    @test any(ti -> 0.29 <= ti <= 0.43, res.pind_good_inS)
    @test any(ti -> 0.99 <= ti <= 1.66, res.pind_good_inS)

    # A threshold above every possible finite score should yield an empty
    # detection vector without slicing away or inventing frames.
    empty_res = DS.detect_tonal(("synthetic_chirp_trill", data, fs); thresh_tonal=Inf)
    @test isempty(empty_res.tonal_indices)

    DS.set_device__hydromoth_dugong(call_type=:trill, sensitivity=:high_recall)
    trill_res = DS.detect_tonal(("synthetic_chirp_trill", data, fs))
    @test any(ti -> 0.99 <= ti <= 1.66, trill_res.pind_good_inS)
end
