# Run directly: julia --project=. test/test_hydromoth_dugong_preset.jl
using Test
using DolphinSpeaker

const DS = DolphinSpeaker

@testset "HydroMoth dugong tonal preset" begin
    # Reject bad input before mutating configuration.
    @test_throws ArgumentError DS.set_device__hydromoth_dugong(sensitivity=:invalid)

    @test DS.set_device__hydromoth_dugong() === nothing
    @test DS.tonal_band_pass == [800, 16_000]
    @test DS.nfft_inS == 0.02
    @test DS.window_tonal == 0.02
    @test DS.ref_channel == 1
    @test DS.threshold_tonal === nothing
    @test DS.tonal_threshold_median == 5
    @test DS.percent_quiet == 0.01
    @test isinf(DS.freq_maxbandwidth)
    @test DS.DEFAULT_fname2timestamp_func === DS.fname2dt_date_time

    @test DS.set_device__hydromoth_dugong(sensitivity=:balanced) === nothing
    @test DS.tonal_band_pass == [800, 16_000]
    @test DS.nfft_inS == 0.02
    @test DS.tonal_threshold_median == 8
    @test DS.freq_width_db == 3
    @test DS.freq_maxbandwidth == 15_000

    # Repeated calls should not leave an earlier sensitivity in effect.
    DS.set_device__hydromoth_dugong()
    @test DS.tonal_threshold_median == 5
    @test isinf(DS.freq_maxbandwidth)
end
