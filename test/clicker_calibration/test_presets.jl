using Test
if !isdefined(@__MODULE__,:ClickerCalibration)
    include(joinpath(@__DIR__,"../../src/clicker_calibration/ClickerCalibration.jl"))
end

module CalibratedLegacyTest
using LinearAlgebra, Statistics, Dates
# The preset's timestamp callback is irrelevant to these array-only controls.
fname2dt_zoomf6(path)=nothing
include(joinpath(@__DIR__,"../../src/localization.jl"))
end

@testset "calibrated presets apply delay/camera corrections exactly once" begin
    mktempdir() do dir
        L=CalibratedLegacyTest; C=ClickerCalibration
        qa=C._pf_initial(2704,1520,"centered_pinhole")
        qa[1]=-0.3;qa[2:4]=[0.025,-0.035,0.018]./pi;qa[7:8]=[0.03,-0.02]
        ca=C._pf_calibration(qa,2704,1520,"centered_pinhole";device="A")
        cb=deepcopy(ca);cb["device"]="B";cb["sound_speed_m_s"]=1560.0
        pa=joinpath(dir,"A.json");pb=joinpath(dir,"B.json")
        C.write_json(pa,ca);C.write_json(pb,cb)
        L.set_device__ophk_acoustic_D3_clicker_A(calibration_path=pa)
        @test L.c==ca["sound_speed_m_s"]
        @test L.default_c==L.c
        @test L.imsize==(1520,2704)
        @test L.ref_channel==3
        geometry=C.clicker_geometry();fs=192000.0
        angles=[0.15 -0.1; -0.18 0.2; 0.0 0.0]
        directions=reduce(vcat,[permutedims([sin(a)*cos(e),sin(e),cos(a)*cos(e)]) for (a,e) in eachrow(angles)])
        tdoas=reduce(vcat,[permutedims((-(geometry.-permutedims(geometry[3,:]))*u/ca["sound_speed_m_s"]+ca["relative_delay_offsets_s"]).*fs) for u in eachrow(directions)])
        recovered,residual=L.tdoa2dir(tdoas,permutedims(geometry),fs;return_residual=true)
        @test recovered≈angles atol=1e-12
        @test maximum(abs.(residual))<1e-20
        direct=C.project_tdoas(ca,tdoas;fs=fs)
        expected=hcat(getproperty.(direct,:x),getproperty.(direct,:y))
        @test L.angle2px(recovered)≈expected atol=1e-9
        @test L.angle2pxd(rad2deg.(recovered))≈expected atol=1e-9
        @test L.tdoa2dir_ls(tdoas,permutedims(geometry),fs)≈angles atol=1e-12
        @test L.angle2px(recovered,L.fov_angle,(760,1352))≈expected./2 atol=1e-9
        @test all(isnan,L.tdoa2dir([900.0 0 0],permutedims(geometry),fs))
        @test_throws ArgumentError L.tdoa2dir(tdoas,permutedims(geometry).*2,fs)
        # Failed activation must preserve the previous active preset.
        @test_throws Exception L.set_device__ophk_acoustic_D3_clicker_B(calibration_path=joinpath(dir,"missing.json"))
        @test L.active_clicker_calibration()["device"]=="A"
        L.set_device__ophk_acoustic_D3_clicker_B(calibration_path=pb)
        @test L.c==1560.0
        @test L.active_clicker_calibration()["device"]=="B"
        L.set_device__calf_hk()
        @test L.active_clicker_calibration()===nothing
        @test L.c==1540
        @test L.default_c==1540
        @test L.default_getTDOA_func===L.get_tdoa_minmax
        @test L.default_window==(-5000:10000)
        @test L.angle2px([0.0 0.0])==[1920.0 1080.0]
        @test L.angle2px([0.0 0.0];calibration=nothing)==[1920.0 1080.0]
        bad=deepcopy(ca);bad["rotation_array_to_camera"]=[[1.,0,0],[0,-1.,0],[0,0,1.]]
        C.write_json(pa,bad)
        @test_throws ErrorException C.load_calibration(pa;device="A")
        # The selected fitted estimator must be the one used by legacy callers.
        time=collect(0:1000).-300.0
        waveform=hcat([exp.(-((time.-delay)./5).^2).*cos.(2pi*0.16.*(time.-delay)) for delay in (4.0,-2.0,0.0)]...)
        for method in ("gcc","srp")
            matched=deepcopy(ca);matched["method"]=method
            matched["estimator_spec"]=C._mr_estimator_spec(method;fs=fs)
            C.write_json(pa,matched);L.set_device__ophk_acoustic_D3_clicker_A(calibration_path=pa)
            @test L.impulsive_band_pass==[1000,Inf]
            @test L.default_window==(-300:700)
            actual_angles,actual_delays=L.detection2angle(waveform,[301],permutedims(geometry);fs=fs)
            expected_delays=C.estimate_method_tdoas(matched,waveform,[301];fs=fs)
            @test actual_delays==expected_delays
            direct=C.project_tdoas(matched,expected_delays;fs=fs)
            expected=hcat(getproperty.(direct,:x),getproperty.(direct,:y))
            @test L.angle2px(actual_angles)≈expected atol=1e-9
            @test_throws ErrorException L.detection2angle(waveform,[301],permutedims(geometry);fs=96000)
        end
        L.set_device__calf_hk()
        @test L.active_clicker_calibration()===nothing
    end
end
