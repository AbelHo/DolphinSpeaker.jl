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

const RESTORED_GLOBALS=(:c,:default_c,:default_getTDOA_func,:default_window,:imsize,:fov_angle,:FILE_device_ID,:impulsive_band_pass)
snapshot(L)=Dict(n=>getglobal(L,n) for n in RESTORED_GLOBALS)
same(L,snap)=L.active_clicker_calibration()===nothing && all(getglobal(L,n)===v for (n,v) in snap)

@testset "packaged independent A/B method-matched presets" begin
    L=CalibratedLegacyTest; C=ClickerCalibration
    for (device,method) in (("A","gcc"),("B","srp"))
        before=snapshot(L)
        L.with_clicker_calibration(device) do
            calibration=L.active_clicker_calibration()
            @test calibration["device"]==device
            @test calibration["method"]==method
            @test calibration["fit_status"]=="fitted_weak_identifiability"
            @test L.impulsive_band_pass==[1000,Inf]
            @test L.imsize==(1520,2704)
            tdoas=[0.0 0.0 0.0]
            angles=L.tdoa2dir(tdoas,permutedims(C.clicker_geometry()),192000.0)
            expected=only(C.project_tdoas(calibration,tdoas;fs=192000.0))
            @test L.angle2px(angles)≈[expected.x expected.y] atol=1e-9
        end
        @test same(L,before)
    end
    @test L.with_clicker_calibration(()->42,"A")==42
    @test L.with_clicker_calibration(()->L.active_clicker_calibration()["device"],:B)=="B"
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
        geometry=C.clicker_geometry();fs=192000.0
        before=snapshot(L)
        L.with_clicker_calibration("A";calibration_path=pa) do
            @test L.c==1540 && L.default_c==1540   # uncalibrated globals are left alone
            @test L.imsize==(1520,2704)
            @test L.ref_channel==3
            angles=[0.15 -0.1; -0.18 0.2; 0.0 0.0]
            directions=reduce(vcat,[permutedims([sin(a)*cos(e),sin(e),cos(a)*cos(e)]) for (a,e) in eachrow(angles)])
            tdoas=reduce(vcat,[permutedims((-(geometry.-permutedims(geometry[3,:]))*u/ca["sound_speed_m_s"]+ca["relative_delay_offsets_s"]).*fs) for u in eachrow(directions)])
            recovered,residual=L.tdoa2dir(tdoas,permutedims(geometry),fs;return_residual=true)
            @test recovered≈angles atol=1e-12
            @test maximum(abs.(residual))<1e-20
            direct=C.project_tdoas(ca,tdoas;fs=fs)
            expected=hcat(getproperty.(direct,:x),getproperty.(direct,:y))
            @test L.angle2px(recovered)≈expected atol=1e-9
            # angle2pxd and tdoa2dir_ls are not on the pipeline path: they stay uncalibrated while a calibration is active.
            # (tdoa2dir_ls on a non-planar array falls back to tdoa2dir, which would hit the D3 geometry check.)
            plain=rad2deg.(recovered) ./ [L.fov_angle[1]/2 -L.fov_angle[2]/2] .* [L.imsize[2]/2 L.imsize[1]/2] .+ [L.imsize[2]/2 L.imsize[1]/2]
            @test L.angle2pxd(rad2deg.(recovered))==plain
            @test !isapprox(L.angle2pxd(rad2deg.(recovered)),expected;atol=1e-3)
            @test !isapprox(L.tdoa2dir_ls(tdoas,permutedims(geometry),fs),angles;atol=1e-6)
            @test L.angle2px(recovered,L.fov_angle,(760,1352))≈expected./2 atol=1e-9
            @test all(isnan,L.tdoa2dir([900.0 0 0],permutedims(geometry),fs))
            @test_throws ArgumentError L.tdoa2dir(tdoas,permutedims(geometry).*2,fs)
            # Nested block: B inside A, and A is active again afterwards.
            L.with_clicker_calibration("B";calibration_path=pb) do
                @test L.active_clicker_calibration()["sound_speed_m_s"]==1560.0
                @test L.active_clicker_calibration()["device"]=="B"
            end
            @test L.active_clicker_calibration()["device"]=="A"
        end
        @test same(L,before)
        @test L.c==1540
        @test L.default_getTDOA_func===L.get_tdoa_minmax
        @test L.default_window==(-5000:10000)
        @test L.angle2px([0.0 0.0])==[1920.0 1080.0]
        @test L.angle2px([0.0 0.0];calibration=nothing)==[1920.0 1080.0]

        # An exception inside the block still restores everything.
        @test_throws ErrorException L.with_clicker_calibration("A";calibration_path=pa) do
            @test L.active_clicker_calibration()["device"]=="A"
            error("boom")
        end
        @test same(L,before)

        # A bad path fails before any global changes.
        @test_throws Exception L.with_clicker_calibration(()->nothing,"A";calibration_path=joinpath(dir,"missing.json"))
        @test same(L,before)
        bad=deepcopy(ca);bad["rotation_array_to_camera"]=[[1.,0,0],[0,-1.,0],[0,0,1.]]
        C.write_json(pa,bad)
        @test_throws ErrorException L.with_clicker_calibration(()->nothing,"A";calibration_path=pa)
        @test same(L,before)

        # The selected fitted estimator must be the one used by legacy callers.
        time=collect(0:1000).-300.0
        waveform=hcat([exp.(-((time.-delay)./5).^2).*cos.(2pi*0.16.*(time.-delay)) for delay in (4.0,-2.0,0.0)]...)
        for method in ("gcc","srp")
            matched=deepcopy(ca);matched["method"]=method
            matched["estimator_spec"]=C._mr_estimator_spec(method;fs=fs)
            C.write_json(pa,matched)
            L.with_clicker_calibration("A";calibration_path=pa) do
                @test L.impulsive_band_pass==[1000,Inf]
                @test L.window_impulsive==(-300:700)
                # The pipeline passes window=window_impulsive (aspod.jl); the fitted estimator checks it.
                actual_angles,actual_delays=L.detection2angle(waveform,[301],permutedims(geometry);fs=fs,window=L.window_impulsive)
                expected_delays=C.estimate_method_tdoas(matched,waveform,[301];fs=fs)
                @test actual_delays==expected_delays
                direct=C.project_tdoas(matched,expected_delays;fs=fs)
                expected=hcat(getproperty.(direct,:x),getproperty.(direct,:y))
                @test L.angle2px(actual_angles)≈expected atol=1e-9
                @test_throws ErrorException L.detection2angle(waveform,[301],permutedims(geometry);fs=96000,window=L.window_impulsive)
            end
            @test same(L,before)
        end
    end
end

# Runs last: it leaves the test module on a plain D3 setter.
@testset "persistent A/B setters switch themselves off on any later setter" begin
    L=CalibratedLegacyTest; C=ClickerCalibration
    tdoas=[0.0 0.0 0.0]; receivers=permutedims(C.clicker_geometry()); fs=192000.0
    original_tdoa=L.default_getTDOA_func
    for (device,setter) in (("A",L.set_device__ophk_acoustic_D3_clicker_A),("B",L.set_device__ophk_acoustic_D3_clicker_B))
        setter()
        calibration=L.active_clicker_calibration()
        @test calibration["device"]==device
        @test L.FILE_device_ID=="ophk_acoustic_D3_clicker_$(device)"
        @test L.default_getTDOA_func===L._clicker_tdoas
        expected=only(C.project_tdoas(calibration,tdoas;fs=fs))
        @test L.angle2px(L.tdoa2dir(tdoas,receivers,fs))≈[expected.x expected.y] atol=1e-9
        L.set_device__ophk_acoustic_D3()
        @test L.active_clicker_calibration()===nothing
        @test L.impulsive_band_pass==[5000,Inf]
        @test L.angle2px([0.0 0.0])==[L.imsize[2]/2 L.imsize[1]/2]
    end
    # A switched-off wrapper calls the estimator that was in place before the first activation.
    @test L.CLICKER_UNCALIBRATED_TDOA===original_tdoa
    data=sin.(0.01 .* (1:1001)) .* [1.0 0.5 0.2]
    @test L.default_getTDOA_func(data,[400];window=-50:50,ref_channel=3)==original_tdoa(data,[400];window=-50:50,ref_channel=3)
    # A block inside a persistent setting restores it afterwards.
    L.set_device__ophk_acoustic_D3_clicker_A()
    L.with_clicker_calibration(:B) do
        @test L.active_clicker_calibration()["device"]=="B"
    end
    @test L.active_clicker_calibration()["device"]=="A"
    # A bad path changes nothing.
    L.set_device__ophk_acoustic_D3()
    @test_throws Exception L.set_device__ophk_acoustic_D3_clicker_A(calibration_path="missing.json")
    @test L.active_clicker_calibration()===nothing
    @test L.impulsive_band_pass==[5000,Inf]
end
