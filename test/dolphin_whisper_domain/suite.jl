# Source bodies are evaluated unchanged, only selected declarations avoid unsafe includes.
module DomainSource
using DSP, FFTW, Optim, SignalAnalysis, SignalAnalysis.Units
using WAV, FFMPEG, JSON, SHA, Dates, LinearAlgebra, Statistics
const ROOT=Main.DOMAIN_ROOT
const DECLARATIONS=Any[]
function declaration_name(original)
    ex=original
    ex isa Expr && ex.head==:macrocall && ex.args[end] isa Expr && (ex=ex.args[end])
    ex isa Expr && ex.head in (:function,:(=)) || return nothing
    signature=ex.args[1]
    while signature isa Expr && signature.head in (:where,:(::))
        signature=signature.args[1]
    end
    signature isa Expr && signature.head==:call ? signature.args[1] : nothing
end
function select_declarations(file,names)
    path=joinpath(ROOT,"src",file);loaded=Symbol[]
    for original in Meta.parseall(read(path,String);filename=path).args
        name=declaration_name(original)
        name in names || continue
        digest=bytes2hex(sha256(repr(original)))
        Core.eval(@__MODULE__,original)
        push!(loaded,name)
        push!(DECLARATIONS,Dict("file"=>file,"name"=>string(name),
            "file_sha256"=>bytes2hex(sha256(read(path))),"original_ast_sha256"=>digest,
            "evaluated_ast_sha256"=>digest,"adaptation"=>"none"))
    end
    Set(loaded)==Set(names) || error("Selected declaration missing in $file: $(setdiff(names,loaded))")
end
select_declarations("media_info.jl",[:get_ffmpeg_metadata,:get_fps,:get_framerate,:get_whatever,:get_videos_audiodata])
select_declarations("audio.jl",[:fname2dt_zoomf6,:readAudio,:extract_clips])
include(joinpath(ROOT,"src","localization.jl")) # audited full dsp.jl/config.jl only
set_device__ophk_acoustic_D3()
# Device filename parsing is disabled in the fixture module; generated media carry explicit INFO names.
DEFAULT_fname2timestamp_func=nothing
select_declarations("synchronization.jl",[:findTrigger,:find_vid_vs_audio_syncdiff_timesegment])
BLAS.set_num_threads(1);FFTW.set_num_threads(1)
end

using DSP, FFTW, LinearAlgebra, Statistics, Logging
BLAS.set_num_threads(1);FFTW.set_num_threads(1)
const DOMAIN_MEDIA_RESULTS=Any[]

function write_pcm24_fixture(path)
    counts=Int32[-8388608,-4194304,-1,0,1,4194304,8388607]
    bytes=3length(counts);padding=isodd(bytes) ? 1 : 0
    open(path,"w") do io
        write(io,"RIFF");write(io,htol(UInt32(36+bytes+padding)));write(io,"WAVEfmt ")
        write(io,htol(UInt32(16)));write(io,htol(UInt16(1)));write(io,htol(UInt16(1)))
        write(io,htol(UInt32(192000)));write(io,htol(UInt32(192000*3)))
        write(io,htol(UInt16(3)));write(io,htol(UInt16(24)));write(io,"data")
        write(io,htol(UInt32(bytes)))
        for count in counts
            unsigned=reinterpret(UInt32,count)
            for shift in (0,8,16);write(io,UInt8((unsigned>>shift)&0xff));end
        end
        padding==1 && write(io,UInt8(0))
    end
    counts
end

function encode_known_offset_fixture(directory,fs_video,shift)
    fs_audio=48000;duration=14;n=fs_audio*duration
    times=collect(0:n-1)./fs_audio
    x=zeros(n)
    # Distinct carrier frequencies prevent another burst from matching a template.
    for (center,frequency,amplitude) in ((1.45,1733.,.45),(2.35,2371.,.35),
            (4.6,3119.,.5),(5.4,4223.,.4),(10.5,5397.,1.))
        x .+= amplitude .* exp.(-((times.-center)./.012).^2) .* cos.(2pi*frequency.*(times.-center))
    end
    pcm=round.(Int16,20000 .* x)
    stereo=hcat(pcm,-pcm)
    audio=joinpath(directory,"reference.wav")
    DomainSource.wavwrite(stereo,audio;Fs=fs_audio,nbits=16,
        chunks=DomainSource.wav_info_write(Dict(:INAM=>"synthetic known-offset reference")))
    resampled=DSP.resample(Float64.(stereo),fs_video/fs_audio;dims=1)
    base=round.(Int16,resampled)
    video_pcm=shift>=0 ? vcat(zeros(Int16,shift,2),base[1:end-shift,:]) :
        vcat(base[1-shift:end,:],zeros(Int16,-shift,2))
    video_wav=joinpath(directory,"video_audio.wav")
    DomainSource.wavwrite(video_pcm,video_wav;Fs=fs_video,nbits=16)
    video=joinpath(directory,"encoded.mp4")
    command=`-hide_banner -loglevel error -nostdin -y -filter_threads 1 -filter_complex_threads 1 -f lavfi -i color=c=black:s=64x48:r=30000/1001:d=14 -i $video_wav -map 0:v:0 -map 1:a:0 -c:v libx264 -threads 1 -pix_fmt yuv420p -c:a alac -write_tmcd 0 $video`
    try
        DomainSource.FFMPEG.exe(command)
    catch error
        Main.domain_failure(error,catch_backtrace();operation="encode known-offset sync media: $(repr(command))")
        rethrow()
    end
    (;audio,video,pcm=video_pcm,fs_audio,fs_video,shift)
end

const DOMAIN_TESTSET=@testset "Original Julia DolphinSpeaker domain" begin
    @testset "Pinned runtime and current source" begin
        @test VERSION==v"1.12.6"
        @test Threads.nthreads()==1
        @test BLAS.get_num_threads()==1
        @test pkgversion(DSP)==v"0.7.10"
        @test pkgversion(DomainSource.Optim)==v"1.13.3"
        for (file,expected) in DOMAIN_SOURCE_HASHES
            @test bytes2hex(sha256(read(joinpath(DOMAIN_ROOT,file))))==expected
        end
        @test length(DomainSource.DECLARATIONS)==11
        @test all(x->x["original_ast_sha256"]==x["evaluated_ast_sha256"],DomainSource.DECLARATIONS)
    end
    @testset "Original existing TDOA suite" begin
        @eval module OriginalTDOA
            # Exact existing 24-check test file, independently scoped to avoid globals.
            include(joinpath(Main.DOMAIN_ROOT,"test","test_tdoa.jl"))
        end
    end
    @testset "Native WAV units and source edge classifications" begin
        file=joinpath(DOMAIN_OUTPUT,"pcm24_edges.wav")
        counts=write_pcm24_fixture(file)
        native,fs,bits,_,_=DomainSource.readAudio(file;fname2timestamp_func=nothing)
        @test fs==192000 && bits==24
        @test size(native)==(7,1)
        @test vec(native)==counts
        normalized,_,_,_=DomainSource.wavread(file)
        @test vec(normalized)==counts./(2.0^23-1)
        data=reshape(Float64.(1:900),300,3)
        clips=DomainSource.extract_clips(data,[1,151],100)
        @test size.(clips)==[(51,3),(100,3)]
        @test_throws BoundsError DomainSource.extract_clips(data,[300],100)
        @test_throws BoundsError DomainSource.extract_clips(data,Int[],100)
    end
    @testset "DSP signed correlation, ties, rate and precision controls" begin
        @test DSP.finddelay([0.,0.,1.,2.,3.],[1.,2.,3.])==2
        @test DSP.finddelay([1.,2.,3.],[0.,0.,1.,2.,3.])==-2
        @test DSP.finddelay([1.,0.,1.],[0.,1.,0.])==1
        @test DSP.finddelay(zeros(5),zeros(5))==0
        delay,peak,_=DomainSource.finddelay2([-1.,-2.,-3.],[1.,2.,3.];flag_norm_rms=true)
        @test delay==0 && peak≈-1
        _,peak,correlation=DomainSource.finddelay2(zeros(2),zeros(3);flag_norm_rms=true)
        @test isnan(peak) && all(isnan,correlation)
        n=4096;fs=192000
        x=sin.(2pi*2000 .* (0:n-1)./fs)+sin.(2pi*40000 .* (0:n-1)./fs)
        filtered=DomainSource.filter_simple(x,[0,10000];fs=fs)
        @test eltype(filtered)==Float64 && size(filtered)==size(x)
        @test filtered≈DomainSource.filter_simple(x,[0,10000];fs=fs) rtol=0 atol=0
        implicit=DomainSource.filter_simple(x,[0,10000];fs=500000)
        @test maximum(abs,filtered.-implicit)>0.1
        signal=sin.(2pi.*(0:1023)./32)
        @test maximum(abs,abs.(DSP.hilbert(signal)).-1)<1e-12
        @test length(DSP.resample(signal,.5))==512
        @test DSP.finddelay(-signal,signal)==0 # integer estimator uses absolute correlation
    end
    @testset "Frozen extrema vote and projection origins" begin
        data=zeros(1000,3);p=500;window=-50:50
        for (ch,offset) in enumerate((-4,0,3))
            data[p+offset,ch]=-10
            data[p+offset+10,ch]=1
        end
        @test DomainSource.get_tdoa_minmax(data,[p];window=window)==[57 61 64]
        @test DomainSource.get_tdoa_minmax2(data,[p];window=window)==[47 51 54]
        @test DomainSource.angle2pxd([0. 0.],[80.,60.],(1520,2740))==[1370. 760.]
        @test DomainSource.angle2pxd([-40. 30.;40. -30.],[80.,60.],(1520,2740))≈[0. 0.;2740. 1520.]
        @test round(Int,1601/192000*(60000/1001))+1==1
        @test round(Int,1602/192000*(60000/1001))+1==2
        @test length(1/192000:1/(60000/1001):48049/192000)==16
    end
    @testset "Real encoded media, known sync offsets and legacy controls" begin
        for (fs_video,shift) in ((48000,1776),(24000,888),(24000,-504))
            directory=joinpath(DOMAIN_OUTPUT,"sync_$(fs_video)_$(shift)");mkpath(directory)
            fixture=encode_known_offset_fixture(directory,fs_video,shift)
            metadata=DomainSource.get_ffmpeg_metadata(fixture.video)
            @test Set(x["codec_type"] for x in metadata["streams"])==Set(["audio","video"])
            @test length(metadata["streams"])==2 # -write_tmcd 0 excludes a data/timecode stream
            @test DomainSource.get_fps(fixture.video)==30000//1001
            decoded,fs_read=DomainSource.get_videos_audiodata(fixture.video)
            @test fs_read==fs_video
            @test size(decoded)==size(fixture.pcm)
            @test decoded==fixture.pcm # lossless audio survives actual MP4 encoding/decoding
            result,confidence=DomainSource.find_vid_vs_audio_syncdiff_timesegment(
                fixture.video,fixture.audio;segment_inS=(1,3),fs=500000,
                flag_return_conf=true,flag_norm_rms=true)
            expected=(shift-1)/fs_video
            @test result==expected # source subtracts template_start, preserving one-sample origin
            @test result+1/fs_video≈shift/fs_video atol=1e-15
            @test isfinite(confidence) && confidence>1 # source Int16 abs2 wraps before widening
            native,_,_,_,_=DomainSource.readAudio(fixture.audio;fname2timestamp_func=nothing)
            down=DSP.resample(@view(native[1:fixture.fs_audio*3,1]),fs_video/fixture.fs_audio)
            template=down[fs_video:3fs_video]
            lag=fs_video+shift-1
            known_dot=sum(Float64.(decoded[lag+1:lag+length(template),1]).*template)
            promoted_energy=sum(abs2,Float64.(decoded[:,1]))
            # Independent signed Int16 modular arithmetic, avoiding source's abs2 implementation.
            wrapped_energy=sum(Int64(mod(Int64(v)^2+32768,65536)-32768) for v in decoded[:,1])
            @test wrapped_energy==sum(abs2,decoded[:,1])
            expected_legacy=known_dot/sqrt(wrapped_energy*sum(abs2,template))
            expected_promoted=known_dot/sqrt(promoted_energy*sum(abs2,template))
            @test confidence≈expected_legacy atol=1e-10 rtol=1e-10
            promoted_result,promoted_confidence=DomainSource.find_vid_vs_audio_syncdiff_timesegment(
                fixture.video,fixture.audio;segment_inS=(1,3),flag_return_conf=true,
                flag_norm_rms=true,norm_func=x->Float64.(x))
            @test promoted_result==expected
            @test 0<promoted_confidence<=1
            @test promoted_confidence≈expected_promoted atol=1e-10 rtol=1e-10
            if fs_video==24000 && shift==888
                logger=Test.TestLogger(min_level=Logging.Info)
                automatic=with_logger(logger) do
                    DomainSource.find_vid_vs_audio_syncdiff_timesegment(fixture.video,fixture.audio;
                        auto_segment_len=3,flag_return_conf=true,flag_norm_rms=true)
                end
                messages=[string(x.message) for x in logger.logs if occursin("Auto segment_inS",string(x.message))]
                @test messages==["Auto segment_inS: (4, 7)"] # legacy uses audio rate for video peak
                physical_start=round(Int,argmax(abs.(decoded[:,1]))/fs_video-1)
                @test physical_start==10 && physical_start!=4
                @test automatic[1]==expected
            end
            push!(DOMAIN_MEDIA_RESULTS,Dict("video"=>fixture.video,"audio"=>fixture.audio,
                "video_sample_rate_hz"=>fs_video,"audio_sample_rate_hz"=>fixture.fs_audio,
                "injected_shift_samples"=>shift,"injected_offset_seconds"=>shift/fs_video,
                "legacy_offset_seconds"=>result,"legacy_confidence"=>confidence,
                "promoted_confidence"=>promoted_confidence,"wrapped_Int16_energy"=>wrapped_energy,
                "promoted_Float64_energy"=>promoted_energy,
                "video_sha256"=>bytes2hex(sha256(read(fixture.video))),
                "audio_sha256"=>bytes2hex(sha256(read(fixture.audio)))))
        end
    end
end
counts=Test.get_test_counts(DOMAIN_TESTSET)
report=Dict("julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
    "passed"=>counts.passes+counts.cumulative_passes,"failed"=>counts.fails+counts.cumulative_fails,
    "errors"=>counts.errors+counts.cumulative_errors,"broken"=>counts.broken+counts.cumulative_broken,
    "source_sha256"=>DOMAIN_SOURCE_HASHES,"selected_declarations"=>DomainSource.DECLARATIONS,
    "complete_includes"=>["src/localization.jl","src/dsp.jl","src/config.jl","test/test_tdoa.jl (isolated module)"],
    "excluded_includes"=>["audio.jl/PAM.jl/include chain","synchronization.jl/include chain","run scripts","pinger/calibration work","test/test_synchronization.jl invalid empty-MP4 fixture","test_multitdoa broad memory probes"],
    "operational_fixture_default"=>"Generated sync WAVs carry synthetic INFO names; DEFAULT_fname2timestamp_func=nothing in fixture module only; manual PCM24 fixture intentionally has no timestamp",
    "media"=>DOMAIN_MEDIA_RESULTS)
write(joinpath(DOMAIN_OUTPUT,"report.json"),JSON.json(report,2))
println("DOMAIN_REPORT=",joinpath(DOMAIN_OUTPUT,"report.json"))
