# Native preparation for ClickerCalibration. No import-time work.
# All samples, frame indices and times exposed in tables are zero-origin.
# Float64 PCM normalized by WAV.jl; delays in recorder samples, positive=later.
const PREPARE_SCHEMA = "clicker-prepare-v1"
const PREPARE_COLUMNS = [:event_id,:device,:audio_file,:sample_index,:audio_time_s,:recording_time_s,:fs,
    :video_file,:video_time_s,:video_global_time_s,:frame_index,:frame_pts_s,:frame_id,:image_path,:width,:height,
    :sync_error_s,:sync_status,:minmax2_d1,:minmax2_d2,:minmax2_d3,:minmax_d1,:minmax_d2,:minmax_d3,
    :gcc_d1,:gcc_d2,:gcc_d3,:gcc_closure_samples,:peak_envelope,:threshold,:window_status,:window_source_segments,:snippet_file,:snippet_byte_offset,
    :snippet_samples,:snippet_channels,:snippet_peak_index,:snippet_sha256]
const PREPARE_FRAME_COLUMNS = [:frame_id,:video_file,:frame_index,:frame_pts_s,:video_global_time_s,:image_path,:width,:height]

_prep_table(cols) = DataFrame([c => Any[] for c in cols])
function _prep_atomic_json(path, value)
    tmp = path * ".tmp-$(getpid())"
    write_json(tmp, value)
    mv(tmp, path; force=true)
end
function _prep_atomic_csv(path, value)
    tmp = path * ".tmp-$(getpid())"
    CSV.write(tmp, value)
    open(tmp,"r+") do io
        ccall(:fsync,Cint,(Cint,),Base.fd(io))
    end
    mv(tmp, path; force=true)
end
_prep_run(cmd; operation, capture=false) = run_logged(cmd; operation=operation, capture=capture)

"""Read RIFF chunks without loading PCM. BEXT time_reference is diagnostic only.
Data length, not BEXT increments, defines the sample-concatenated timeline.
"""
function prep_wav_info(path)
    open(path, "r") do io
        String(read(io,4)) == "RIFF" || error("Only RIFF WAV supported: $path")
        read(io, UInt32); String(read(io,4)) == "WAVE" || error("Not WAVE")
        d = Dict{String,Any}("path"=>abspath(path), "bext_time_reference"=>nothing)
        while !eof(io) && position(io)+8 <= filesize(path)
            id = String(read(io,4)); n = Int(ltoh(read(io,UInt32))); at = position(io)
            if id == "fmt "
                d["format"] = Int(ltoh(read(io,UInt16)))
                d["channels"] = Int(ltoh(read(io,UInt16)))
                d["fs"] = Int(ltoh(read(io,UInt32)))
                read(io,UInt32); d["block_align"] = Int(ltoh(read(io,UInt16)))
                d["bits"] = Int(ltoh(read(io,UInt16)))
            elseif id == "data"
                d["data_offset"] = at; d["data_bytes"] = min(n,filesize(path)-at)
            elseif id == "bext" && n >= 346
                seek(io,at+338)
                lo=UInt64(ltoh(read(io,UInt32))); hi=UInt64(ltoh(read(io,UInt32)))
                d["bext_time_reference"] = lo | (hi << 32)
            end
            seek(io,min(at+n+isodd(n),filesize(path)))
        end
        haskey(d,"data_bytes") && haskey(d,"block_align") || error("Missing WAV data/fmt: $path")
        d["data_bytes"] % d["block_align"] == 0 || error("Incomplete PCM frame: $path")
        d["samples"] = d["data_bytes"] ÷ d["block_align"]
        d["duration_s"] = d["samples"] / d["fs"]
        d
    end
end

"""Read only selected channels into samples×channels Float64, in bounded WAV chunks."""
function prep_read_selected(path, info; channels=[3,4,5], first_sample=0, samples=info["samples"]-first_sample, chunk_samples=1_920_000)
    first_sample >= 0 && first_sample+samples <= info["samples"] || error("WAV range outside PCM")
    out=zeros(Float64,samples,length(channels))
    for s in 0:chunk_samples:samples-1
        n=min(chunk_samples,samples-s)
        a,fs= WAV.wavread(path; subrange=(first_sample+s+1):(first_sample+s+n))
        fs == info["fs"] || error("WAV sample rate changed")
        out[s+1:s+n,:] .= @view a[:,channels]
    end
    out
end

function prep_boundaries(infos)
    result=Any[]; offset=0
    for i in eachindex(infos)
        infos[i]["start_sample"] = offset
        if i>1
            p=infos[i-1]; q=infos[i]
            a=prep_read_selected(p["path"],p; channels=1:p["channels"], first_sample=p["samples"]-1,samples=1)
            b=prep_read_selected(q["path"],q; channels=1:q["channels"],samples=1)
            delta = isnothing(p["bext_time_reference"]) || isnothing(q["bext_time_reference"]) ? nothing : Int128(q["bext_time_reference"])-Int128(p["bext_time_reference"])-p["samples"]
            push!(result,Dict("previous"=>p["path"],"next"=>q["path"],"pcm_start_sample"=>offset,
                "bext_minus_pcm_samples"=>delta,"boundary_pcm_equal"=>a==b,"samples_removed"=>0))
        end
        offset+=infos[i]["samples"]
    end
    result
end

function prep_minmax(snip; corrected=true, reference=3)
    usemax=count(j->maximum(@view(snip[:,j])) > (corrected ? -minimum(@view(snip[:,j])) : minimum(@view(snip[:,j]))),axes(snip,2)) > size(snip,2)/2
    t=[Float64(usemax ? argmax(@view(snip[:,j])) : argmin(@view(snip[:,j]))) for j in axes(snip,2)]
    t .- t[reference]
end

# FFT correlation with genuine fractional lag grid, strictly inside physical bounds.
function prep_gcc_delay(x,y; fs, maxlag, phat=true, band=(1000.0,0.98fs/2), guide=nothing, interp=8)
    n=nextpow(2,2length(x)); X=rfft(vcat(Float64.(x),zeros(n-length(x)))); Y=rfft(vcat(Float64.(y),zeros(n-length(y))))
    R=X .* conj.(Y)
    if phat
        R ./= abs.(R) .+ (1e-12maximum(abs,R)+floatmin(Float64))
        for i in eachindex(R)
            f=(i-1)*fs/n
            (f<band[1] || f>band[2]) && (R[i]=0)
        end
    end
    nu=n*interp; z=zeros(ComplexF64,nu÷2+1); z[1:length(R)]=R
    r=irfft(z,nu); m=min(floor(Int,maxlag*interp),nu÷2-1)
    lags=collect(-m:m)./interp; curve=r[mod.(-m:m,nu).+1]
    k=argmax(curve)
    if !isnothing(guide)
        choices=[i for i in 2:length(curve)-1 if curve[i]>=curve[i-1] && curve[i]>=curve[i+1] && curve[i]>=0.7curve[k]]
        isempty(choices) || (k=choices[argmin(abs.(lags[choices].-guide))])
    end
    frac=0.0
    if 1<k<length(curve)
        den=curve[k-1]-2curve[k]+curve[k+1]
        den==0 || (frac=clamp(0.5(curve[k-1]-curve[k+1])/den,-0.5,0.5))
    end
    clamp(lags[k]+frac/interp,-maxlag,maxlag),curve[k]
end
# Constrained all-pair least squares on the D3 physical-lag hexagon.
function prep_reconcile_pairs(tau,bound)
    A=[1.0 -1.0;1.0 0.0;0.0 1.0]
    z=A\tau
    maximum(abs,A*z)<=bound && return [z[1],z[2],0.0]
    vertices=[(bound,0.0),(bound,bound),(0.0,bound),(-bound,0.0),(-bound,-bound),(0.0,-bound)]
    best=zeros(2); score=Inf
    for i in 1:6
        p=collect(vertices[i]);q=collect(vertices[mod1(i+1,6)]);v=q-p
        av=A*v; u=clamp(dot(av,tau-A*p)/dot(av,av),0.0,1.0)
        candidate=p+u*v; cost=sum(abs2,A*candidate-tau)
        if cost<score
            best=candidate;score=cost
        end
    end
    [best[1],best[2],0.0]
end
function prep_gcc(snip; fs=192000, c_min=1400.0, reference=3)
    reference==3 || error("D3 calibration uses selected reference3, physical5")
    # D3 is equilateral 0.18m; no empirical margin is added to the physical bound.
    bound=0.18fs/c_min; pairs=((1,2),(1,3),(2,3)); tau=zeros(3); t=zeros(3)
    env=[abs.(hilbert(snip[:,j])) for j in 1:3]
    for (p,(j,k)) in enumerate(pairs)
        guide=prep_gcc_delay(env[j].-mean(env[j]),env[k].-mean(env[k]);fs=fs,maxlag=bound,phat=false)[1]
        tau[p]=prep_gcc_delay(snip[:,j],snip[:,k];fs=fs,maxlag=bound,guide=guide)[1]
        t[j]+=tau[p]/3; t[k]-=tau[p]/3
    end
    d=prep_reconcile_pairs(tau,bound)
    d,abs(tau[1]+tau[3]-tau[2])
end

"""D3 click candidates: fourth-order zero-phase HP1000, analytic amplitude Q75×12,
15000-sample spacing, no train-length rejection. Returns all threshold candidates,
including edge windows. Peaks.findmaxima supplies the original peak-spacing rule.
"""
function prep_detect(data; fs=192000, ratio=12.0, distance=15000, threshold=nothing)
    size(data,2)==3 || error("Expected physical channels 3–5")
    f=digitalfilter(Highpass(1000;fs=fs),Butterworth(4))
    filtered=similar(data)
    for j in 1:3
        filtered[:,j]=filtfilt(f,@view(data[:,j]))
    end
    envelope=abs.(hilbert(filtered[:,3]))
    limit=isnothing(threshold) ? quantile(envelope,0.75)*ratio : Float64(threshold)
    inds,vals=Peaks.findmaxima(envelope,distance)
    keep=findall(>(limit),vals)
    # Endpoint clicks remain candidates if an endpoint is itself above threshold
    # and dominates its physical spacing neighborhood.
    peaks=collect(Int,inds[keep]); heights=collect(Float64,vals[keep])
    for k in (1,length(envelope))
        if envelope[k]>limit && all(i->abs(i-k)>=distance,peaks) && envelope[k]>=maximum(@view envelope[max(1,k-distance+1):min(end,k+distance-1)])
            push!(peaks,k); push!(heights,envelope[k])
        end
    end
    order=sortperm(peaks)
    (;peaks=peaks[order],heights=heights[order],threshold=limit,filtered)
end

function _prep_window_sources(infos, global_sample, window)
    lo=global_sample+first(window); hi=global_sample+last(window)+1
    [Dict("audio_file"=>q["path"],"first_sample"=>max(lo,q["start_sample"])-q["start_sample"],
        "samples"=>min(hi,q["start_sample"]+q["samples"])-max(lo,q["start_sample"])) for q in infos
        if min(hi,q["start_sample"]+q["samples"])>max(lo,q["start_sample"])]
end

function _prep_cross_snip(filtered, peak, infos, file_index; window=-300:700)
    lo=peak+first(window); hi=peak+last(window); n=size(filtered,1)
    snip=filtered[max(1,lo):min(n,hi),:]
    status=(lo>=1 && hi<=n) ? "complete" : n<infos[file_index]["samples"] && hi>n ? "analysis_edge_clipped" : "recording_edge_clipped"
    f=digitalfilter(Highpass(1000;fs=infos[file_index]["fs"]),Butterworth(4))
    if lo<1 && file_index>1
        p=infos[file_index-1]; nread=min(p["samples"],4096+1-lo)
        a=prep_read_selected(p["path"],p;first_sample=p["samples"]-nread,samples=nread)
        a=hcat([filtfilt(f,a[:,j]) for j in 1:3]...)
        snip=vcat(a[end-(1-lo)+1:end,:],snip); status="cross_file_context"
    end
    if hi>n && n==infos[file_index]["samples"] && file_index<length(infos)
        p=infos[file_index+1]; nread=min(p["samples"],4096+hi-n)
        a=prep_read_selected(p["path"],p;samples=nread)
        a=hcat([filtfilt(f,a[:,j]) for j in 1:3]...)
        snip=vcat(snip,a[1:hi-n,:]); status="cross_file_context"
    end
    snip,status
end

function _prep_decode_audio(path,fs;ffmpeg="/usr/bin/ffmpeg",start=nothing,duration=nothing)
    args=String[ffmpeg,"-v","error","-nostdin"]
    isnothing(start) || append!(args,["-ss",string(max(0,start))])
    append!(args,["-i",path]); isnothing(duration) || append!(args,["-t",string(duration)])
    append!(args,["-map","0:a:0","-ac","1","-ar",string(fs),"-f","f64le","pipe:1"])
    bytes=_prep_run(Cmd(args);operation="decode camera audio $(basename(path))",capture=true)
    # Shared run_logged returns String or bytes. Preserve arbitrary binary bytes.
    reinterpret(Float64,bytes isa AbstractString ? collect(codeunits(bytes)) : bytes) |> collect
end
function _prep_probe(path; ffprobe="/usr/bin/ffprobe")
    raw=_prep_run(`$ffprobe -v error -show_streams -show_format -of json $path`;operation="probe $(basename(path))",capture=true)
    JSON.parse(String(raw))
end
function prep_video_info(path;ffprobe="/usr/bin/ffprobe")
    p=_prep_probe(path;ffprobe=ffprobe); v=only(filter(s->s["codec_type"]=="video",p["streams"]))
    v["width"]==2704 && v["height"]==1520 || error("Expected original 2704×1520 video: $path")
    a=only(filter(s->s["codec_type"]=="audio",p["streams"])); parse(Int,a["sample_rate"])==48000 || error("Camera must be 48kHz")
    raw=_prep_run(`$ffprobe -v error -select_streams v:0 -show_frames -show_entries frame=best_effort_timestamp_time -of csv=p=0 $path`;operation="frame PTS $(basename(path))",capture=true)
    pts=Float64[]
    for line in split(String(raw),'\n')
        x=tryparse(Float64,strip(first(split(line,','))))
        isnothing(x) || push!(pts,x)
    end
    isempty(pts) && error("No actual frame timestamps: $path")
    issorted(pts) || error("Nonmonotonic PTS")
    Dict{String,Any}("path"=>abspath(path),"width"=>v["width"],"height"=>v["height"],"pts"=>pts,
        "duration_s"=>parse(Float64,v["duration"]),"audio_duration_s"=>parse(Float64,a["duration"]),
        "audio_start_s"=>parse(Float64,get(a,"start_time","0")),"video_start_s"=>parse(Float64,get(v,"start_time","0")))
end

# Normalized sliding waveform correlation. Time returned is template's start in x.
function prep_wave_match(x,y; fs, min_separation_s=0.5)
    n=length(x); m=length(y); n>=m>=16 || error("Insufficient waveform search")
    yy=Float64.(y).-mean(y); ey=sum(abs2,yy)
    ey>eps(Float64) || return (time_s=NaN,score=0.0,ambiguity_ratio=0.0,alternate_time_s=NaN,alternate_score=0.0,polarity=0.0)
    N=nextpow(2,n+m-1)
    corr=irfft(rfft(vcat(x,zeros(N-n))).*rfft(vcat(reverse(yy),zeros(N-m))),N)[m:n]
    csum=vcat(0.0,cumsum(x)); c2=vcat(0.0,cumsum(abs2.(x)))
    for i in eachindex(corr)
        ex=max(c2[i+m]-c2[i]-(csum[i+m]-csum[i])^2/m,0.0)
        corr[i]=abs(corr[i])/sqrt(max(ex*ey,floatmin(Float64)))
    end
    k=argmax(corr); score=corr[k]; radius=round(Int,min_separation_s*fs)
    other=0.0;other_k=0
    for i in eachindex(corr)
        if abs(i-k)>radius && corr[i]>other
            other=corr[i];other_k=i
        end
    end
    (time_s=(k-1)/fs,score=score,ambiguity_ratio=score/max(other,eps(Float64)),
        alternate_time_s=other_k==0 ? NaN : (other_k-1)/fs,alternate_score=other,
        polarity=sign(dot(@view(x[k:k+m-1]),yy)))
end

"""Robust affine recorder→camera mapping from waveform anchors, with recorded
residuals and explicit failure. Clock metadata never determines synchronization.
"""
function prep_fit_sync(anchors; residual_limit_s=0.015,min_anchors=4,max_drift_ppm=2000.0)
    good=[a for a in anchors if get(a,"accepted",false)]
    length(good)>=min_anchors || return Dict{String,Any}("status"=>"insufficient_waveform_anchors","anchors"=>anchors)
    x=Float64[a["recording_time_s"] for a in good]; y=Float64[a["video_global_time_s"] for a in good]
    slopes=[(y[j]-y[i])/(x[j]-x[i]) for i in eachindex(x) for j in i+1:length(x) if abs(x[j]-x[i])>10]
    isempty(slopes) && return Dict{String,Any}("status"=>"insufficient_anchor_span","anchors"=>anchors)
    b=median(slopes); a=median(y.-b.*x); keep=abs.(y.-a.-b.*x).<=residual_limit_s
    count(keep)>=min_anchors || return Dict{String,Any}("status"=>"inconsistent_waveform_anchors","anchors"=>anchors)
    xx=x[keep]; yy=y[keep]; b=sum((xx.-mean(xx)).*(yy.-mean(yy)))/sum(abs2,xx.-mean(xx)); a=mean(yy)-b*mean(xx)
    residual=abs.(yy.-a.-b.*xx); valid=maximum(residual)<=residual_limit_s && abs(b-1)*1e6<=max_drift_ppm
    for (i,q) in enumerate(good)
        q["fit_residual_s"]=y[i]-a-b*x[i]; q["fit_inlier"]=keep[i]
    end
    Dict{String,Any}("status"=>valid ? "valid" : "drift_or_residual_failed","offset_s"=>a,"scale"=>b,
        "drift_ppm"=>(b-1)*1e6,"max_residual_s"=>maximum(residual),"rms_residual_s"=>sqrt(mean(abs2,residual)),
        "anchor_min_time_s"=>minimum(xx),"anchor_max_time_s"=>maximum(xx),"anchors"=>anchors)
end

# DSP 0.7 resample can return half an input sample fewer at high integer
# interpolation. Explicit zero extension and exact slicing preserve t=0 and the
# requested duration; only padding after the actual signal is discarded.
function prep_resample_exact(x,output_fs,input_fs)
    count=ceil(Int,length(x)*(output_fs//input_fs))
    y=resample(vcat(x,zeros(eltype(x),64)),output_fs//input_fs)
    length(y)>=count || error("Resampler returned insufficient output samples")
    y[1:count]
end

# Stretch recorder waveform onto the provisional camera clock. Start stays t=0;
# omitting this correction biases start anchors by half a template's drift.
function prep_warp_template(x,scale)
    isfinite(scale) && 0.99<=scale<=1.01 || error("Clock scale outside waveform-warp range")
    n=ceil(Int,length(x)*scale)
    y=resample(vcat(x,zeros(eltype(x),64)),Float64(scale))
    length(y)>=n || error("Clock-warp resampler returned insufficient output")
    y[1:n]
end

function _prep_recorder_segment(infos,start_s,duration_s; target_fs=48000,channel=1)
    fs=infos[1]["fs"]; lo=round(Int,start_s*fs); hi=lo+round(Int,duration_s*fs)
    pieces=Vector{Float64}[]
    for q in infos
        a=max(lo,q["start_sample"]); b=min(hi,q["start_sample"]+q["samples"])
        b>a || continue
        x=prep_read_selected(q["path"],q;channels=[channel],first_sample=a-q["start_sample"],samples=b-a)[:,1]
        push!(pieces,x)
    end
    isempty(pieces) && return Float64[]
    prep_resample_exact(reduce(vcat,pieces),target_fs,fs)
end

function prep_sync_diagnostics(anchors)
    valid=sort([q for q in anchors if get(q,"accepted",false) && haskey(q,"fit_residual_s")];by=q->q["recording_time_s"])
    groups=Any[]
    for field in ("audio_file","video_file")
        ids=unique([get(q,field,"") for q in valid])
        for id in ids
            group=[q for q in valid if get(q,field,"")==id]
            r=Float64[q["fit_residual_s"] for q in group]
            push!(groups,Dict("field"=>field,"source"=>id,"anchor_count"=>length(r),
                "median_residual_s"=>median(r),"max_abs_residual_s"=>maximum(abs,r),
                "inlier_count"=>count(q->get(q,"fit_inlier",false),group)))
        end
    end
    # Diagnose a coherent timing step; do not apply a piecewise correction.
    # At least four anchors on either side avoid a two-point affine overfit.
    best=nothing; mismatch=false
    if length(valid)>=8
        x=Float64[q["recording_time_s"] for q in valid]; y=Float64[q["video_global_time_s"] for q in valid]
        global_rss=sum(q->abs2(q["fit_residual_s"]),valid)
        linearfit(xx,yy)=begin
            b=sum((xx.-mean(xx)).*(yy.-mean(yy)))/sum(abs2,xx.-mean(xx))
            a=mean(yy)-b*mean(xx); (a,b,sum(abs2,yy.-a.-b.*xx))
        end
        for k in 4:length(valid)-4
            x[k]-x[1]>10 && x[end]-x[k+1]>10 || continue
            a1,b1,r1=linearfit(x[1:k],y[1:k]);a2,b2,r2=linearfit(x[k+1:end],y[k+1:end])
            midpoint=(x[k]+x[k+1])/2;gap=(a2+b2*midpoint)-(a1+b1*midpoint)
            ratio=(r1+r2)/max(global_rss,eps(Float64)); noise=sqrt((r1+r2)/length(valid))
            flag=ratio<0.2 && abs(gap)>max(0.005,6noise)
            if isnothing(best) || ratio<best["rss_ratio"]
                best=Dict("split_recording_time_s"=>midpoint,"offset_step_s"=>gap,"rss_ratio"=>ratio,
                    "piecewise_rms_s"=>noise,"coherent_mismatch"=>flag)
            end
            mismatch |= flag
        end
    end
    Dict("source_groups"=>groups,"best_chronological_step"=>best,"coherent_segment_mismatch"=>mismatch,
        "interpretation"=>"diagnostic only; no medium classification or piecewise clock correction")
end

function prep_synchronize(infos,videos; ffmpeg="/usr/bin/ffmpeg",sync_fs=2000,anchor_window_s=6.0,
    sync_min_score=0.08,sync_min_ambiguity=1.15,sync_residual_limit_s=0.015,kwargs...)
    # Coarse signals are decoded one source at a time. Camera PCM is never written
    # and video pictures are never concatenated. Recorder resampling stays native.
    camera=Float64[]
    for v in videos
        v["global_start_s"]=length(camera)/sync_fs
        y=_prep_decode_audio(v["path"],sync_fs;ffmpeg=ffmpeg)
        # AAC decoders may emit a final padded codec frame beyond the track
        # duration. Trim only declared camera track padding, never recorder PCM.
        expected=round(Int,get(v,"audio_duration_s",v["duration_s"])*sync_fs)
        length(y)>=expected || error("Camera decoder returned shorter than declared track duration")
        resize!(y,expected)
        append!(camera,y)
        v["decoded_audio_duration_s"]=length(y)/sync_fs
    end
    bp=digitalfilter(Bandpass(200,800;fs=sync_fs),Butterworth(3))
    camera=filtfilt(bp,camera)
    anchors=Dict{String,Any}[]
    # Each WAV gets independent early/mid/late waveform searches. More anchors
    # are then placed in every chapter using the provisional waveform mapping.
    starts=Float64[]
    for q in infos
        base=q["start_sample"]/q["fs"]; dur=q["duration_s"]
        append!(starts,base .+ [max(0.0,min(dur-anchor_window_s,f*dur)) for f in (0.1,0.35,0.65,0.9)])
    end
    coarse=Dict{String,Any}[]
    for t in starts
        x=_prep_recorder_segment(infos,t,anchor_window_s;target_fs=sync_fs)
        length(x)>=sync_fs || continue
        x=filtfilt(bp,x)
        hit=prep_wave_match(camera,x;fs=sync_fs)
        accepted=hit.score>=sync_min_score && hit.ambiguity_ratio>=sync_min_ambiguity
        push!(coarse,Dict("recording_time_s"=>t,"video_global_time_s"=>hit.time_s,
            "score"=>hit.score,"ambiguity_ratio"=>hit.ambiguity_ratio,
            "alternate_time_s"=>hit.alternate_time_s,"alternate_score"=>hit.alternate_score,"waveform_polarity"=>hit.polarity,"accepted"=>accepted))
    end
    provisional=prep_fit_sync(coarse;residual_limit_s=max(0.1,sync_residual_limit_s),min_anchors=3)
    if !haskey(provisional,"offset_s")
        return Dict{String,Any}("status"=>"coarse_sync_failed","coarse_anchors"=>coarse,"anchors"=>anchors)
    end
    a=provisional["offset_s"]; b=provisional["scale"]
    for v in videos
        dur=v["decoded_audio_duration_s"]
        for f in (0.15,0.5,0.85)
            t=(v["global_start_s"]+f*dur-a)/b
            0<=t<=sum(q["duration_s"] for q in infos)-anchor_window_s && push!(starts,t)
        end
    end
    # Refine twice: warp each recorder template onto the provisional camera
    # clock, refit affine timing, then re-warp/recheck once. Fixed cached camera
    # search windows preserve exact source timestamps across the two iterations.
    abs(b-1)<=0.01 || return Dict{String,Any}("status"=>"coarse_clock_scale_failed","coarse_anchors"=>coarse,"anchors"=>anchors)
    waveform_cache=Dict{Any,Any}();history=Any[];result=Dict{String,Any}();applied_scale=b
    finebp=digitalfilter(Bandpass(200,10000;fs=48000),Butterworth(3))
    for iteration in 1:2
        applied_scale=b;anchors=Dict{String,Any}[]
        for t in sort(unique(starts))
            expected=a+b*t
            vid=findlast(v->v["global_start_s"]<=expected,videos)
            isnothing(vid) && continue
            v=videos[vid]; localtime=expected-v["global_start_s"]
            margin=0.3; begin_s=localtime-margin
            begin_s>=0 && begin_s+anchor_window_s+2margin<=v["decoded_audio_duration_s"] || continue
            cached=get!(waveform_cache,(t,vid)) do
                cam=_prep_decode_audio(v["path"],48000;ffmpeg=ffmpeg,start=begin_s,duration=anchor_window_s+2margin)
                rec=_prep_recorder_segment(infos,t,anchor_window_s;target_fs=48000)
                (;cam=filtfilt(finebp,cam),rec=rec,begin_s=begin_s)
            end
            warped=prep_warp_template(cached.rec,applied_scale)
            length(cached.cam)>=length(warped) || continue
            hit=prep_wave_match(cached.cam,filtfilt(finebp,warped);fs=48000,min_separation_s=0.05)
            qidx=findlast(q->q["start_sample"]/q["fs"]<=t,infos)
            accepted=hit.score>=sync_min_score && hit.ambiguity_ratio>=sync_min_ambiguity && isfinite(hit.time_s)
            push!(anchors,Dict("recording_time_s"=>t,"video_global_time_s"=>v["global_start_s"]+cached.begin_s+hit.time_s,
                "audio_file"=>infos[qidx]["path"],"video_file"=>v["path"],"window_s"=>anchor_window_s,
                "score"=>hit.score,"ambiguity_ratio"=>hit.ambiguity_ratio,
                "alternate_time_s"=>hit.alternate_time_s,"alternate_score"=>hit.alternate_score,"waveform_polarity"=>hit.polarity,"accepted"=>accepted,
                "coarse_expected_video_s"=>expected,"sample_precision_s"=>1/48000,
                "template_clock_scale"=>applied_scale,"refinement_iteration"=>iteration))
        end
        result=prep_fit_sync(anchors;residual_limit_s=sync_residual_limit_s)
        push!(history,Dict("iteration"=>iteration,"template_clock_scale"=>applied_scale,
            "fit_status"=>result["status"],"offset_s"=>get(result,"offset_s",nothing),"scale"=>get(result,"scale",nothing)))
        haskey(result,"scale") || break
        a=result["offset_s"];b=result["scale"]
        abs(b-1)<=0.01 || break
    end
    result["refinement_history"]=history
    if haskey(result,"scale")
        span=max(result["anchor_max_time_s"]-result["anchor_min_time_s"],1.0)
        # Include unresolved clock-scale change and anchor quantization, because
        # a shared template-stretch bias can be invisible in fit residuals.
        scale_error=max(abs(result["scale"]-applied_scale),2*(result["max_residual_s"]+1/48000)/span)
        result["clock_scale_uncertainty"]=scale_error
        result["template_scale_bias_bound_s"]=anchor_window_s*scale_error/2
    end
    result["segment_diagnostics"]=prep_sync_diagnostics(anchors)
    if result["status"]=="valid" && result["segment_diagnostics"]["coherent_segment_mismatch"]
        result["status"]="coherent_segment_mismatch"
    end
    result["coarse_anchors"]=coarse; result["method"]="multi-window normalized waveform correlation, robust affine recorder-to-camera fit"
    result["global_video_timeline"]= "sum of declared camera audio-track durations; codec padding excluded; each source video retains its original local PTS"
    result["camera_rate_hz"]=48000; result["recorder_sync_physical_channel"]=1
    good=[q for q in anchors if get(q,"fit_inlier",false)]
    covered_audio=unique([q["audio_file"] for q in good]); covered_video=unique([q["video_file"] for q in good])
    result["covered_audio_files"]=covered_audio;result["covered_video_files"]=covered_video
    result["minimum_inlier_anchors_per_source"]=2
    coverage_ok=all(q->count(a->a["audio_file"]==q["path"],good)>=2,infos) && all(v->count(a->a["video_file"]==v["path"],good)>=2,videos)
    if result["status"]=="valid" && !coverage_ok
        result["status"]="incomplete_source_coverage"
    end
    result
end

function prep_match_frame(videos,t)
    i=findlast(v->v["global_start_s"]<=t,videos)
    isnothing(i) && return nothing
    v=videos[i]; localtime=t-v["global_start_s"]+get(v,"audio_start_s",0.0)
    0<=localtime<v["duration_s"] || return nothing
    pts=v["pts"]; k=searchsortedfirst(pts,localtime)
    k=k>length(pts) ? length(pts) : k==1 ? 1 : abs(pts[k-1]-localtime)<=abs(pts[k]-localtime) ? k-1 : k
    (;video=v,index=k-1,pts=pts[k],global_pts=v["global_start_s"]+pts[k]-get(v,"audio_start_s",0.0),localtime)
end

function prep_sync_uncertainty(sync,recording_time_s)
    residual=sync["max_residual_s"]; span=max(sync["anchor_max_time_s"]-sync["anchor_min_time_s"],1.0)
    outside=max(sync["anchor_min_time_s"]-recording_time_s,recording_time_s-sync["anchor_max_time_s"],0.0)
    # Empirical bound, not a calibrated statistical confidence interval.
    # Quantization and extrapolated slope uncertainty are always included.
    residual+1/48000+get(sync,"template_scale_bias_bound_s",0.0)+outside*max(get(sync,"clock_scale_uncertainty",0.0),2residual/span,1/(48000span))
end
function prep_frame_ambiguous(videos,t,error_s)
    center=prep_match_frame(videos,t); isnothing(center) && return true
    key(hit)=isnothing(hit) ? nothing : (hit.video["path"],hit.index)
    key(prep_match_frame(videos,t-error_s))!=key(center) || key(prep_match_frame(videos,t+error_s))!=key(center)
end

"""Extract exact zero-origin source frame indices in batches. One seek and bounded
sequential decode per batch, original 2704×1520 PNG, no concatenated video.
"""
function prep_png_dimensions(path)
    open(path,"r") do io
        read(io,8)==UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a] || error("Invalid PNG: $path")
        seek(io,16)
        (Int(ntoh(read(io,UInt32))),Int(ntoh(read(io,UInt32))))
    end
end
function prep_extract_frames(frames;ffmpeg="/usr/bin/ffmpeg",batch_size=24,max_batch_span_s=20.0)
    for video in unique(frames.video_file)
        rows=sort(findall(==(video),frames.video_file);by=i->frames.frame_index[i]); pos=1
        while pos<=length(rows)
            lastpos=pos
            while lastpos<length(rows) && lastpos-pos+1<batch_size && frames.frame_pts_s[rows[lastpos+1]]-frames.frame_pts_s[rows[pos]]<=max_batch_span_s
                lastpos+=1
            end
            chunk=rows[pos:lastpos]; firstpts=frames.frame_pts_s[first(chunk)]
            # Accurate input seek preserves original PTS via copyts. Filter uses
            # source timestamps, not a new local n after the seek.
            expr=join(["eq(pts\\,$(round(Int,frames.frame_pts_s[i]*1e6)))" for i in chunk],"+")
            # settb=1/1e6 makes the exact probed decimal timestamps comparable;
            # GoPro PTS are integral microseconds after rounding. Check counts.
            destdir=dirname(frames.image_path[first(chunk)]);mkpath(destdir)
            prefix=joinpath(destdir,".batch-$(getpid())-$(frames.frame_id[first(chunk)])")
            pattern=prefix*"-%04d.png"
            vf="settb=1/1000000,select='$expr'"
            cmd=`$ffmpeg -v error -nostdin -y -ss $(string(max(0,firstpts-1))) -copyts -i $video -vf $vf -fps_mode passthrough -frames:v $(length(chunk)) -start_number 0 $pattern`
            _prep_run(cmd;operation="extract source frames $video indices $(frames.frame_index[chunk])")
            for (j,i) in enumerate(chunk)
                src=prefix*"-"*lpad(j-1,4,'0')*".png"
                isfile(src) || error("Frame extraction missing expected source frame $(frames.frame_index[i])")
                prep_png_dimensions(src)==(frames.width[i],frames.height[i]) || error("Decoded source frame dimensions differ")
                mv(src,frames.image_path[i];force=true)
            end
            pos=lastpos+1
        end
    end
    nothing
end

# Exact filtered-snippet cache: little-endian Float64, column-major, no header.
_prep_snippet_bytes(snip)=collect(reinterpret(UInt8,htol.(reinterpret(UInt64,vec(Matrix{Float64}(snip))))))
const _PREP_VERIFIED_FILES=Dict{Any,Bool}()
const _PREP_VERIFIED_EVENT_TABLES=Dict{String,Any}()
function _prep_verify_file(path,expected_sha256)
    isfile(path) || error("Missing authenticated file: $path")
    before=stat(path); key=(abspath(path),expected_sha256,before.size,before.mtime,before.ctime,before.inode)
    if !haskey(_PREP_VERIFIED_FILES,key)
        sha256_file(path)==expected_sha256 || error("SHA256 mismatch: $path")
        after=stat(path)
        (before.size,before.mtime,before.ctime,before.inode)==(after.size,after.mtime,after.ctime,after.inode) || error("Authenticated file changed while hashing")
        _PREP_VERIFIED_FILES[key]=true
    end
    key
end
_prep_event_value(row,key::Symbol)=row isa AbstractDict ? (haskey(row,key) ? row[key] : row[string(key)]) : row[key]

"""
    read_event_snippet(output_dir,event_row;verify=true) -> Matrix{Float64}

Bounded read of the exact prepared filtered samples×3 snippet. Storage is
Float64 little-endian column-major with no header. Zero-origin byte offset,
shape, zero-origin peak position and payload SHA256 are in the event row.
verify=true checks authenticated whole-cache/events.csv digests, canonical row
identity/sample/fs/window/cache metadata, and each payload digest. Whole-file
checks are cached only while inode/size/mtime/ctime and expected hash agree.
No source PCM is re-filtered; cross-file/clipped context is exactly as prepared.
"""
function read_event_snippet(output_dir,event_row;verify=true)
    filename=String(_prep_event_value(event_row,:snippet_file))
    filename=="event_snippets.f64" || error("Unexpected snippet cache filename")
    path=joinpath(abspath(output_dir),filename)
    offset=_prep_event_value(event_row,:snippet_byte_offset)
    ns=_prep_event_value(event_row,:snippet_samples);nc=_prep_event_value(event_row,:snippet_channels)
    peak=_prep_event_value(event_row,:snippet_peak_index)
    offset isa Integer && offset>=0 && offset%8==0 || error("Malformed snippet byte offset")
    ns isa Integer && 1<=ns<=1001 && nc isa Integer && nc==3 || error("Malformed snippet shape")
    peak isa Integer && 0<=peak<ns || error("Malformed snippet peak index")
    count=ns*nc*8
    isfile(path) && offset<=filesize(path)-count || error("Snippet range outside cache")
    if verify
        manifest=read_json(joinpath(output_dir,"prepare_manifest.json"))
        cache=manifest["snippet_cache"]
        cache["dtype"]=="Float64" && cache["endianness"]=="little" && cache["layout"]=="column-major" || error("Unsupported snippet cache encoding")
        cache["bytes"]==filesize(path) || error("Snippet cache size mismatch")
        _prep_verify_file(path,cache["sha256"])
        output_entry=only(q for q in manifest["outputs"] if q["path"]==path)
        output_entry["sha256"]==cache["sha256"] || error("Inconsistent cache manifest identity")
        events_path=joinpath(abspath(output_dir),"events.csv")
        table_entry=only(q for q in manifest["outputs"] if q["path"]==events_path)
        _prep_verify_file(events_path,table_entry["sha256"])
        canonical=get!(_PREP_VERIFIED_EVENT_TABLES,table_entry["sha256"]) do
            table=CSV.read(events_path,DataFrame)
            length(unique(table.event_id))==nrow(table) || error("Duplicate event IDs")
            Dict(String(r.event_id)=>NamedTuple(r) for r in eachrow(table))
        end
        id=String(_prep_event_value(event_row,:event_id))
        haskey(canonical,id) || error("Unknown event identity")
        expected=canonical[id]
        for field in (:event_id,:device,:audio_file,:sample_index,:audio_time_s,:recording_time_s,:fs,
            :window_status,:window_source_segments,:snippet_file,:snippet_byte_offset,:snippet_samples,:snippet_channels,:snippet_peak_index,:snippet_sha256)
            isequal(_prep_event_value(event_row,field),expected[field]) || error("Event identity/cache metadata mismatch: $field")
        end
    end
    before=stat(path)
    bytes=open(path,"r") do io
        seek(io,offset);read(io,count)
    end
    length(bytes)==count || error("Truncated snippet payload")
    if verify
        after=stat(path)
        (before.size,before.mtime,before.ctime,before.inode)==(after.size,after.mtime,after.ctime,after.inode) || error("Snippet cache changed during read")
        bytes2hex(SHA.sha256(bytes))==_prep_event_value(event_row,:snippet_sha256) || error("Snippet payload SHA256 mismatch")
    end
    words=ltoh.(reinterpret(UInt64,bytes))
    Matrix(reshape(reinterpret(Float64,words),ns,nc))
end

function _prep_inventory(paths)
    result=Any[]
    for p in paths
        before=stat(p); digest=sha256_file(p); after=stat(p)
        (before.size,before.mtime,before.ctime,before.inode)==(after.size,after.mtime,after.ctime,after.inode) || error("File changed while hashing: $p")
        push!(result,Dict("path"=>abspath(p),"bytes"=>after.size,"sha256"=>digest,
            "mtime"=>after.mtime,"ctime"=>after.ctime,"inode"=>after.inode))
    end
    result
end
function _prep_resume_valid(path,request)
    isfile(path) || return nothing
    m=read_json(path)
    get(m,"status","")=="complete" && get(m,"request",nothing)==request || return nothing
    for q in get(m,"outputs",Any[])
        p=q["path"]
        isfile(p) && filesize(p)==q["bytes"] && sha256_file(p)==q["sha256"] || return nothing
    end
    m
end

"""
    prepare_device(input_dir,output_dir;device,resume=true,kwargs...) -> Dict

Native Julia D3 preparation. Physical channels 3–5/ref5; 192kHz five-channel
PCM. Each WAV is loaded sequentially, no click-train rejection. Main event table
keeps unmapped candidates. Camera synchronization uses waveform evidence across
all source WAVs/chapters; camera metadata alone is never accepted. Input SHA256,
actual source PTS and exact frame IDs are persisted. PNG extraction is bounded.

Useful bounded execution: `audio_limit_s=10.0` analyzes first ten seconds per WAV;
`extract_frames=false` writes mapping without pictures and cannot be complete.
`detection_only=true` skips camera decoding/probing for detector feasibility.
These are explicit request parameters and cannot reuse a full preparation receipt.
Resume authenticates source/content/config and every output; any mismatch reruns.
"""
function prepare_device(input_dir,output_dir;device,resume=true,ffmpeg="/usr/bin/ffmpeg",ffprobe="/usr/bin/ffprobe",
    audio_limit_s=nothing,extract_frames=true,detection_only=false,ratio=12.0,distance=15000,
    sync_min_score=0.08,sync_min_ambiguity=1.15,sync_residual_limit_s=0.015,kwargs...)
    try
        _prepare_device(input_dir,output_dir;device=device,resume=resume,ffmpeg=ffmpeg,ffprobe=ffprobe,
            audio_limit_s=audio_limit_s,extract_frames=extract_frames,detection_only=detection_only,
            ratio=ratio,distance=distance,sync_min_score=sync_min_score,sync_min_ambiguity=sync_min_ambiguity,
            sync_residual_limit_s=sync_residual_limit_s,kwargs...)
    catch e
        log_failure("prepare_device device=$device input=$input_dir output=$output_dir",e,catch_backtrace())
        rethrow()
    end
end
function _prepare_device(input_dir,output_dir;device,resume,ffmpeg,ffprobe,audio_limit_s,extract_frames,detection_only,
    ratio,distance,sync_min_score,sync_min_ambiguity,sync_residual_limit_s,kwargs...)
    isempty(kwargs) || error("Unknown preparation keywords: $(keys(kwargs))")
    input_dir=abspath(input_dir); output_dir=abspath(output_dir);mkpath(output_dir)
    islink(output_dir) && error("Output directory must not be symlink")
    wavs=sort(filter(p->lowercase(splitext(p)[2])==".wav",readdir(input_dir;join=true)))
    vids=sort(filter(p->lowercase(splitext(p)[2])==".mp4",readdir(input_dir;join=true)))
    isempty(wavs) && error("No recorder WAVs")
    infos=prep_wav_info.(wavs)
    all(q->q["channels"]==5 && q["fs"]==192000,infos) || error("D3 requires five channels at 192kHz")
    boundaries=prep_boundaries(infos)
    inputs=_prep_inventory(detection_only ? wavs : vcat(wavs,vids))
    config=Dict("device"=>string(device),"ratio"=>ratio,"distance"=>distance,"highpass_hz"=>1000,
        "window"=>[-300,700],"physical_channels"=>[3,4,5],"reference_physical_channel"=>5,
        "c_min_m_s"=>1400,"audio_limit_s"=>audio_limit_s,"extract_frames"=>extract_frames,
        "detection_only"=>detection_only,"sync_min_score"=>sync_min_score,"sync_min_ambiguity"=>sync_min_ambiguity,
        "sync_residual_limit_s"=>sync_residual_limit_s)
    environment=Dict("julia_version"=>string(VERSION),"DSP_version"=>string(Base.pkgversion(DSP)),
        "FFTW_version"=>string(Base.pkgversion(FFTW)),"WAV_version"=>string(Base.pkgversion(WAV)))
    request=Dict("environment"=>environment,"schema"=>PREPARE_SCHEMA,"inputs"=>inputs,"config"=>config,"source_sha256"=>sha256_file(@__FILE__))
    manifest=joinpath(output_dir,"prepare_manifest.json")
    if resume
        prior=_prep_resume_valid(manifest,request)
        isnothing(prior) || return prior
    end
    _prep_atomic_json(manifest,Dict("status"=>"incomplete","request"=>request,"started_utc"=>string(now(UTC))))
    events=_prep_table(PREPARE_COLUMNS); frames=_prep_table(PREPARE_FRAME_COLUMNS)
    thresholds=Any[]
    cachepath=joinpath(output_dir,"event_snippets.f64");cachetemp=cachepath*".tmp-$(getpid())"
    open(cachetemp,"w") do snippet_io
    for (i,q) in enumerate(infos)
        count=isnothing(audio_limit_s) ? q["samples"] : min(q["samples"],round(Int,audio_limit_s*q["fs"]))
        data=prep_read_selected(q["path"],q;samples=count)
        d=prep_detect(data;fs=q["fs"],ratio=ratio,distance=distance)
        data=nothing
        push!(thresholds,Dict("audio_file"=>q["path"],"threshold"=>d.threshold,"analyzed_samples"=>count,"candidate_count"=>length(d.peaks)))
        for (j,p) in enumerate(d.peaks)
            snip,edge=_prep_cross_snip(d.filtered,p,infos,i)
            t2=prep_minmax(snip); t1=prep_minmax(snip;corrected=false); tg,closure=prep_gcc(snip;fs=q["fs"])
            payload=_prep_snippet_bytes(snip);byte_offset=position(snippet_io);write(snippet_io,payload)
            sample=p-1; id="$(device)-$(splitext(basename(q["path"]))[1])-s$sample"
            row=Dict{Symbol,Any}(c=>missing for c in PREPARE_COLUMNS)
            merge!(row,Dict(:event_id=>id,:device=>string(device),:audio_file=>q["path"],:sample_index=>sample,
                :audio_time_s=>sample/q["fs"],:recording_time_s=>(q["start_sample"]+sample)/q["fs"],:fs=>q["fs"],
                :sync_status=>"sync_failed",:peak_envelope=>d.heights[j],:threshold=>d.threshold,:window_status=>edge,
                :window_source_segments=>JSON.json(_prep_window_sources(infos,q["start_sample"]+sample,isnothing(audio_limit_s) ? (-300:700) : (-300:min(700,count-1-sample)))),
                :snippet_file=>"event_snippets.f64",:snippet_byte_offset=>byte_offset,:snippet_samples=>size(snip,1),:snippet_channels=>3,
                :snippet_peak_index=>min(300,q["start_sample"]+sample),:snippet_sha256=>bytes2hex(SHA.sha256(payload)),
                :minmax2_d1=>t2[1],:minmax2_d2=>t2[2],:minmax2_d3=>t2[3],
                :minmax_d1=>t1[1],:minmax_d2=>t1[2],:minmax_d3=>t1[3],
                :gcc_d1=>tg[1],:gcc_d2=>tg[2],:gcc_d3=>tg[3],:gcc_closure_samples=>closure))
            push!(events,row)
        end
        d=nothing;GC.gc()
    end
    flush(snippet_io);ccall(:fsync,Cint,(Cint,),Base.fd(snippet_io))
    end
    mv(cachetemp,cachepath;force=true)
    # Cross-boundary peak spacing: retain the stronger candidate in each conflict;
    # preserve the original per-file candidate in rejected_candidates.csv.
    rejected=_prep_table(PREPARE_COLUMNS); keep=trues(nrow(events))
    order=sortperm(events.peak_envelope;rev=true)
    accepted=Int[]
    for i in order
        if any(j->abs(events.recording_time_s[i]-events.recording_time_s[j])*events.fs[i]<distance,accepted)
            keep[i]=false;push!(rejected,events[i,:])
        else
            push!(accepted,i)
        end
    end
    # Every audio candidate remains in events.csv; spacing conflicts are explicit
    # labels and remain available to reviewers, with no duplicate silent deletion.
    for i in findall(.!keep)
        events.window_status[i]=string(events.window_status[i],";cross_boundary_spacing_rejected")
    end
    videos=Dict{String,Any}[]; sync=Dict{String,Any}("status"=>"detection_only")
    if !detection_only
        isempty(vids) && error("No camera chapters")
        videos=prep_video_info.(vids;ffprobe=ffprobe)
        sync=prep_synchronize(infos,videos;ffmpeg=ffmpeg,sync_min_score=sync_min_score,
            sync_min_ambiguity=sync_min_ambiguity,sync_residual_limit_s=sync_residual_limit_s)
        framekeys=Dict{Tuple{String,Int},String}()
        for i in 1:nrow(events)
            if sync["status"]!="valid"
                events.sync_status[i]=(occursin("inconsistent",sync["status"]) || occursin("coherent",sync["status"])) ? "ambiguous" : "sync_failed"
                continue
            end
            t=sync["offset_s"]+sync["scale"]*events.recording_time_s[i]
            hit=prep_match_frame(videos,t)
            if isnothing(hit)
                events.sync_status[i]="outside_video";continue
            end
            v=hit.video;key=(v["path"],hit.index)
            id=get!(framekeys,key) do
                vidsha=only(q["sha256"] for q in inputs if q["path"]==v["path"])
                fid="$(device)-$(vidsha[1:12])-f$(hit.index)"
                push!(frames,(fid,v["path"],hit.index,hit.pts,hit.global_pts,joinpath(output_dir,"frames",fid*".png"),v["width"],v["height"]))
                fid
            end
            events.video_file[i]=v["path"];events.video_time_s[i]=hit.localtime;events.video_global_time_s[i]=t
            events.frame_index[i]=hit.index;events.frame_pts_s[i]=hit.pts;events.frame_id[i]=id
            events.image_path[i]=joinpath(output_dir,"frames",id*".png");events.width[i]=v["width"];events.height[i]=v["height"]
            err=prep_sync_uncertainty(sync,events.recording_time_s[i])
            events.sync_error_s[i]=err;events.sync_status[i]=prep_frame_ambiguous(videos,t,err) ? "ambiguous" : "mapped"
        end
        extract_frames && prep_extract_frames(frames;ffmpeg=ffmpeg)
    end
    # Store exact source identities for windows straddling split boundaries.
    qpath=joinpath(output_dir,"events.csv");fpath=joinpath(output_dir,"frames.csv")
    _prep_atomic_csv(qpath,events);_prep_atomic_csv(fpath,frames)
    _prep_atomic_csv(joinpath(output_dir,"rejected_candidates.csv"),rejected)
    _prep_atomic_json(joinpath(output_dir,"synchronization.json"),sync)
    meta=Dict("wav_files"=>infos,"boundaries"=>boundaries,"detector"=>thresholds,
        "videos"=>videos,"time_origin"=>"zero; PCM concatenated by exact byte count",
        "tdoa_units"=>"recorder samples; positive=later; relative selected3 physical5",
        "snippet_cache_encoding"=>"Float64 little-endian, column-major, zero-origin byte offsets, exact existing filtered snippets; no header",
        "sync_error_interpretation"=>"empirical max anchor residual plus48kHz quantization, unresolved template clock-scale bias and extrapolated slope error; ambiguous if interval selects different exact source frames",
        "edge_context"=>"neighbor PCM filtered independently with 4096 extra samples; edge candidates retained",
        "image_coordinates"=>"original 2704x1520, zero-origin x right y down")
    _prep_atomic_json(joinpath(output_dir,"prepare_metadata.json"),meta)
    paths=[cachepath,qpath,fpath,joinpath(output_dir,"rejected_candidates.csv"),joinpath(output_dir,"synchronization.json"),joinpath(output_dir,"prepare_metadata.json")]
    extract_frames && append!(paths,String.(frames.image_path))
    # Full successful preparation requires authenticated source-wide waveform sync.
    complete=!detection_only && extract_frames && sync["status"]=="valid" && isnothing(audio_limit_s)
    for source in inputs
        current=stat(source["path"])
        (current.size,current.mtime,current.ctime,current.inode)==(source["bytes"],source["mtime"],source["ctime"],source["inode"]) || error("Input changed during preparation: $(source["path"])")
    end
    sha256_file(@__FILE__)==request["source_sha256"] || error("Preparation source changed during execution")
    m=Dict("status"=>complete ? "complete" : "bounded_or_unmapped","request"=>request,"completed_utc"=>string(now(UTC)),
        "events_csv"=>qpath,"frames_csv"=>fpath,"event_count"=>nrow(events),"frame_count"=>nrow(frames),"mapped_count"=>count(==("mapped"),events.sync_status),
        "sync_status"=>sync["status"],"snippet_cache"=>Dict("path"=>cachepath,"sha256"=>sha256_file(cachepath),"bytes"=>filesize(cachepath),
            "dtype"=>"Float64","endianness"=>"little","layout"=>"column-major"),"outputs"=>_prep_inventory(paths))
    _prep_atomic_json(manifest,m)
    m
end
