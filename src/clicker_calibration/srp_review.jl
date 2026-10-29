# Post-first-fit estimator review. Included inside ClickerCalibration after
# physical_fit.jl and prepare.jl. Definitions only: no import-time I/O/work.

"""Fixed original Julia SRP formulation, explicitly frozen rather than taken
from a potentially modified live parent source. No heldout tuning is allowed.
"""
function _srv_spec()
    return Dict{String,Any}("schema_version"=>"clicker-estimator-review/1.0",
        "original_source"=>"source_snapshot/src/localization.jl",
        "original_source_sha256"=>"c35143f112a31dafd2cfdcbefa14f291c87377a93077d61a5aa34baae886f663",
        "original_functions"=>["gcc_curve","curve_at","get_tdoa_srp"],
        "sound_speed_m_s"=>1540.0,"geometry"=>"D3 equilateral side0.18m; selected physical3,4,5; ref selected3/physical5",
        "weighting"=>"PHAT; epsilon=1e-12*maximum(abs,cross_spectrum)+floatmin(Float64)",
        "band_hz"=>[1000.0,"0.98*fs/2"],"interp"=>8,"taper"=>0.0,
        "grid_deg"=>1.0,"grid_equation"=>"u,v in -1:sind(1):1, u^2+v^2<=1; positive array z branch",
        "refine"=>10,"refinement"=>"21x21 offsets range(-sind(1),sind(1)); strict > updates around coarse winner",
        "fft"=>"DSP.nextfastfft(2L); R=X*conj(Y); high-frequency spectrum zero padding to N*8; original normalized irfft",
        "lag_limit"=>"sensor pair distance/c*fs + 2 samples; ceil lag*8, capped at Nup/2-1",
        "lag_interpolation"=>"original curve_at linear interpolation with endpoint range guard; grid steering always physically feasible",
        "score"=>"mean of three normalized pair curves at winning physical direction; diagnostic, not probability or sound-speed truth",
        "delay_sign"=>"positive means channel later than reference, dt_i=-(r_i-r_ref) dot u /1540 * fs",
        "snippets"=>"read_event_snippet(...;verify=true): exact preparation-filtered cache bytes; original HP1kHz fourth-order filtfilt; inclusive -300:700 window; exact cached clipping/cross-file context retained",
        "projection"=>"Every method uses frozen evaluation_calibration through project_tdoas once; its relative-delay correction and rotation are applied once. No method-specific fit or heldout tuning.",
        "departure_from_original"=>"Silent/nonfinite/missing-channel snippets are explicit invalid estimates rather than arbitrary first-grid directions. Valid nonzero finite snippets use original fixed formulation.")
end

function _srv_curve(x,y;fs,maxlag)
    L=length(x);N=DSP.nextfastfft(2L)
    X=rfft([Float64.(x);zeros(N-L)]);Y=rfft([Float64.(y);zeros(N-L)])
    R=X.*conj.(Y);epsilon=1e-12*maximum(abs,R)+floatmin(Float64)
    R .*= 1.0 ./ (abs.(R).+epsilon)
    for i in eachindex(R)
        frequency=(i-1)*fs/N
        (frequency < 1000.0 || frequency > 0.98fs/2) && (R[i]=0)
    end
    Nup=N*8;Rup=zeros(ComplexF64,Nup÷2+1);Rup[1:length(R)]=R
    curve=irfft(Rup,Nup)./max(2sum(abs,R)/Nup,floatmin(Float64))
    M=min(ceil(Int,maxlag*8),Nup÷2-1)
    return ((-M:M)./8,curve[mod.(-M:M,Nup).+1])
end

function _srv_curve_at(curve,tau)
    lags,r=curve
    x=clamp((tau-first(lags))/step(lags)+1,1,length(r))
    k=min(floor(Int,x),length(r)-1)
    return r[k]+(x-k)*(r[k+1]-r[k])
end

function _srv_plan(fs)
    _pf_finite(fs) && fs > 1000.0/0.49 || throw(ArgumentError("SRP needs finite sample rate with a nonempty 1kHz..0.98Nyquist band"))
    geometry=clicker_geometry();pairs=[(1,2),(1,3),(2,3)]
    maxlags=[norm(geometry[j,:]-geometry[k,:])/1540.0*fs+2 for (j,k) in pairs]
    h=sind(1.0);grid=[(u,v) for u in -1:h:1 for v in -1:h:1 if u^2+v^2<=1]
    delays=zeros(3,length(grid))
    for (i,(u,v)) in enumerate(grid)
        arrival=-(geometry*[u,v,sqrt(max(1-u^2-v^2,0.0))])./1540.0.*fs
        for (p,(j,k)) in enumerate(pairs);delays[p,i]=arrival[j]-arrival[k];end
    end
    return (fs=Float64(fs),geometry=geometry,pairs=pairs,maxlags=maxlags,h=h,
        grid=grid,delays=delays,fine=range(-h,h;length=21))
end

function _srv_invalid(status)
    return (tdoas_samples=fill(NaN,3),uv=fill(NaN,2),score=NaN,
        pair_scores=fill(NaN,3),valid=false,status=string(status),grid_boundary=false)
end

"""
    review_srp_phat(snippet; fs=192000, plan=nothing)

Native Float64 original-formulation SRP-PHAT on an exact prepared samples×3
snippet (selected physical3–5). Fixed c1540m/s, PHAT1kHz..0.98Nyquist, interp8,
sine-direction grid1degree/refine10, no taper. Positive lag means later arrival.
Returns reference-relative sample delays, planar uv, mean/pair scores and
explicit validity. Silence/nonfinite or a silent selected channel is invalid;
no annotation or fitted parameter enters estimation. The fitted calibration is
applied separately by project_tdoas, exactly once. Scores do not measure source
sound-speed truth, correctness, range, or probability.
"""
function review_srp_phat(snippet::AbstractMatrix;fs=192000,plan=nothing)
    size(snippet,2)==3 || throw(DimensionMismatch("SRP snippet must have three selected channels"))
    size(snippet,1)>=4 || return _srv_invalid("too_short")
    all(isfinite,snippet) || return _srv_invalid("nonfinite_snippet")
    any(j->maximum(abs,view(snippet,:,j))==0.0,1:3) && return _srv_invalid("silent_channel")
    search=isnothing(plan) ? _srv_plan(fs) : plan
    search.fs==Float64(fs) || throw(ArgumentError("SRP plan sample-rate mismatch"))
    curves=[_srv_curve(view(snippet,:,j),view(snippet,:,k);fs=fs,maxlag=search.maxlags[p]) for (p,(j,k)) in enumerate(search.pairs)]
    all(curve->all(isfinite,last(curve)),curves) || return _srv_invalid("nonfinite_correlation")
    any(curve->maximum(abs,last(curve))==0.0,curves) && return _srv_invalid("no_inband_pair_energy")
    score(d)=sum(p->_srv_curve_at(curves[p],d[p]),1:3)
    winner=argmax(i->score(view(search.delays,:,i)),eachindex(search.grid))
    best=score(view(search.delays,:,winner));bu,bv=search.grid[winner];u0,v0=bu,bv
    for du in search.fine,dv in search.fine
        u,v=u0+du,v0+dv
        u^2+v^2<=1 || continue
        arrival=-(search.geometry*[u,v,sqrt(max(1-u^2-v^2,0.0))])./1540.0.*fs
        d=[arrival[j]-arrival[k] for (j,k) in search.pairs]
        value=score(d)
        value>best && ((best,bu,bv)=(value,u,v))
    end
    arrival=-(search.geometry*[bu,bv,sqrt(max(1-bu^2-bv^2,0.0))])./1540.0.*fs
    delays=arrival.-arrival[3]
    pair_scores=[_srv_curve_at(curves[p],arrival[j]-arrival[k]) for (p,(j,k)) in enumerate(search.pairs)]
    return (tdoas_samples=delays,uv=[bu,bv],score=best/3,pair_scores=pair_scores,
        valid=isfinite(best),status="ok",grid_boundary=bu^2+bv^2>0.98^2)
end

function _srv_metric(rows)
    valid=[row for row in eachrow(rows) if row.projection_valid]
    errors=Float64[row.error_px for row in valid]
    scores=Float64[value for value in rows.srp_score if isfinite(value)]
    return merge(_pf_error_summary(errors),Dict("n_population"=>nrow(rows),
        "n_estimator_valid"=>count(identity,rows.estimator_valid),
        "coverage"=>nrow(rows)==0 ? nothing : length(valid)/nrow(rows),
        "estimator_coverage"=>nrow(rows)==0 ? nothing : count(identity,rows.estimator_valid)/nrow(rows),
        "in_image_coverage"=>nrow(rows)==0 ? nothing : count(identity,rows.in_frame)/nrow(rows),
        "srp_score_summary"=>isempty(scores) ? nothing : Dict("n"=>length(scores),"minimum"=>minimum(scores),"median"=>median(scores),"p90"=>quantile(scores,0.9),"maximum"=>maximum(scores)),
        "srp_grid_boundary_count"=>count(identity,rows.srp_grid_boundary),
        "projection_status_counts"=>Dict(status=>count(==(status),rows.projection_status) for status in unique(rows.projection_status)),
        "estimator_status_counts"=>Dict(status=>count(==(status),rows.estimator_status) for status in unique(rows.estimator_status))))
end

function _srv_summarize(table)
    summaries=Dict{String,Any}()
    for scope in ("train","heldout","all_descriptive")
        selected=scope=="all_descriptive" ? table : table[table.split.==scope,:]
        methods=Dict{String,Any}()
        for method in ("minmax2","minmax","gcc","srp")
            methods[method]=_srv_metric(selected[selected.method.==method,:])
        end
        paired=Dict{String,Any}()
        for left in ("minmax2","minmax","gcc","srp"),right in ("minmax2","minmax","gcc","srp")
            left>=right && continue
            lrows=selected[(selected.method.==left).&selected.projection_valid,:]
            rrows=selected[(selected.method.==right).&selected.projection_valid,:]
            rmap=Dict(_pf_text(r.event_id)=>Float64(r.error_px) for r in eachrow(rrows))
            le=Float64[];re=Float64[]
            for row in eachrow(lrows)
                haskey(rmap,_pf_text(row.event_id)) || continue
                push!(le,Float64(row.error_px));push!(re,rmap[_pf_text(row.event_id)])
            end
            population=length(unique(selected.event_id))
            paired[left*"_vs_"*right]=Dict("n_common"=>length(le),"n_population"=>population,
                "common_coverage"=>population==0 ? nothing : length(le)/population,
                left=>_pf_error_summary(le),right=>_pf_error_summary(re),
                "median_right_minus_left_error_px"=>isempty(le) ? nothing : median(re-le))
        end
        all_methods=Dict{String,Dict{String,Float64}}()
        for method in ("minmax2","minmax","gcc","srp")
            all_methods[method]=Dict(_pf_text(row.event_id)=>Float64(row.error_px) for row in eachrow(selected) if row.method==method && row.projection_valid)
        end
        common=sort(collect(intersect((Set(keys(lookup)) for lookup in values(all_methods))...)))
        population=length(unique(selected.event_id))
        summaries[scope]=Dict("methods"=>methods,"paired_common_events"=>paired,
            "all_four_common_coverage"=>population==0 ? nothing : length(common)/population,
            "all_four_common_event_ids"=>common,"all_four_common_event_errors"=>Dict(method=>_pf_error_summary([lookup[id] for id in common]) for (method,lookup) in all_methods))
    end
    return summaries
end

function _srv_bindings(output_dir)
    files=["heldout_evaluation.json","calibration.json","fit_report.json","events.csv","frames.csv",
        "annotations.csv","annotations.complete.json","prepare_manifest.json"]
    all(name->isfile(joinpath(output_dir,name)),files) || throw(ArgumentError("review requires completed fit and authenticated preparation inputs"))
    hashes=Dict(name=>sha256_file(joinpath(output_dir,name)) for name in files)
    frozen=read_json(joinpath(output_dir,"heldout_evaluation.json"));report=read_json(joinpath(output_dir,"fit_report.json"))
    report["heldout_evaluation_sha256"]==hashes["heldout_evaluation.json"] || error("Frozen heldout evaluation digest mismatch")
    report["calibration_sha256"]==hashes["calibration.json"] || error("First-fit deployment calibration digest mismatch")
    for (name,key) in (("events.csv","events_sha256"),("frames.csv","frames_sha256"),("annotations.csv","annotations_sha256"),("annotations.complete.json","annotations_complete_sha256"))
        frozen["hashes"][key]==hashes[name] || error("Frozen evaluation input changed: $name")
    end
    manifest=read_json(joinpath(output_dir,"prepare_manifest.json"))
    cache=manifest["snippet_cache"];cache_path=joinpath(abspath(output_dir),"event_snippets.f64")
    abspath(cache["path"])==cache_path || error("Unexpected authenticated snippet-cache path")
    cache["bytes"]==filesize(cache_path) || error("Snippet-cache byte count mismatch")
    sha256_file(cache_path)==cache["sha256"] || error("Snippet-cache digest mismatch")
    hashes["event_snippets.f64"]=cache["sha256"]
    hashes["srp_review_source_sha256"]=sha256_file(@__FILE__)
    hashes["prepare_source_sha256"]=sha256_file(joinpath(@__DIR__,"prepare.jl"))
    hashes["physical_fit_source_sha256"]=sha256_file(joinpath(@__DIR__,"physical_fit.jl"))
    return hashes,frozen
end

"""
    review_estimators(output_dir; resume=true) -> Dict

Post-first-fit comparison on the exact frozen eligible annotated underwater,
visible, mapped, nonpurged population. Existing minmax2/original minmax/GCC
delays are read without alteration; fixed native SRP-PHAT uses authenticated
exact preparation snippets. Every projection uses the training-only frozen
evaluation_calibration, including heldout rows. No fitting, method choice or
hyperparameter tuning occurs. Writes estimator_review.csv/json ONLY; first-fit
heldout_evaluation, deployment calibration and all input hashes are preserved.
"""
function review_estimators(output_dir;resume=true)
    try
        output_dir=abspath(output_dir);hashes,frozen=_srv_bindings(output_dir)
        json_path=joinpath(output_dir,"estimator_review.json");csv_path=joinpath(output_dir,"estimator_review.csv")
        if resume && isfile(json_path) && isfile(csv_path)
            prior=read_json(json_path)
            get(prior,"hashes",nothing)==hashes && get(prior,"estimator_review_csv_sha256",nothing)==sha256_file(csv_path) && return prior
        end
        input=_pf_load(output_dir,"minmax2");events=input.events;data=input.data
        expected_train=Set(_pf_text.(frozen["split"]["train_event_ids"]));expected_heldout=Set(_pf_text.(frozen["split"]["heldout_event_ids"]))
        isempty(intersect(expected_train,expected_heldout)) || error("Frozen event appears in both outer splits")
        expected=union(expected_train,expected_heldout)
        actual=Set(_pf_text.(events.event_id[data.event_row]))
        actual==expected || error("Eligible population changed after first-fit evaluation freeze")
        length(unique(events.event_id))==nrow(events) || error("Duplicate event_id in estimator review")
        cal=frozen["evaluation_calibration"];plans=Dict{Float64,Any}();records=NamedTuple[]
        for (i,row_index) in enumerate(data.event_row)
            row=events[row_index,:];id=_pf_text(row.event_id)
            split=id in expected_heldout ? "heldout" : "train"
            input.split[row_index]==split || error("Event split changed after first-fit freeze")
            snippet=read_event_snippet(output_dir,row;verify=true)
            plan=get!(plans,Float64(row.fs)) do;_srv_plan(row.fs);end
            srp=review_srp_phat(snippet;fs=row.fs,plan=plan)
            for method in ("minmax2","minmax","gcc","srp")
                delays=method=="srp" ? srp.tdoas_samples : [row[Symbol(method*"_d"*string(j))] for j in 1:3]
                estimator_valid=method=="srp" ? srp.valid : all(_pf_finite,delays) && abs(delays[3]/row.fs)<=1e-12
                status=method=="srp" ? srp.status : estimator_valid ? "stored_preparation_estimate" : "invalid_stored_delays"
                projected=estimator_valid ? project_tdoas(cal,delays;fs=row.fs) : _pf_invalid("invalid_estimator")
                pixel_error=projected.valid ? hypot(projected.x-data.x[i],projected.y-data.y[i]) : NaN
                push!(records,(event_id=id,frame_id=_pf_text(row.frame_id),device=_pf_text(row.device),
                    audio_file=_pf_text(row.audio_file),sample_index=row.sample_index,recording_time_s=row.recording_time_s,
                    fs=row.fs,split=split,method=method,annotation_x=data.x[i],annotation_y=data.y[i],
                    d1=_pf_number(delays[1]),d2=_pf_number(delays[2]),d3=_pf_number(delays[3]),
                    estimator_valid=estimator_valid,estimator_status=status,x=projected.x,y=projected.y,
                    projection_valid=projected.valid,projection_status=projected.status,in_frame=projected.in_frame,error_px=pixel_error,
                    srp_score=method=="srp" ? srp.score : NaN,srp_pair12=method=="srp" ? srp.pair_scores[1] : NaN,
                    srp_pair13=method=="srp" ? srp.pair_scores[2] : NaN,srp_pair23=method=="srp" ? srp.pair_scores[3] : NaN,
                    srp_u=method=="srp" ? srp.uv[1] : NaN,srp_v=method=="srp" ? srp.uv[2] : NaN,
                    srp_grid_boundary=method=="srp" && srp.grid_boundary,snippet_sha256=_pf_text(row.snippet_sha256),
                    window_status=_pf_text(row.window_status),window_source_segments=_pf_text(row.window_source_segments)))
            end
        end
        isempty(records) && throw(ArgumentError("No frozen eligible events to review"))
        table=DataFrame(records);summary=_srv_summarize(table)
        # Validate every protected identity again before committing review-only
        # outputs, including cache/table/source mutations during computation.
        after,_=_srv_bindings(output_dir);after==hashes || error("Protected fit/input/source identity changed during estimator review")
        temporary=csv_path*".tmp."*string(getpid());CSV.write(temporary,table)
        open(temporary,"r+") do io;flush(io);ccall(:fsync,Cint,(Cint,),Base.fd(io));end
        mv(temporary,csv_path;force=true)
        report=Dict{String,Any}("schema_version"=>"clicker-estimator-review/1.0","status"=>"complete",
            "created_utc"=>string(Dates.now(Dates.UTC))*"Z","hashes"=>hashes,"specification"=>_srv_spec(),
            "estimator_review_csv_sha256"=>sha256_file(csv_path),"population_event_ids"=>sort(collect(expected)),
            "n_population"=>length(expected),"n_unique_frames"=>length(unique(data.frame)),"summaries"=>summary,
            "evaluation_model"=>"frozen evaluation_calibration trained without heldout labels, applied identically to all four methods",
            "heldout_interpretation"=>"Fixed post-first-fit estimator comparison, without method-specific calibration or heldout tuning. No estimator is selected for deployment from these results.",
            "quality_limits"=>["SRP score and GCC closure are internal consistency diagnostics, not source truth or sound-speed evidence.",
                "Far-field positive-front planar model; no range/front-back resolution; multipath and near-field effects remain.",
                "SRP steers at predeclared c1540 and returns raw reference-relative delays; the frozen fitted delay offsets, c and rotation then apply once through project_tdoas.",
                "Same event and annotation population is retained even when an estimator/projection is invalid; errors require accompanying coverage and paired common-event counts."],
            "recommendation"=>"Review heldout paired common-event errors together with per-method coverage. Keep the frozen first-fit calibration; recommend any method-specific calibration or broader estimator changes only as a new training-block validation experiment. No broader algorithms are implemented here.")
        write_json(json_path,report)
        return report
    catch exception
        log_failure("review_estimators $(output_dir)",exception,catch_backtrace())
        rethrow()
    end
end
