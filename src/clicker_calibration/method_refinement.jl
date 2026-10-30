# Versioned physical refinement after the frozen first-fit estimator review.
# Include after physical_fit.jl, prepare.jl and srp_review.jl. Definitions only.

function _mr_estimator_spec(method;fs=192000.0)
    method in ("gcc","srp") || throw(ArgumentError("refinement supports gcc or srp only"))
    common=Dict{String,Any}("method"=>method,"sample_rate_hz"=>Float64(fs),
        "input"=>"Float64 samples x selected physical3,4,5; already fourth-order zero-phase HP1000 filtered, matching preparation",
        "window_samples_inclusive"=>[-300,700],"reference_channel"=>3,
        "output"=>"raw reference-relative sample delays; fitted offsets and rotation applied later by project_tdoas exactly once",
        "prepare_source_sha256"=>sha256_file(joinpath(@__DIR__,"prepare.jl")),
        "srp_review_source_sha256"=>sha256_file(joinpath(@__DIR__,"srp_review.jl")))
    if method=="gcc"
        common["function"]="prep_gcc"
        common["parameters"]=Dict("c_min_m_s"=>1400.0,"reference"=>3,
            "band_hz"=>[1000.0,"0.98*fs/2"],"interp"=>8,"guide"=>"nearest envelope peak among peaks >=0.7 strongest",
            "pair_reconciliation"=>"exact preparation constrained all-pair least squares on D3 physical-lag hexagon")
    else
        common["function"]="review_srp_phat"
        common["parameters"]=_srv_spec()
    end
    return common
end

"""
    estimate_method_tdoas(calibration, data, peaks; fs, window=-300:700, ref_channel=3)

Return N×3 raw reference-relative sample delays for the calibration's method.
`data` is samples×3 for selected physical3–5, ALREADY fourth-order filtfilt
HP1000 processed exactly as preparation. This function does not filter again.
`peaks` are ordinary Julia ONE-origin row indices; `window` is an inclusive
integer offset range, clipped to the provided data just like parent TDOA APIs.
Use the same cached cross-file context if estimating a split-boundary event.
GCC calls prep_gcc with fixed c_min=1400; SRP calls the original fixed review
estimator (c1540/grid1degree/refine10/interp8). Neither uses fitted sound speed
to steer/restrict estimation. No fitted delay correction or rotation is applied
here; project_tdoas applies them once later. Five-channel/raw or HP5000 input
must be selected/preprocessed by the caller before this operation.
"""
function estimate_method_tdoas(calibration::AbstractDict,data::AbstractMatrix,peaks::AbstractVector;
        fs,window=-300:700,ref_channel=3)
    ref_channel==3 || throw(ArgumentError("selected reference3/physical5 is required"))
    size(data,2)==3 || throw(DimensionMismatch("method-matched data must contain selected physical3-5 only"))
    _pf_finite(fs) && fs>0 || throw(ArgumentError("positive finite fs required"))
    window isa AbstractRange{<:Integer} && !isempty(window) && step(window)==1 || throw(ArgumentError("window must be a nonempty contiguous integer offset range"))
    method=_pf_text(get(calibration,"method",""))
    method in ("gcc","srp") || throw(ArgumentError("calibration has no supported method-matched gcc/srp estimator"))
    specification=get(calibration,"estimator_spec",nothing)
    if specification!==nothing
        get(specification,"method",nothing)==method || error("calibration/estimator method mismatch")
        Float64(specification["sample_rate_hz"])==Float64(fs) || error("sample rate differs from fitted estimator")
        expected=_mr_estimator_spec(method;fs=fs)
        specification==expected || error("fitted estimator specification/source differs from loaded implementation")
        [first(window),last(window)]==specification["window_samples_inclusive"] || error("window differs from fitted estimator")
    end
    delays=fill(NaN,length(peaks),3)
    search=method=="srp" ? _srv_plan(fs) : nothing
    for (i,peak) in enumerate(peaks)
        peak isa Integer && 1<=peak<=size(data,1) || throw(ArgumentError("peaks must be one-origin rows within data"))
        lo=max(1,peak+first(window));hi=min(size(data,1),peak+last(window))
        lo<=hi || continue
        snippet=Matrix{Float64}(data[lo:hi,:])
        all(isfinite,snippet) || continue
        if method=="gcc"
            size(snippet,1)>=4 || continue
            delay,_=prep_gcc(snippet;fs=fs,c_min=1400.0,reference=3)
            delays[i,:]=delay
        else
            estimate=review_srp_phat(snippet;fs=fs,plan=search)
            estimate.valid && (delays[i,:]=estimate.tdoas_samples)
        end
    end
    return delays
end

function _mr_plan(plan)
    get(plan,"version","")=="clicker-method-refinement-v1" || error("unsupported refinement plan version")
    get(plan,"methods",nothing)==["gcc","srp"] || error("plan must freeze GCC then SRP tie priority")
    get(plan,"model",nothing)=="centered_pinhole" && get(plan,"parameters",nothing)==8 || error("only the frozen centered eight-parameter physical model is permitted")
    get(plan,"sound_speed_bounds_m_s",nothing)==[1400,1600] || error("sound-speed bounds differ from frozen model")
    get(plan,"nstarts",nothing)==6 && get(plan,"maxiters",nothing)==2500 && get(plan,"robust_delta_px",nothing)==12.0 || error("plan differs from frozen six-start/2500-iteration/12px settings")
    length(unique(plan["train_blocks"]))>=2 || error("blocked selection needs at least two frozen training blocks")
    return plan
end

_mr_canonical_names()=["calibration.json","fit_report.json","heldout_evaluation.json","predictions.csv",
    "estimator_review.json","estimator_review.csv","events.csv","frames.csv","annotations.csv",
    "annotations.complete.json","event_snippets.f64"]

"""Freeze the same predeclared protocol for a future device before any fit.
Existing plans (including A's user-reviewed plan) are read without modification.
Only event identity/timing metadata determine training blocks; no image errors,
heldout labels, method scores or fitted values enter automatic plan creation.
"""
function _mr_ensure_plan(output_dir)
    destination=joinpath(output_dir,"refinement_v1");path=joinpath(destination,"plan.json")
    isfile(path) && return path
    names=_mr_canonical_names()
    all(name->isfile(joinpath(output_dir,name)),names) || error("first fit and fixed estimator review must complete before freezing a refinement plan")
    frozen=read_json(joinpath(output_dir,"heldout_evaluation.json"))
    train_ids=Set(_pf_text.(frozen["split"]["train_event_ids"]))
    events=CSV.read(joinpath(output_dir,"events.csv"),DataFrame;stringtype=String)
    rows=[row for row in eachrow(events) if _pf_text(row.event_id) in train_ids]
    length(rows)==length(train_ids) || error("cannot freeze refinement plan: missing/duplicate training identities")
    all(row->_pf_finite(row.recording_time_s),rows) || error("invalid frozen training times")
    blocks=sort(unique(floor(Int,row.recording_time_s/30) for row in rows))
    length(blocks)>=2 || error("cannot freeze leave-block-out refinement with fewer than two training blocks")
    plan=Dict{String,Any}("version"=>"clicker-method-refinement-v1","created_utc"=>string(Dates.now(Dates.UTC))*"Z",
        "reason"=>"Fixed post-first-fit GCC/SRP physical refinement; automatic plan frozen before fitting",
        "canonical_input_sha256"=>Dict(name=>sha256_file(joinpath(output_dir,name)) for name in names),
        "methods"=>["gcc","srp"],"model"=>"centered_pinhole","parameters"=>8,
        "sound_speed_bounds_m_s"=>[1400,1600],"fixed_geometry"=>"D3 side0.18m physical3-5 ref5",
        "nstarts"=>6,"maxiters"=>2500,"robust_delta_px"=>12.0,
        "training_event_source"=>"heldout_evaluation.json split.train_event_ids","train_blocks"=>blocks,
        "selection"=>"leave-one-training-block-out; minimize equal-block mean of existing invalid-penalized RMSE; deterministic tie priority GCC then SRP; no heldout-label selection",
        "fold_limits"=>"Report small training folds and parameter/convergence/identifiability instability; never change criteria from reused evaluation",
        "evaluation"=>"Freeze selection and full-training fit before descriptive replay of already-inspected outer blocks; not untouched confirmation",
        "deployment"=>"Refit selected method to all frozen eligible annotations after evaluation; preserve canonical first-fit artifacts",
        "automatic_defaults"=>"No defaults installed by refinement; inspect errors, coverage and physical uncertainty first")
    _mr_plan(plan)
    mkpath(destination)
    isfile(path) && return path
    write_json(path,plan)
    return path
end

function _mr_bindings(output_dir)
    destination=joinpath(output_dir,"refinement_v1");plan_path=_mr_ensure_plan(output_dir)
    plan=_mr_plan(read_json(plan_path))
    allowed=Set(_mr_canonical_names())
    Set(keys(plan["canonical_input_sha256"]))==allowed || error("plan must bind every canonical first-fit/review/input/cache artifact")
    hashes=Dict{String,Any}("plan_sha256"=>sha256_file(plan_path))
    for name in sort(collect(allowed))
        path=joinpath(output_dir,name)
        isfile(path) && sha256_file(path)==plan["canonical_input_sha256"][name] || error("frozen canonical input changed: $name")
        hashes[name]=plan["canonical_input_sha256"][name]
    end
    for name in ("physical_fit.jl","prepare.jl","srp_review.jl","method_refinement.jl","support.jl")
        hashes["source/"*name]=sha256_file(joinpath(@__DIR__,name))
    end
    firstfit=read_json(joinpath(output_dir,"heldout_evaluation.json"))
    baseline_report=read_json(joinpath(output_dir,"fit_report.json"))
    baseline_report["calibration_sha256"]==hashes["calibration.json"] &&
        baseline_report["heldout_evaluation_sha256"]==hashes["heldout_evaluation.json"] || error("canonical first-fit bindings are inconsistent")
    review=read_json(joinpath(output_dir,"estimator_review.json"))
    review["estimator_review_csv_sha256"]==hashes["estimator_review.csv"] || error("review table digest mismatch")
    for name in ("calibration.json","heldout_evaluation.json","events.csv","frames.csv","annotations.csv","annotations.complete.json","event_snippets.f64")
        review["hashes"][name]==hashes[name] || error("review is not bound to frozen input $name")
    end
    for name in ("physical_fit.jl","prepare.jl","srp_review.jl")
        key=replace(name,".jl"=>"_source_sha256")
        review["hashes"][key]==hashes["source/"*name] || error("review estimator/projector source changed: $name")
    end
    return (destination=destination,plan=plan,hashes=hashes,firstfit=firstfit)
end

function _mr_load(output_dir,binding)
    input=_pf_load(output_dir,"minmax2");events=input.events;data=input.data
    split=binding.firstfit["split"]
    train_ids=Set(_pf_text.(split["train_event_ids"]));heldout_ids=Set(_pf_text.(split["heldout_event_ids"]))
    isempty(intersect(train_ids,heldout_ids)) || error("frozen outer event IDs overlap")
    expected=union(train_ids,heldout_ids)
    length(unique(events.event_id))==nrow(events) || error("duplicate event IDs")
    Set(_pf_text.(events.event_id[data.event_row]))==expected || error("eligible population changed after first-fit freeze")
    train=findall(i->_pf_text(events.event_id[data.event_row[i]]) in train_ids,eachindex(data.x))
    heldout=findall(i->_pf_text(events.event_id[data.event_row[i]]) in heldout_ids,eachindex(data.x))
    sort(unique(data.block[train]))==sort(Int.(binding.plan["train_blocks"])) || error("training blocks differ from frozen plan")
    for i in eachindex(data.x)
        id=_pf_text(events.event_id[data.event_row[i]])
        input.split[data.event_row[i]]==(id in heldout_ids ? "heldout" : "train") || error("outer split identity changed")
    end
    review=CSV.read(joinpath(output_dir,"estimator_review.csv"),DataFrame;stringtype=String)
    all(in(names(review)),["event_id","frame_id","method","fs","split","d1","d2","d3","estimator_valid","snippet_sha256"]) || error("review table lacks method/identity fields")
    methods=Dict{String,Any}()
    for method in ("gcc","srp")
        selected=review[review.method.==method,:]
        nrow(selected)==length(expected) && length(unique(selected.event_id))==length(expected) && Set(_pf_text.(selected.event_id))==expected || error("method $method does not cover the same frozen event population")
        lookup=Dict(_pf_text(row.event_id)=>row for row in eachrow(selected))
        tau=fill(NaN,length(data.x),2)
        for (i,index) in enumerate(data.event_row)
            event=events[index,:];row=lookup[_pf_text(event.event_id)]
            _pf_text(row.frame_id)==_pf_text(event.frame_id) && row.fs==event.fs &&
                _pf_text(row.split)==input.split[index] && _pf_text(row.snippet_sha256)==_pf_text(event.snippet_sha256) || error("method review identity mismatch for $(event.event_id)")
            delays=[row.d1,row.d2,row.d3]
            if row.estimator_valid && all(_pf_finite,delays) && abs(Float64(delays[3])/event.fs)<=1e-12
                tau[i,:]=Float64.(delays[1:2])./Float64(event.fs)
            end
        end
        methods[method]=merge(data,(tau=tau,))
    end
    return (input=input,methods=methods,train=train,heldout=heldout,
        train_ids=train_ids,heldout_ids=heldout_ids)
end

function _mr_convergence(fit)
    matching=[row for row in fit.starts if row["q"]==fit.q]
    best=isempty(matching) ? Dict{String,Any}() : last(matching)
    return Dict("best_start"=>get(best,"start",nothing),"best_converged"=>get(best,"converged",false),
        "best_iterations"=>get(best,"iterations",nothing),
        "converged_runs"=>count(row->row["converged"],fit.starts),"n_runs"=>length(fit.starts))
end

function _mr_calibration(fit,input,method)
    calibration=_pf_calibration(fit.q,input.width,input.height,"centered_pinhole";device=input.device,method=method)
    calibration["estimator_spec"]=_mr_estimator_spec(method;fs=Float64(input.events.fs[first(input.data.event_row)]))
    calibration["parameter_bound_hits_normalized"]=findall(abs.(fit.q).>=0.98)
    calibration["optimization_convergence"]=_mr_convergence(fit)
    return calibration
end

function _mr_fold_instability(folds)
    successful=[fold for fold in folds if get(fold,"status","")=="fitted"]
    isempty(successful) && return Dict("status"=>"no_successful_fold")
    c=[fold["calibration"]["sound_speed_m_s"] for fold in successful]
    fov=reduce(vcat,[permutedims(Float64.(fold["calibration"]["fov_degrees"])) for fold in successful])
    pose=0.0
    for left in successful,right in successful
        A=_pf_matrix(left["calibration"]["rotation_array_to_camera"])
        B=_pf_matrix(right["calibration"]["rotation_array_to_camera"])
        pose=max(pose,rad2deg(acos(clamp((tr(A'B)-1)/2,-1.0,1.0))))
    end
    return Dict("sound_speed_range_m_s"=>[minimum(c),maximum(c)],
        "fov_min_degrees"=>vec(minimum(fov;dims=1)),"fov_max_degrees"=>vec(maximum(fov;dims=1)),
        "max_pairwise_rotation_difference_degrees"=>pose,
        "minimum_fold_training_events"=>minimum(fold["n_training"] for fold in folds),
        "n_nonconverged_best"=>count(fold->!fold["convergence"]["best_converged"],successful),
        "n_bound_hit_folds"=>count(fold->!isempty(fold["calibration"]["parameter_bound_hits_normalized"]),successful),
        "interpretation"=>"Fold parameter variation is a stability diagnostic, not an uncertainty confidence interval; small training folds and missing angular/elevation support can make physical parameters nonunique.")
end

# Receives training data ONLY; the signature has no outer-evaluation labels.
function _mr_select(training,width,height,plan,device)
    evidence=Dict{String,Any}();scores=Float64[]
    for method in ("gcc","srp")
        data=training[method];folds=Any[];fold_scores=Float64[]
        for block in sort(Int.(plan["train_blocks"]))
            validation_indices=findall(==(block),data.block)
            validation_frames=Set(data.frame[validation_indices])
            training_indices=findall(i->data.block[i]!=block && !(data.frame[i] in validation_frames),eachindex(data.x))
            train=_pf_subset(data,training_indices);validation=_pf_subset(data,validation_indices)
            entry=Dict{String,Any}("validation_block"=>block,"training_blocks"=>sort(unique(train.block)),
                "n_training"=>length(train.x),"n_validation"=>length(validation.x),
                "training_event_rows"=>train.event_row,"validation_event_rows"=>validation.event_row,
                "small_fold_warning"=>length(train.x)<=9 ? "At most nine training points for eight parameters; interpret validation and parameter stability cautiously" : nothing)
            try
                fit=_pf_fit(train,width,height,"centered_pinhole";nstarts=plan["nstarts"],maxiters=plan["maxiters"],delta=plan["robust_delta_px"])
                cal=_pf_calibration(fit.q,width,height,"centered_pinhole";device=device,method=method)
                score=_pf_validation_score(cal,validation,width,height)
                merge!(entry,Dict("status"=>"fitted","calibration"=>cal,"q"=>fit.q,
                    "training_objective"=>fit.objective,"penalized_validation_rmse_px"=>score,
                    "validation_metrics"=>_pf_metrics(cal,validation),"training_input_screen"=>_pf_screen_summary(train),
                    "convergence"=>_mr_convergence(fit),"optimization"=>fit.starts,
                    "diagnostics_training_only"=>_pf_diagnostics(fit.q,train,width,height,"centered_pinhole";profiles=false,delta=plan["robust_delta_px"])))
                cal["parameter_bound_hits_normalized"]=findall(abs.(fit.q).>=0.98)
                push!(fold_scores,score)
            catch exception
                log_failure("method refinement $method leave-block-$block",exception,catch_backtrace())
                merge!(entry,Dict("status"=>"fit_failed","error"=>sprint(showerror,exception),
                    "penalized_validation_rmse_px"=>Inf));push!(fold_scores,Inf)
            end
            push!(folds,entry)
        end
        aggregate=mean(fold_scores);push!(scores,aggregate)
        evidence[method]=Dict("equal_block_mean_penalized_rmse_px"=>aggregate,"folds"=>folds,
            "fold_instability"=>_mr_fold_instability(folds))
    end
    all(value->!isfinite(value),scores) && error("all predeclared methods failed blocked selection")
    # Strict improvement only: an exact tie deterministically retains GCC.
    selected=scores[2]<scores[1] ? "srp" : "gcc"
    return selected,evidence
end

function _mr_predictions(loaded,evaluation,calibration,method)
    data=loaded.methods[method];events=loaded.input.events
    rows=NamedTuple[]
    for (i,index) in enumerate(data.event_row)
        event=events[index,:];delays=[data.tau[i,1],data.tau[i,2],0.0]
        e=project_tdoas(evaluation,delays;units="seconds");d=project_tdoas(calibration,delays;units="seconds")
        baseline=_pf_d3_baseline(delays,loaded.input.width,loaded.input.height;units="seconds")
        deploy_error=d.valid ? hypot(d.x-data.x[i],d.y-data.y[i]) : NaN
        push!(rows,(event_id=_pf_text(event.event_id),frame_id=_pf_text(event.frame_id),device=loaded.input.device,
            method=method,split=loaded.input.split[index],annotation_x=data.x[i],annotation_y=data.y[i],
            d1_samples=data.tau[i,1]*event.fs,d2_samples=data.tau[i,2]*event.fs,d3_samples=0.0,
            evaluation_x=e.x,evaluation_y=e.y,evaluation_valid=e.valid,evaluation_status=e.status,evaluation_in_frame=e.in_frame,
            evaluation_error_px=e.valid ? hypot(e.x-data.x[i],e.y-data.y[i]) : NaN,
            predicted_x=d.x,predicted_y=d.y,predicted_valid=d.valid,predicted_status=d.status,predicted_in_frame=d.in_frame,
            predicted_error_px=deploy_error,final_x=d.x,final_y=d.y,final_valid=d.valid,final_status=d.status,final_in_frame=d.in_frame,final_error_px=deploy_error,
            baseline_x=baseline.x,baseline_y=baseline.y,baseline_valid=baseline.valid,baseline_status=baseline.status,baseline_in_frame=baseline.in_frame,
            baseline_error_px=baseline.valid ? hypot(baseline.x-data.x[i],baseline.y-data.y[i]) : NaN))
    end
    return DataFrame(rows)
end

"""
    refine_calibration(output_dir;resume=true) -> Dict

Read the predeclared output_dir/refinement_v1/plan.json, preserving every
canonical first-fit/input/review/cache artifact. Compare method-specific GCC
and fixed native SRP centered eight-parameter physical fits using ONLY the
frozen training IDs. Selection is leave-one-training-block-out equal-block mean
existing invalid-penalized RMSE, with deterministic GCC-first exact ties.
Persist selection AND the full-training fit before replaying the already-seen
outer block. Persist that reused evaluation before the selected-method all-ID
deployment refit. Writes only refinement_v1/{selection,evaluation,calibration,
predictions,report}.json/csv; never installs a default. No empirical pixel map.
If no plan exists for a future device, freeze these same fixed settings and
input hashes/metadata-derived training blocks before any fitting takes place.
"""
function refine_calibration(output_dir;resume=true)
    try
        output_dir=abspath(output_dir);binding=_mr_bindings(output_dir);destination=binding.destination
        paths=Dict(name=>joinpath(destination,name) for name in ("selection.json","evaluation.json","calibration.json","predictions.csv","report.json"))
        if isfile(paths["selection.json"])
            previous=read_json(paths["selection.json"])
            get(previous,"hashes",nothing)==binding.hashes || error("refinement source/plan/input changed; preserve this version and predeclare a new one")
        end
        if resume && all(isfile,values(paths))
            previous=read_json(paths["report.json"])
            previous["hashes"]==binding.hashes || error("refinement report binding changed")
            for (name,digest) in previous["outputs_sha256"]
                sha256_file(paths[name])==digest || error("refinement output integrity failure: $name")
            end
            return previous
        end
        loaded=_mr_load(output_dir,binding);input=loaded.input;plan=binding.plan
        training=Dict(method=>_pf_subset(loaded.methods[method],loaded.train) for method in ("gcc","srp"))
        method,candidates=_mr_select(training,input.width,input.height,plan,input.device)
        train=training[method]
        trainfit=_pf_fit(train,input.width,input.height,"centered_pinhole";nstarts=plan["nstarts"],maxiters=plan["maxiters"],delta=plan["robust_delta_px"])
        traincal=_mr_calibration(trainfit,input,method)
        diagnostics=_pf_diagnostics(trainfit.q,train,input.width,input.height,"centered_pinhole";
            profiles=true,maxiters=min(plan["maxiters"],1500),delta=plan["robust_delta_px"])
        selection=Dict{String,Any}("schema_version"=>"clicker-refinement-selection/1.0","created_utc"=>string(Dates.now(Dates.UTC))*"Z",
            "hashes"=>binding.hashes,"plan"=>plan,"selected_method"=>method,"uses_heldout_labels"=>false,
            "train_event_ids"=>collect(input.events.event_id[train.event_row]),"train_blocks"=>sort(unique(train.block)),
            "rule"=>"equal-block mean existing invalid-penalized RMSE; strict improvement required to replace GCC",
            "candidates"=>candidates,"full_training_calibration"=>traincal,"full_training_q"=>trainfit.q,
            "full_training_objective"=>trainfit.objective,"full_training_metrics"=>_pf_metrics(traincal,train),
            "full_training_optimization"=>trainfit.starts,"diagnostics_training_only"=>diagnostics,
            "sequence"=>"This selection and full-training model are persisted before any reused outer-evaluation metric is computed")
        _mr_bindings(output_dir).hashes==binding.hashes || error("protected identities changed during training selection")
        write_json(paths["selection.json"],selection);selection_sha=sha256_file(paths["selection.json"])
        # Only now may the already-inspected outer labels enter evaluation.
        heldout=_pf_subset(loaded.methods[method],loaded.heldout)
        evaluation=Dict{String,Any}("schema_version"=>"clicker-refinement-evaluation/1.0","created_utc"=>string(Dates.now(Dates.UTC))*"Z",
            "hashes"=>binding.hashes,"selection_sha256"=>selection_sha,"selected_method"=>method,
            "evaluation_calibration"=>traincal,"heldout_event_ids"=>collect(input.events.event_id[heldout.event_row]),
            "heldout_blocks"=>sort(unique(heldout.block)),"reused_evaluation_metrics"=>_pf_metrics(traincal,heldout),
            "paired_uncalibrated_same_method_baseline"=>_pf_paired_baseline(heldout,traincal,input.width,input.height),
            "interpretation"=>"Reused evaluation after inspected first-fit/review outcomes; not untouched confirmatory validation. No heldout labels entered candidate fits, selection, or training diagnostics.")
        write_json(paths["evaluation.json"],evaluation);evaluation_sha=sha256_file(paths["evaluation.json"])
        # Deployment training may include all frozen eligible labels only after
        # selection and reused evaluation have both been durably frozen.
        all_data=loaded.methods[method]
        finalfit=_pf_fit(all_data,input.width,input.height,"centered_pinhole";nstarts=plan["nstarts"],maxiters=plan["maxiters"],delta=plan["robust_delta_px"],initial=trainfit.q)
        calibration=_mr_calibration(finalfit,input,method)
        final_diagnostics=_pf_diagnostics(finalfit.q,all_data,input.width,input.height,"centered_pinhole";profiles=false,delta=plan["robust_delta_px"])
        calibration["hashes"]=binding.hashes;calibration["selection_sha256"]=selection_sha
        calibration["evaluation_sha256"]=evaluation_sha;calibration["refinement_plan_sha256"]=binding.hashes["plan_sha256"]
        calibration["fit_status"]=(get(diagnostics,"weak_identifiability",true) || get(final_diagnostics,"weak_identifiability",true)) ? "fitted_weak_identifiability" : "fitted"
        calibration["deployment_fit_scope"]="selected method, all frozen eligible annotations after frozen training selection and reused evaluation; in-sample deployment errors"
        write_json(paths["calibration.json"],calibration)
        predictions=_mr_predictions(loaded,traincal,calibration,method)
        temporary=paths["predictions.csv"]*".tmp."*string(getpid());CSV.write(temporary,predictions)
        open(temporary,"r+") do io;ccall(:fsync,Cint,(Cint,),Base.fd(io));end
        mv(temporary,paths["predictions.csv"];force=true)
        _mr_bindings(output_dir).hashes==binding.hashes || error("protected identities changed during refinement")
        sha256_file(paths["selection.json"])==selection_sha && sha256_file(paths["evaluation.json"])==evaluation_sha || error("frozen refinement selection/evaluation changed")
        outputs=Dict(name=>sha256_file(path) for (name,path) in paths if name!="report.json")
        report=Dict{String,Any}("schema_version"=>"clicker-method-refinement-report/1.0","status"=>"complete",
            "hashes"=>binding.hashes,"outputs_sha256"=>outputs,"selected_method"=>method,"fit_status"=>calibration["fit_status"],
            "selection_rule"=>selection["rule"],"selection_candidates"=>candidates,
            "training_metrics"=>selection["full_training_metrics"],"reused_evaluation_metrics"=>evaluation["reused_evaluation_metrics"],
            "deployment_metrics_descriptive"=>_pf_metrics(calibration,all_data),
            "uncalibrated_same_method_baseline_spec"=>merge(_pf_baseline_spec(input.width,input.height),Dict("method"=>method)),
            "paired_uncalibrated_baseline_training"=>_pf_paired_baseline(train,traincal,input.width,input.height),
            "paired_uncalibrated_baseline_reused_evaluation"=>evaluation["paired_uncalibrated_same_method_baseline"],
            "paired_uncalibrated_baseline_deployment_descriptive"=>_pf_paired_baseline(all_data,calibration,input.width,input.height),
            "diagnostics_training_only"=>diagnostics,"diagnostics_deployment_descriptive"=>final_diagnostics,
            "optimization_training"=>_mr_convergence(trainfit),"optimization_deployment"=>_mr_convergence(finalfit),
            "n_training"=>length(train.x),"n_reused_evaluation"=>length(heldout.x),"n_deployment"=>length(all_data.x),
            "inference_contract"=>calibration["estimator_spec"],"automatic_defaults_installed"=>false,
            "limits"=>["Eight physical parameters can be weakly identifiable with sparse vertical/directional support; inspect bound hits, speed profiles and leave-block-out pose/FOV variation.",
                "Some folds may train on only nine points; optimizer nonconvergence and fold instability are retained, not hidden or repaired using outer labels.",
                "The outer block has already been inspected. These are reused descriptive evaluation results, not fresh confirmatory accuracy.",
                "No fitted defaults are installed. Far-field front-branch assumptions remain; delay offsets/rotation must be applied once after method-matched raw estimation."])
        write_json(paths["report.json"],report)
        return report
    catch exception
        log_failure("refine_calibration $(output_dir)",exception,catch_backtrace())
        rethrow()
    end
end
