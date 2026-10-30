# Same-session projection of every prepared candidate. No fitting or import-time I/O.

function _ep_track!(inputs, path; expected=nothing)
    path=abspath(path)
    isfile(path) || error("Missing all-event projection input: $path")
    digest=sha256_file(path)
    expected===nothing || digest==expected || error("All-event projection input digest mismatch: $path")
    inputs[path]=digest
    return digest
end

_ep_binding(hashes,name,alternate)=get(hashes,name,get(hashes,alternate,nothing))

function _ep_table(path,key)
    table=CSV.read(path,DataFrame;stringtype=String)
    key in names(table) || error("Missing $key in $path")
    ids=_pf_text.(table[!,key])
    all(!isempty,ids) && length(unique(ids))==length(ids) || error("Empty/duplicate $key in $path")
    return table,Dict(id=>i for (i,id) in enumerate(ids))
end

function _ep_check_cohort(table,events,event_rows)
    "frame_id" in names(table) || error("Frozen prediction table lacks frame_id")
    for row in eachrow(table)
        id=_pf_text(row.event_id)
        haskey(event_rows,id) || error("Frozen prediction event is outside the original session: $id")
        _pf_text(row.frame_id)==_pf_text(events[event_rows[id],:frame_id]) || error("Frozen prediction frame differs: $id")
    end
end

function _ep_inputs(original_dir,model_dir;python=nothing)
    inputs=Dict{String,String}();hashes=Dict{String,Any}()
    for (name,key) in (("events.csv","events_sha256"),("frames.csv","frames_sha256"),
            ("annotations.csv","annotations_sha256"),("annotations.complete.json","annotations_complete_sha256"),
            ("prepare_manifest.json","prepare_manifest_sha256"))
        hashes[key]=_ep_track!(inputs,joinpath(original_dir,name))
    end
    seal=read_json(joinpath(original_dir,"annotations.complete.json"))
    for key in ("frames_sha256","annotations_sha256")
        get(seal,key,nothing)==hashes[key] || error("Original annotation seal is stale: $key")
    end
    prepared=read_json(joinpath(original_dir,"prepare_manifest.json"))
    outputs=Dict(abspath(row["path"])=>row["sha256"] for row in prepared["outputs"])
    for name in ("events.csv","frames.csv")
        path=joinpath(original_dir,name)
        get(outputs,path,nothing)==inputs[path] || error("Original preparation does not authenticate $name")
    end
    calibration_path=joinpath(model_dir,"calibration.json")
    hashes["calibration_sha256"]=_ep_track!(inputs,calibration_path)
    hashes["cohort_predictions_sha256"]=_ep_track!(inputs,joinpath(model_dir,"predictions.csv"))
    calibration=load_calibration(calibration_path)
    binding=calibration["hashes"]
    refined=isfile(joinpath(model_dir,"report.json"))
    fit_dir=refined ? dirname(model_dir) : model_dir
    report_name=refined ? "report.json" : "fit_report.json"
    hashes["model_report_sha256"]=_ep_track!(inputs,joinpath(model_dir,report_name))
    report=read_json(joinpath(model_dir,report_name))
    if refined
        for name in ("calibration.json","predictions.csv","evaluation.json","selection.json")
            expected=get(report["outputs_sha256"],name,nothing)
            expected===nothing && error("Refinement report lacks $name binding")
            _ep_track!(inputs,joinpath(model_dir,name);expected=expected)
        end
        report["hashes"]==binding || error("Refinement report/calibration provenance differs")
        evaluation=read_json(joinpath(model_dir,"evaluation.json"))
        calibration["evaluation_sha256"]==inputs[joinpath(model_dir,"evaluation.json")] || error("Evaluation calibration binding differs")
        evaluation_calibration=evaluation["evaluation_calibration"]
        hashes["baseline_predictions_sha256"]=_ep_track!(inputs,joinpath(fit_dir,"predictions.csv");expected=binding["predictions.csv"])
    else
        report["calibration_sha256"]==hashes["calibration_sha256"] || error("First-fit calibration digest differs")
        report["predictions_sha256"]==hashes["cohort_predictions_sha256"] || error("First-fit prediction digest differs")
        _ep_track!(inputs,joinpath(model_dir,"heldout_evaluation.json");expected=report["heldout_evaluation_sha256"])
        evaluation_calibration=read_json(joinpath(model_dir,"heldout_evaluation.json"))["evaluation_calibration"]
    end
    evaluation_calibration["method"]==calibration["method"] || error("Evaluation and deployment estimators differ")
    for (name,key) in (("events.csv","events_sha256"),("frames.csv","frames_sha256"))
        _ep_binding(binding,name,key)==hashes[key] || error("Model belongs to different original $name")
    end
    model_ann=_ep_binding(binding,"annotations.csv","annotations_sha256")
    model_seal=_ep_binding(binding,"annotations.complete.json","annotations_complete_sha256")
    model_ann===nothing && error("Model lacks annotation binding")
    model_seal===nothing && error("Model lacks annotation completion binding")
    bridge=nothing
    if model_ann!=hashes["annotations_sha256"] || model_seal!=hashes["annotations_complete_sha256"]
        helper=joinpath(@__DIR__,"fit_annotation_policy.py")
        exe=isnothing(python) ? get(ENV,"CLICKER_PYTHON","python3") : string(python)
        bridge=JSON.parse(String(run_logged(`$exe $helper validate $original_dir $fit_dir`;
            operation="validate original annotation / frozen fit-input bridge",capture=true)))
        checks=Dict("valid"=>true,"original_annotations_sha256"=>hashes["annotations_sha256"],
            "original_completion_sha256"=>hashes["annotations_complete_sha256"],"frames_sha256"=>hashes["frames_sha256"],
            "fit_annotations_sha256"=>model_ann,"fit_completion_sha256"=>model_seal)
        all(get(bridge,k,nothing)==v for (k,v) in checks) || error("Validated annotation bridge does not bind this model")
        for (path,key) in ((helper,"helper_sha256"),(joinpath(@__DIR__,"annotation_revision.py"),"clone_helper_sha256"),
                (joinpath(fit_dir,"annotation_fit_policy.json"),"policy_sha256"),
                (joinpath(fit_dir,"annotations.csv"),"fit_annotations_sha256"),
                (joinpath(fit_dir,"annotations.complete.json"),"fit_completion_sha256"),
                (joinpath(original_dir,"human_review_selection.json"),"selection_sha256"),
                (joinpath(fit_dir,"prepare_manifest.json"),"fit_prepare_manifest_sha256"))
            _ep_track!(inputs,path;expected=bridge[key])
        end
    end
    # Native B policy is provenance only. It does not admit/select projection rows.
    first_frozen=joinpath(fit_dir,"heldout_evaluation.json")
    if isfile(first_frozen)
        _ep_track!(inputs,first_frozen;expected=get(binding,"heldout_evaluation.json",nothing))
        firstfit=read_json(first_frozen)
        policy_sha=get(firstfit["hashes"],"fit_observation_policy_sha256",nothing)
        policy_sha===nothing || _ep_track!(inputs,joinpath(fit_dir,"fit_observation_policy.json");expected=policy_sha)
    end
    events,event_rows=_ep_table(joinpath(original_dir,"events.csv"),"event_id")
    frames,frame_rows=_ep_table(joinpath(original_dir,"frames.csv"),"frame_id")
    annotations,annotation_rows=_ep_table(joinpath(original_dir,"annotations.csv"),"frame_id")
    Set(keys(frame_rows))==Set(keys(annotation_rows)) || error("Original annotations do not cover original frames exactly")
    required=["frame_id","device","fs","recording_time_s","sync_status"]
    all(in(names(events)),required) || error("Original events lack required identity/timing fields")
    all(_pf_text(d)==_pf_text(calibration["device"]) for d in events.device) || error("Device differs from model")
    for event in eachrow(events)
        fid=_pf_text(event.frame_id)
        isempty(fid) && continue
        haskey(frame_rows,fid) || error("Event refers to unknown original frame $fid")
        frame=frames[frame_rows[fid],:]
        frame.width==calibration["width"] && frame.height==calibration["height"] || error("Original frame resolution differs from model")
    end
    cohort,cohort_rows=_ep_table(joinpath(model_dir,"predictions.csv"),"event_id")
    _ep_check_cohort(cohort,events,event_rows)
    baseline,baseline_rows=refined ? _ep_table(joinpath(fit_dir,"predictions.csv"),"event_id") : (cohort,cohort_rows)
    _ep_check_cohort(baseline,events,event_rows)
    method=calibration["method"]
    if method in ("gcc","srp")
        spec=get(calibration,"estimator_spec",nothing)
        spec===nothing && error("Method-specific calibration lacks estimator specification")
        spec==_mr_estimator_spec(method;fs=spec["sample_rate_hz"]) || error("Frozen estimator/source specification differs")
        all(!_pf_finite(fs) || fs==spec["sample_rate_hz"] for fs in events.fs) || error("Original sample rate differs from estimator")
        get(prepared["request"],"source_sha256",nothing)==spec["prepare_source_sha256"] || error("Original delay producer differs from estimator")
    end
    saved_srp=Dict{String,Any}();cache_dir=nothing
    if method=="srp"
        review_csv=joinpath(fit_dir,"estimator_review.csv");review_json=joinpath(fit_dir,"estimator_review.json")
        _ep_track!(inputs,review_json;expected=binding["estimator_review.json"])
        review=read_json(review_json)
        _ep_track!(inputs,review_csv;expected=binding["estimator_review.csv"])
        review["estimator_review_csv_sha256"]==inputs[review_csv] || error("Saved SRP review CSV binding differs")
        for (name,key) in (("events.csv","events_sha256"),("frames.csv","frames_sha256"))
            review["hashes"][name]==hashes[key] || error("Saved SRP review belongs to different $name")
        end
        for row in eachrow(CSV.read(review_csv,DataFrame;stringtype=String))
            _pf_text(row.method)=="srp" || continue
            id=_pf_text(row.event_id)
            haskey(saved_srp,id) && error("Duplicate saved SRP event $id")
            haskey(event_rows,id) || error("Saved SRP event outside original session")
            event=events[event_rows[id],:]
            _pf_text(row.frame_id)==_pf_text(event.frame_id) && row.fs==event.fs &&
                _pf_text(row.snippet_sha256)==_pf_text(event.snippet_sha256) || error("Saved SRP snippet/event identity differs")
            saved_srp[id]=row
        end
        cache_path=joinpath(original_dir,"event_snippets.f64")
        if isfile(cache_path)
            cache=prepared["snippet_cache"]
            cache["dtype"]=="Float64" && cache["endianness"]=="little" && cache["layout"]=="column-major" || error("Unsupported original cache encoding")
            cache["bytes"]==filesize(cache_path) || error("Original snippet cache size differs")
            digest=_ep_track!(inputs,cache_path;expected=cache["sha256"])
            get(outputs,cache_path,nothing)==digest || error("Original cache producer binding differs")
            get(binding,"event_snippets.f64",digest)==digest || error("Original cache differs from fitted estimator cache")
            hashes["event_snippets_sha256"]=digest;cache_dir=original_dir
        end
    else
        all(in(names(events)),[method*"_d$j" for j in 1:3]) || error("Original events lack saved $method delays")
    end
    source_hashes=Dict(name=>sha256_file(joinpath(@__DIR__,name)) for name in
        ("event_projection.jl","physical_fit.jl","prepare.jl","srp_review.jl","method_refinement.jl","support.jl","workflow.jl"))
    return (;inputs,hashes,source_hashes,calibration,evaluation_calibration,events,annotations,annotation_rows,
        cohort,cohort_rows,baseline,baseline_rows,saved_srp,cache_dir,bridge)
end

"""Acoustic-only phase. Annotation data cannot enter this function."""
function _ep_project(events,calibration,evaluation_calibration,saved_srp;cache_dir=nothing)
    method=calibration["method"];n=nrow(events);table=copy(events)
    plans=Dict{Float64,Any}()
    table.method=fill(method,n)
    table.d1_samples=fill(NaN,n);table.d2_samples=fill(NaN,n);table.d3_samples=fill(NaN,n)
    table.estimator_status=fill("",n)
    table.td_source=fill("stored_events",n);table.estimator_error=fill("",n)
    for prefix in ("final","evaluation","baseline")
        table[!,prefix*"_x"]=fill(NaN,n);table[!,prefix*"_y"]=fill(NaN,n)
        table[!,prefix*"_valid"]=falses(n);table[!,prefix*"_status"]=fill("",n);table[!,prefix*"_in_frame"]=falses(n)
    end
    for (i,event) in enumerate(eachrow(events))
        delays=fill(NaN,3);estimate_status="stored_preparation_estimate"
        if method=="srp"
            saved=get(saved_srp,_pf_text(event.event_id),nothing)
            if saved===nothing
                table.td_source[i]="unavailable"
                estimate_status="missing_cached_snippet"
                if cache_dir!==nothing
                    try
                        # read_event_snippet caches whole-file authentication by
                        # inode/size/mtime/ctime; every payload remains verified.
                        snippet=read_event_snippet(cache_dir,event;verify=true)
                        plan=get!(plans,Float64(event.fs)) do;_srv_plan(event.fs);end
                        estimate=review_srp_phat(snippet;fs=event.fs,plan=plan)
                        delays=estimate.tdoas_samples
                        estimate_status=estimate.valid ? "fixed_srp_ok" : "fixed_srp_"*estimate.status
                        table.td_source[i]="fixed_estimator_cached_snippet"
                    catch exception
                        (exception isa InterruptException || exception isa OutOfMemoryError) && rethrow()
                        estimate_status="invalid_cached_snippet"
                        table.estimator_error[i]=sprint(showerror,exception)
                        log_failure("all-event fixed SRP inference $(_pf_text(event.event_id))",exception,catch_backtrace())
                    end
                end
            elseif !saved.estimator_valid
                estimate_status="invalid_saved_srp"
                table.td_source[i]="cached_review"
            else
                delays=_pf_number.([saved.d1,saved.d2,saved.d3]);estimate_status="saved_exact_cache_srp"
                table.td_source[i]="cached_review"
            end
        else
            delays=_pf_number.([event[Symbol(method*"_d$j")] for j in 1:3])
        end
        table[i,[:d1_samples,:d2_samples,:d3_samples]]=delays
        table.estimator_status[i]=estimate_status
        for (prefix,cal) in (("final",calibration),("evaluation",evaluation_calibration),("baseline",calibration))
            result=if estimate_status in ("missing_cached_snippet","invalid_cached_snippet","invalid_saved_srp") ||
                    (startswith(estimate_status,"fixed_srp_") && estimate_status!="fixed_srp_ok")
                _pf_invalid(estimate_status)
            elseif !_pf_finite(event.fs) || event.fs<=0
                _pf_invalid("invalid_fs")
            elseif prefix=="baseline"
                _pf_d3_baseline(delays,cal["width"],cal["height"];fs=event.fs,units="samples")
            else
                project_tdoas(cal,delays;fs=event.fs)
            end
            for field in ("x","y","valid","status","in_frame")
                table[i,prefix*"_"*field]=getproperty(result,Symbol(field))
            end
        end
    end
    for field in ("x","y","valid","status","in_frame")
        table[!,"predicted_"*field]=copy(table[!,"final_"*field])
    end
    return table
end

"""Join original visual truth and frozen cohort metadata after acoustic projection."""
function _ep_join_truth!(table,input)
    n=nrow(table)
    table.annotation_x=fill(NaN,n);table.annotation_y=fill(NaN,n);table.uncertainty_px=fill(NaN,n)
    table.visibility=fill("",n);table.medium=fill("",n)
    table.split=fill("not_in_frozen_cohort",n);table.eligibility_status=fill("not_in_frozen_cohort",n)
    for (i,event) in enumerate(eachrow(table))
        fid=_pf_text(event.frame_id);id=_pf_text(event.event_id)
        if haskey(input.annotation_rows,fid)
            ann=input.annotations[input.annotation_rows[fid],:]
            table.annotation_x[i]=_pf_number(ann.x);table.annotation_y[i]=_pf_number(ann.y)
            table.uncertainty_px[i]=_pf_number(ann.uncertainty_px)
            table.visibility[i]=_pf_text(ann.visibility);table.medium[i]=_pf_text(ann.medium)
        end
        for (cohort,lookup) in ((input.baseline,input.baseline_rows),(input.cohort,input.cohort_rows))
            haskey(lookup,id) || continue
            row=cohort[lookup[id],:]
            "split" in names(cohort) && (table.split[i]=_pf_text(row.split))
            "eligibility_status" in names(cohort) && (table.eligibility_status[i]=_pf_text(row.eligibility_status))
        end
    end
    table.original_visibility=copy(table.visibility);table.original_uncertainty_px=copy(table.uncertainty_px)
    for prefix in ("final","evaluation","baseline","predicted")
        table[!,prefix*"_error_px"]=[table[i,prefix*"_valid"] && isfinite(table.annotation_x[i]) && isfinite(table.annotation_y[i]) ?
            hypot(table[i,prefix*"_x"]-table.annotation_x[i],table[i,prefix*"_y"]-table.annotation_y[i]) : NaN for i in 1:n]
    end
    return table
end

"""
    project_calibration_events(original_dir, model_dir; resume=true, python=nothing)

Project every original candidate using its frozen same-session model and saved
raw minmax2/GCC delays, or authenticated saved SRP review delays. Missing SRP
review rows use the identical fixed estimator on authenticated exact prepared
snippets. No refitting, original-media reads/refiltering, annotation selection,
or changes to frozen products. Missing/invalid cache rows remain explicit.
Original zero-origin event/frame/time fields are retained. Pixels are original
resolution x-right/y-down; invalid projections are NaN. Visual truth and frozen
cohort split/status are joined only after all rays are produced. Errors on these
joined estimates are descriptive, never new heldout accuracy.

Writes only model_dir/inference_events.csv and inference_events.json. Matching
resume verifies CSV integrity; stale/incomplete sidecars require resume=false.
Source/model tampering always fails, including with resume=false.
"""
function project_calibration_events(original_dir,model_dir;resume=true,python=nothing)
    original_dir=realpath(original_dir);model_dir=realpath(model_dir)
    try
        input=_ep_inputs(original_dir,model_dir;python=python)
        csv_path=joinpath(model_dir,"inference_events.csv");manifest_path=joinpath(model_dir,"inference_events.json")
        if resume && (isfile(csv_path) || isfile(manifest_path))
            isfile(csv_path) && isfile(manifest_path) || error("Incomplete inference sidecars; explicitly regenerate with resume=false")
            previous=read_json(manifest_path)
            get(previous,"schema",nothing)=="clicker-all-event-projection-v1" && get(previous,"status",nothing)=="complete" || error("Unsupported/incomplete projection manifest")
            previous["inputs"]==input.inputs && previous["source_hashes"]==input.source_hashes &&
                previous["hashes"]==input.hashes && previous["output_csv"]==csv_path &&
                previous["original_dir"]==original_dir && previous["model_dir"]==model_dir || error("Stale projection inputs/sources; explicitly regenerate sidecars with resume=false")
            previous["output_csv_sha256"]==sha256_file(csv_path) || error("Inference CSV was modified")
            return previous
        end
        table=_ep_project(input.events,input.calibration,input.evaluation_calibration,input.saved_srp;cache_dir=input.cache_dir)
        _ep_join_truth!(table,input)
        for (path,digest) in input.inputs
            sha256_file(path)==digest || error("Projection input changed during execution: $path")
        end
        for (name,digest) in input.source_hashes
            sha256_file(joinpath(@__DIR__,name))==digest || error("Projection source changed during execution: $name")
        end
        tmp=csv_path*".tmp."*string(getpid())
        CSV.write(tmp,table)
        open(tmp,"r+") do io;ccall(:fsync,Cint,(Cint,),Base.fd(io));end
        mv(tmp,csv_path;force=true)
        counts=Dict("n_original_events"=>nrow(table),"n_projected"=>count(table.final_valid),
            "n_in_frame"=>count(table.final_in_frame),"n_mapped"=>count(==("mapped"),table.sync_status),
            "projection_status"=>Dict(status=>count(==(status),table.final_status) for status in unique(table.final_status)))
        manifest=Dict("schema"=>"clicker-all-event-projection-v1","status"=>"complete",
            "created_utc"=>string(Dates.now(Dates.UTC))*"Z","original_dir"=>original_dir,"model_dir"=>model_dir,
            "device"=>input.calibration["device"],"method"=>input.calibration["method"],"hashes"=>input.hashes,
            "inputs"=>input.inputs,"source_hashes"=>input.source_hashes,"policy_bridge"=>input.bridge,
            "output_csv"=>csv_path,"output_csv_sha256"=>sha256_file(csv_path),"counts"=>counts,
            "semantics"=>Dict("refit"=>false,"annotation_selects_or_changes_ray"=>false,
                "row_scope"=>"every original candidate in original order, including purged/nonannotated/unmapped events",
                "annotation_authority"=>"unchanged original annotations joined only after acoustic projection",
                "evaluation"=>"same-session descriptive replay; frozen cohort evaluation remains unchanged",
                "srp_missing"=>"missing saved review rows use the frozen fixed SRP estimator on authenticated prepared snippets; missing/invalid snippets remain explicit invalid rows",
                "medium"=>"fixed underwater model is replayed for all candidates; air/mixed/unknown predictions are not physically qualified",
                "units"=>"raw relative delays in recorder samples; seconds and frame indices zero-origin; original-resolution pixels x-right/y-down"))
        write_json(manifest_path,manifest)
        return manifest
    catch exception
        log_failure("project_calibration_events $model_dir",exception,catch_backtrace())
        rethrow()
    end
end
