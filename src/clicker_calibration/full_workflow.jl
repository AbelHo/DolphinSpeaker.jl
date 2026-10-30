# One-call composition of existing scientific and delivery stages. Definitions only.

_full_device_option(value,name)=value isa AbstractDict ? get(value,name,get(value,Symbol(name),nothing)) : value
_full_inside(path,root)=path==root || startswith(path,root*string(Base.Filesystem.path_separator))
_full_path(path)=ispath(path) ? realpath(path) : abspath(path)

function _full_video_matches(manifest,original_dir,output_dir;model_dir=nothing,annotation_only=false)
    get(manifest,"schema_version",nothing)=="clicker-video-export/1.0" && get(manifest,"status",nothing)=="complete" || return false
    get(get(manifest,"settings",Dict()),"annotation_only",nothing)===annotation_only || return false
    get(manifest,"mode",nothing)==(annotation_only ? "annotation_only" : "calibrated_comparison") || return false
    hashes=get(manifest,"inputs_sha256",Dict())
    annotation=joinpath(_full_path(original_dir),"annotations.csv")
    isfile(annotation) && get(hashes,annotation,nothing)==sha256_file(annotation) || return false
    if !annotation_only
        model_dir===nothing && return false
        calibration=joinpath(_full_path(model_dir),"calibration.json")
        isfile(calibration) && get(hashes,calibration,nothing)==sha256_file(calibration) || return false
    end
    destination=_full_path(output_dir)
    for key in ("video_path","reference_audio_path")
        path=get(manifest,key,nothing)
        path isa AbstractString && isfile(path) && _full_inside(_full_path(path),destination) || return false
    end
    return true
end

function _full_restore_video_manifest(original_dir,output_dir;model_dir=nothing,annotation_only=false)
    current=joinpath(output_dir,"export_manifest.json")
    if isfile(current) && _full_video_matches(read_json(current),original_dir,output_dir;model_dir=model_dir,annotation_only=annotation_only)
        return false
    end
    saved=joinpath(output_dir,annotation_only ? "review_export_manifest.json" : "comparison_export_manifest.json")
    isfile(saved) || return false
    bytes=read(saved);manifest=JSON.parse(String(copy(bytes)))
    _full_video_matches(manifest,original_dir,output_dir;model_dir=model_dir,annotation_only=annotation_only) || return false
    # Restore only the mode pointer, byte-exact. The actual exporter must still
    # validate ALL sources/settings/output hashes and decide whether to reuse.
    temporary=current*".restore."*string(getpid())
    open(temporary,"w") do io
        write(io,bytes);flush(io);ccall(:fsync,Cint,(Cint,),Base.fd(io))
    end
    mv(temporary,current;force=true)
    return true
end

function _full_export_video(original_dir,output_dir;model_dir=nothing,annotation_only=false,
        python=nothing,resume=true,ffmpeg="/usr/bin/ffmpeg",ffprobe="/usr/bin/ffprobe",
        threads=0,encoder="auto",decoder="auto",crf=20,preset="fast",runner=run_logged)
    exe=_workflow_python(python);helper=joinpath(@__DIR__,"video_export.py")
    command=annotation_only ? `$exe $helper $original_dir $output_dir --annotation-only` :
        `$exe $helper $original_dir $model_dir $output_dir`
    command=`$command --ffmpeg $ffmpeg --ffprobe $ffprobe --threads $threads --encoder $encoder --decoder $decoder --crf $crf --preset $preset`
    resume || (command=`$command --no-resume`)
    resume && _full_restore_video_manifest(original_dir,output_dir;model_dir=model_dir,annotation_only=annotation_only)
    # The exporter prints a resource-selection record before its pretty JSON
    # result. Its atomic manifest, not the mixed stdout stream, is authoritative.
    runner(command;operation="$(annotation_only ? "Terra review" : "calibrated comparison") video",capture=true)
    manifest=read_json(joinpath(output_dir,"export_manifest.json"))
    _full_video_matches(manifest,original_dir,output_dir;model_dir=model_dir,annotation_only=annotation_only) ||
        error("Video exporter did not produce a complete matching mode/source/model/destination manifest")
    return manifest
end

function _full_package(original_dir,destination;model_dir=nothing,video_dir,final=false,python=nothing)
    exe=_workflow_python(python);helper=joinpath(@__DIR__,"package_results.py")
    command=`$exe $helper $original_dir $destination --video-dir $video_dir`
    model_dir===nothing || (command=`$command --model-dir $model_dir`)
    final && (command=`$command --final`)
    result=JSON.parse(String(run_logged(command;operation="$(final ? "final" : "initial review") portable delivery",capture=true)))
    manifest_path=joinpath(destination,"delivery_manifest.json")
    manifest=read_json(manifest_path)
    video_key=final ? "comparison_video" : "review_video"
    video=get(manifest,video_key,nothing)
    video===nothing && error("Portable package lacks $video_key")
    result["delivery_manifest"]=manifest_path
    result["delivery_manifest_sha256"]=sha256_file(manifest_path)
    result["delivered_video_sha256"]=manifest["files_sha256"][video]
    return result
end

function _full_notify_delivery(delivery;phase,notify=:telegram,telegram_sender=nothing,resume=true)
    notify==:none && return Dict("status"=>"disabled")
    notify==:telegram || throw(ArgumentError("notify must be :telegram or :none"))
    phase in (:annotation,:final) || throw(ArgumentError("Unknown delivery notification phase"))
    destination=delivery["destination"]
    video=delivery[phase==:annotation ? "review_video" : "comparison_video"]
    identity=Dict("phase"=>string(phase),"annotations_sha256"=>delivery["annotations_sha256"],
        "calibration_sha256"=>phase==:annotation ? nothing : delivery["calibration_sha256"],
        "destination"=>destination,"video"=>video,"video_sha256"=>delivery["delivered_video_sha256"])
    receipt=joinpath(destination,phase==:annotation ? "annotation_ready_notification.json" : "final_ready_notification.json")
    if resume && isfile(receipt)
        previous=read_json(receipt)
        get(previous,"status",nothing)=="delivered" && get(previous,"identity",nothing)==identity && return previous
    end
    sender=telegram_sender===nothing ? get(ENV,"CLICKER_TELEGRAM_SENDER",joinpath(homedir(),".agents","skills","telegram-webhook","send.sh")) : string(telegram_sender)
    isfile(sender) || error("Telegram sender is unavailable: $sender")
    description=phase==:annotation ? "Terra gallery and synchronized hydrophone-audio video are ready. Calibration processing continues while you review." :
        "Calibration comparison gallery, plots, synchronized video and portable results are ready. See the model reports for coverage and fit limitations."
    message="Clicker $(basename(destination)): $description\nReview: $(joinpath(destination,"index.html"))\nVideo: $(joinpath(destination,video))\nAnnotations: $(joinpath(destination,"annotations.csv"))"
    run_logged(`$sender --text $message`;operation="$(phase) delivery-ready notification",capture=true)
    result=Dict("status"=>"delivered","identity"=>identity,"sent_utc"=>string(Dates.now(Dates.UTC))*"Z")
    write_json(receipt,result)
    return result
end

function _full_operations(overrides)
    # Lazy wrappers keep inclusion side-effect free and permit minimal meaningful
    # orchestration fixtures without importing or impersonating numerical code.
    defaults=(
        prepare=(args...;kwargs...)->prepare_device(args...;kwargs...),
        sync=(args...;kwargs...)->review_sync_eligibility(args...;kwargs...),
        annotate=(args...;kwargs...)->request_annotations(args...;kwargs...),
        annotation_current=(args...)->_workflow_annotation_current(args...),
        validate_annotation=(dir;python=nothing)->_workflow_annotation_tool(dir,"validate";python=python),
        gallery=(args...;kwargs...)->render_review(args...;kwargs...),
        video=_full_export_video,package=_full_package,notify=_full_notify_delivery,
        fit_input=(args...;kwargs...)->prepare_fit_annotation_input(args...;kwargs...),
        policy=(args...;kwargs...)->predeclare_observation_policy(args...;kwargs...),
        fit=(args...;kwargs...)->fit_device(args...;kwargs...),
        review=(args...;kwargs...)->run_estimator_review(args...;kwargs...),
        refine=(args...;kwargs...)->refine_calibration(args...;kwargs...),
        project=(args...;kwargs...)->project_calibration_events(args...;kwargs...),
        diagnostics=(args...;kwargs...)->render_diagnostics(args...;kwargs...),
        clock=time,pause=sleep)
    overrides isa NamedTuple || throw(ArgumentError("operations must be a NamedTuple of stage callbacks"))
    all(in(keys(defaults)),keys(overrides)) || throw(ArgumentError("Unknown operation override"))
    return merge(defaults,overrides)
end

function _full_existing_science(fit_dir,model_dir;delivery_only=false)
    if delivery_only
        all(isfile(joinpath(model_dir,name)) for name in ("calibration.json","predictions.csv")) &&
            (isfile(joinpath(model_dir,"report.json")) || isfile(joinpath(model_dir,"fit_report.json"))) ||
            error("Frozen delivery requires existing model calibration, predictions and report: $model_dir")
        return true
    end
    present=false
    for (directory,name,source_keys) in (
            (fit_dir,"fit_report.json",(("physical_fit_source_sha256","physical_fit.jl"),("support_source_sha256","support.jl"))),
            (model_dir,"report.json",(("source/physical_fit.jl","physical_fit.jl"),("source/prepare.jl","prepare.jl"),
                ("source/srp_review.jl","srp_review.jl"),("source/method_refinement.jl","method_refinement.jl"),("source/support.jl","support.jl"))))
        path=joinpath(directory,name)
        isfile(path) || continue
        present=true;report=read_json(path);hashes=get(report,"hashes",Dict())
        for (key,source) in source_keys
            get(hashes,key,nothing)==sha256_file(joinpath(@__DIR__,source)) ||
                error("Existing scientific results use different source ($source). Use delivery_only=true for frozen delivery, or a new run directory; existing science will not be refitted.")
        end
    end
    if !present && any(isfile(joinpath(fit_dir,name)) for name in ("calibration.json","heldout_evaluation.json"))
        error("Partial frozen scientific results need explicit recovery in their original source version; refusing automatic refit")
    end
    return present
end

"""
    run_full_clicker_calibration(input_root; output_dir, results_dir,
        devices=(:A,:B), annotator=:herdr, notify=:telegram, resume=true,
        observation_policy=nothing, delivery_only=false, ...)

Prepare/synchronize/request independent annotations; deliver each ready device's
Terra gallery and synchronized hydrophone-audio video BEFORE fitting. Then run
first fit, estimator review, refinement, all-candidate acoustic projection, final
gallery/plots/comparison video and `package_results.py --final`. Notifications
announce actual packaged deliveries and never require manual approval to fit.
Heavy stages execute sequentially. Video defaults are threads=0, encoder=auto,
decoder=auto; this function introduces no CPU/GPU/memory caps or nested job queue.
Run it under the shared job launcher when process-level serialization is needed.

`annotation_fit_policy`/`fit_input_root` support the existing explicit A adapter.
`observation_policy` opts into the native policy without changing original labels.
Either policy may be scalar or keyed by string/symbol device; do not combine both
for one device. The native default remains strict. `prepare_options`,
`fit_options`, and `video_options` pass through to their respective stages.

`delivery_only=true` reuses existing originals/models via `original_dirs` and
`model_dirs` (device-keyed dictionaries or scalar paths for a single device).
It never prepares, launches annotation, declares policy, fits, reviews estimators,
or refines. Existing old-source scientific runs are rejected by normal mode;
use frozen delivery mode or a new run directory, never silent refitting.

Returns an explicit `awaiting_annotation` state when `wait_annotations=false`
or the annotation deadline expires, preserving completed devices. Original
annotations remain untouched. Stage failures are recorded and rethrown.
`operations` overrides named stage callbacks for orchestration tests; these
callbacks must retain the existing stage contracts, not weaken science checks.
"""
function run_full_clicker_calibration(input_root::AbstractString;output_dir::AbstractString,results_dir::AbstractString,
        devices=(:A,:B),annotator=:herdr,notify=:telegram,resume=true,observation_policy=nothing,
        annotation_fit_policy=nothing,fit_input_root=nothing,delivery_only=false,original_dirs=nothing,model_dirs=nothing,
        wait_annotations=true,annotation_timeout_s=21600,python=nothing,plot_python=nothing,telegram_sender=nothing,
        diagnostic_plots=true,prepare_options=NamedTuple(),fit_options=NamedTuple(),video_options=NamedTuple(),
        error_log=get(ENV,"CLICKER_ERROR_LOG",joinpath(pwd(),"docs","error.log")),operations=NamedTuple())
    names=string.(collect(devices))
    !isempty(names) && length(unique(names))==length(names) && all(in(("A","B")),names) || throw(ArgumentError("devices must be unique A/B identifiers"))
    annotator in (:herdr,:external) || throw(ArgumentError("annotator must be :herdr or :external"))
    notify in (:telegram,:none) || throw(ArgumentError("notify must be :telegram or :none"))
    annotation_timeout_s>0 || throw(ArgumentError("annotation_timeout_s must be positive"))
    model_dirs===nothing || delivery_only || throw(ArgumentError("Explicit existing model_dirs require delivery_only=true"))
    for value in (original_dirs,model_dirs)
        value===nothing || value isa AbstractDict || length(names)==1 || throw(ArgumentError("Use per-device directory dictionaries for multiple devices"))
    end
    all(in((:ffmpeg,:ffprobe,:threads,:encoder,:decoder,:crf,:preset)),keys(video_options)) || throw(ArgumentError("Unsupported video option"))
    input_root=_full_path(input_root);output_dir=_full_path(output_dir);results_dir=_full_path(results_dir)
    recording_dirs=[_full_path(joinpath(input_root,name)) for name in ("A","B")]
    for path in (output_dir,results_dir)
        any(recording->_full_inside(path,recording),recording_dirs) && error("Generated work/delivery must be outside source recording directories")
    end
    ops=_full_operations(operations);paths=Dict{String,Any}()
    for name in names
        explicit=_full_device_option(original_dirs,name)
        original=explicit===nothing ? joinpath(output_dir,name) : _full_path(explicit)
        adapter=_full_device_option(annotation_fit_policy,name);observed=_full_device_option(observation_policy,name)
        adapter===nothing || observed===nothing || throw(ArgumentError("Choose one observation-admission policy per device"))
        fit=adapter===nothing ? original : joinpath(fit_input_root===nothing ? joinpath(output_dir,"fit_inputs") : _full_path(fit_input_root),name)
        explicit_model=_full_device_option(model_dirs,name)
        model=explicit_model===nothing ? joinpath(fit,"refinement_v1") : _full_path(explicit_model)
        delivery_only && (fit=isfile(joinpath(model,"report.json")) ? dirname(model) : model)
        destination=joinpath(results_dir,name)
        for path in (original,fit,model,destination)
            any(recording->_full_inside(path,recording),recording_dirs) && error("Workflow paths must not write into source recording directories")
        end
        (_full_inside(destination,original) || _full_inside(original,destination)) && error("Portable delivery must be separate from original scientific inputs")
        existing=_full_existing_science(fit,model;delivery_only=delivery_only)
        # The exporter publishes distinct authenticated generations. Package
        # directly from the destination's videos directory: same-file copies
        # validate hashes and avoid a second full-size media copy.
        paths[name]=(;original,fit,model,destination,adapter,observed,existing,
            review_video=joinpath(destination,"videos"),final_video=joinpath(destination,"videos"))
    end
    mkpath(output_dir);mkpath(results_dir);CALIBRATION_ERROR_LOG[]=abspath(error_log)
    state_path=joinpath(output_dir,"full_workflow_manifest.json")
    state=Dict{String,Any}("schema"=>"clicker-full-workflow-v1","status"=>"running",
        "input_root"=>input_root,"output_dir"=>output_dir,"results_dir"=>results_dir,"delivery_only"=>delivery_only,
        "source_sha256"=>sha256_file(@__FILE__),"started_utc"=>string(Dates.now(Dates.UTC))*"Z",
        "devices"=>Dict(name=>Dict{String,Any}("original_dir"=>paths[name].original,"fit_dir"=>paths[name].fit,
            "model_dir"=>paths[name].model,"destination"=>paths[name].destination,"status"=>"pending") for name in names))
    write_json(state_path,state)
    current_device="";current_stage="initialization"
    function stage!(name,label,action)
        current_device=name;current_stage=label
        result=state["devices"][name];result["active_stage"]=label;write_json(state_path,state)
        value=action();result[label]=value;write_json(state_path,state)
        return value
    end
    try
        pending=String[]
        for name in names
            p=paths[name];result=state["devices"][name]
            if !delivery_only && !p.existing
                prepared=stage!(name,"prepare",()->ops.prepare(joinpath(input_root,name),p.original;device=name,resume=resume,prepare_options...))
                if get(prepared,"status",nothing)!="complete"
                    result["status"]="needs_preparation_review";continue
                end
                stage!(name,"sync",()->ops.sync(p.original;resume=resume))
            end
            if !delivery_only
                annotation=stage!(name,"annotation_request",()->ops.annotate(p.original;annotator=annotator,python=python,resume=resume))
                if get(annotation,"status",nothing)=="no_frames"
                    result["status"]="no_frames_to_annotate";continue
                end
            end
            push!(pending,name)
        end
        deadline=ops.clock()+annotation_timeout_s
        while !isempty(pending)
            ready=String[]
            # Publish every currently ready visual review before fitting any of
            # this batch. Independent annotation agents may continue meanwhile.
            for name in pending
                p=paths[name];result=state["devices"][name]
                if !ops.annotation_current(p.original)
                    result["status"]="awaiting_annotation";continue
                end
                stage!(name,"annotation_validation",()->ops.validate_annotation(p.original;python=python))
                stage!(name,"annotation_gallery",()->ops.gallery(p.original;python=python))
                stage!(name,"annotation_video",()->ops.video(p.original,p.review_video;
                    annotation_only=true,python=python,resume=resume,video_options...))
                delivery=stage!(name,"initial_delivery",()->ops.package(p.original,p.destination;
                    video_dir=p.review_video,final=false,python=python))
                stage!(name,"initial_notification",()->ops.notify(delivery;phase=:annotation,notify=notify,telegram_sender=telegram_sender,resume=resume))
                result["status"]="review_delivered";push!(ready,name)
            end
            for name in ready
                p=paths[name];result=state["devices"][name]
                if !delivery_only
                    p.adapter===nothing || stage!(name,"annotation_fit_policy",()->ops.fit_input(p.original,p.fit;policy=p.adapter,python=python))
                    p.observed===nothing || stage!(name,"observation_policy",()->ops.policy(p.fit;policy=p.observed,
                        reason="Exploratory calibration from sealed Terra estimates with original uncertainty; no human-confirmed positions are asserted."))
                    stage!(name,"fit",()->ops.fit(p.fit;resume=resume,fit_options...))
                    stage!(name,"estimator_review",()->ops.review(p.fit;resume=resume))
                    stage!(name,"refinement",()->ops.refine(p.fit;resume=resume))
                end
                stage!(name,"all_event_projection",()->ops.project(p.original,p.model;resume=resume,python=python))
                stage!(name,"comparison_gallery",()->ops.gallery(p.original;fitted=true,model_dir=p.model,python=python))
                diagnostic_plots && stage!(name,"diagnostics",()->ops.diagnostics(p.fit;model_dir=p.model,python=plot_python))
                stage!(name,"comparison_video",()->ops.video(p.original,p.final_video;
                    model_dir=p.model,annotation_only=false,python=python,resume=resume,video_options...))
                delivery=stage!(name,"final_delivery",()->ops.package(p.original,p.destination;
                    model_dir=p.model,video_dir=p.final_video,final=true,python=python))
                get(delivery,"status",nothing)=="complete" || error("Final package did not report completion")
                stage!(name,"final_notification",()->ops.notify(delivery;phase=:final,notify=notify,telegram_sender=telegram_sender,resume=resume))
                result["status"]="complete";result["active_stage"]=nothing
                filter!(!=(name),pending);write_json(state_path,state)
            end
            isempty(pending) && break
            if !wait_annotations || delivery_only || ops.clock()>=deadline
                state["status"]="awaiting_annotation";state["pending_devices"]=copy(pending)
                write_json(state_path,state);return state
            end
            write_json(state_path,state);ops.pause(5)
        end
        state["status"]=all(state["devices"][name]["status"]=="complete" for name in names) ? "complete" : "needs_review"
        state["completed_utc"]=string(Dates.now(Dates.UTC))*"Z";write_json(state_path,state)
        return state
    catch exception
        state["status"]="failed";state["failed_device"]=current_device;state["failed_stage"]=current_stage
        state["error"]=sprint(showerror,exception);write_json(state_path,state)
        log_failure("run_full_clicker_calibration device=$current_device stage=$current_stage",exception,catch_backtrace())
        rethrow()
    end
end
