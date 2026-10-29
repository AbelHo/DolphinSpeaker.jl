# Integration, annotation handoff, and reporting. No import-time I/O.

function load_calibration(path::AbstractString; device=nothing)
    calibration = read_json(path)
    calibration["schema_version"]=="clicker-physical-calibration/1.0" || error("Unsupported calibration schema")
    isnothing(device) || string(calibration["device"]) == string(device) || error("Calibration belongs to a different device")
    calibration["physical_channels"] == [3,4,5] || error("Unexpected physical channel mapping")
    calibration["reference_channel"] == 3 || error("Unexpected reference channel")
    1400 <= calibration["sound_speed_m_s"] <= 1600 || error("Sound speed outside calibrated bounds")
    geometry=_pf_matrix(calibration["geometry_m"])
    size(geometry)==(3,3) && all(isfinite,geometry) && maximum(abs.(geometry.-clicker_geometry()))<1e-12 || error("Calibration geometry is not the confirmed D3 array")
    rotation=_pf_matrix(calibration["rotation_array_to_camera"])
    size(rotation)==(3,3) && all(isfinite,rotation) &&
        norm(rotation' * rotation-I)<1e-9 && abs(det(rotation)-1)<1e-9 || error("Calibration rotation must be proper SO(3)")
    offsets=Float64.(calibration["relative_delay_offsets_s"])
    length(offsets)==3 && all(isfinite,offsets) && abs(offsets[3])<=1e-12 && maximum(abs.(offsets))<=1e-4+eps(Float64) || error("Invalid relative channel delay offsets")
    for key in ("width","height")
        value=calibration[key]
        value isa Integer && value>0 || error("Image dimensions must be positive integers")
    end
    camera=calibration["camera"]
    all(key->camera[key] isa Real && isfinite(camera[key]),("fx","fy","cx","cy","k1","k2")) || error("Camera parameters must be finite")
    camera["fx"]>0 && camera["fy"]>0 || error("Focal lengths must be positive")
    return calibration
end

function _workflow_python(python)
    isnothing(python) ? get(ENV, "CLICKER_PYTHON", "python3") : string(python)
end

function _workflow_annotation_tool(output_dir, action; python=nothing)
    exe = _workflow_python(python)
    helper = joinpath(@__DIR__, "annotation_tools.py")
    bytes = run_logged(`$exe $helper $action $output_dir`;
        operation="annotation $(action) $(basename(output_dir))", capture=true)
    return JSON.parse(String(bytes))
end

function _workflow_annotation_current(output_dir)
    path=joinpath(output_dir,"annotations.complete.json")
    isfile(path) || return false
    a=joinpath(output_dir,"annotations.csv"); f=joinpath(output_dir,"frames.csv")
    isfile(a) && isfile(f) || return false
    m=read_json(path)
    return get(m,"annotations_sha256",nothing)==sha256_file(a) && get(m,"frames_sha256",nothing)==sha256_file(f)
end

function _workflow_annotation_prompt(output_dir, python)
    frames=abspath(joinpath(output_dir,"frames.csv"))
    helper=joinpath(@__DIR__,"annotation_tools.py")
    return """
User-authorized visual clicker annotation task. Parent calibration lead owns all
scientific processing. Use gpt-6.1-terra high. Do not spawn other agents.

Read $frames. Inspect EVERY unique image individually using image tools. You may
use labeled contact sheets to navigate, but inspect native image/crops wherever
needed. Do not infer unseen points, interpolate a track, copy one coordinate
across frames, or use acoustic predictions. No automatic object detector is a
substitute for your requested image inspection. Do not open events.csv,
predictions.csv, fit reports or calibration outputs when annotating.

Target: the CENTER OF THE CLEARLY IDENTIFIABLE CLICKER HOUSING, consistently in
every image. In representative footage it is a small blue handheld clicker with
a red button and a cord. Button visibility is not required when the housing
outline and center are clear; use uncertainty covering the plausible sound
emitter position within that housing and note button occlusion. Never substitute
the hand, swimmer, hydrophone, cord alone, or an acoustic prediction. This visual
cue is not an assumption of presence: a mostly hidden/indistinct housing or
multiple plausible targets must be uncertain/occluded. Do not alternate between
button-center and housing-center labels. Initial frames may be on the pool deck;
later frames are underwater.

Write only within $(abspath(output_dir)): annotations.pending.csv, optional
visual-inspection crops/contact sheets under annotation_work/, then final
annotations.csv and annotations.complete.json. Keep every frames.csv identity,
one annotation row per frame, including non-visible frames. Schema:
frame_id,x,y,visibility,medium,uncertainty_px,note,annotator
visibility: visible/uncertain/occluded/absent.
medium: underwater/air/mixed/unknown; mixed includes an uncertain water boundary.
x,y are ZERO-ORIGIN ORIGINAL-RESOLUTION pixels, x right/y down, never coordinates
of a resized display/crop. Visible rows need finite x,y and positive uncertainty
in original pixels. Others may leave coordinates empty. Use a conservative
uncertainty when the clicker's edges are blurred. Record brief ambiguity notes.
annotator must identify gpt-6.1-terra. Retain inspection evidence/notes for all
frames; never mark all absent without viewing them. Save progress atomically.

When ALL rows are complete, atomically replace annotations.csv, then execute
the Python helper with the argument vector [$python, $helper, "seal",
$(abspath(output_dir))]. It validates coverage, identity, coordinates, dimensions
and writes the completion manifest. Do not seal partial work. Send the parent
a DONE: message with counts, final paths and any visibility limits. The parent
will generate overlays and notify the user on Telegram; you need not send that
notification. No input media edits, source edits, commits or other writes.
"""
end

"""Create the image-only handoff and optionally launch the specified Terra agent."""
function _workflow_live_annotation_agent(name)
    bytes=run_logged(`hagent list`;operation="inspect existing annotation worker",capture=true)
    listing=JSON.parse(String(bytes))
    agents=listing["result"]["agents"]
    matches=filter(a->get(a,"name",nothing)==name,agents)
    return isempty(matches) ? nothing : only(matches)
end

function _workflow_spawn_annotation(name,prompt)
    run_logged(`hagent spawn codex $name --model gpt-6.1-terra --effort high -- $prompt`;
        operation="launch Terra annotation $name",capture=true)
end

_workflow_annotation_name(device) = lowercase("clickann-$(device)-$(Dates.format(now(UTC),"mmddHHMMSS"))-$(randstring(4))")

function request_annotations(output_dir::AbstractString; annotator=:herdr, python=nothing, resume=true,
    agent_lookup=_workflow_live_annotation_agent,agent_spawn=_workflow_spawn_annotation)
    frames=joinpath(output_dir,"frames.csv")
    isfile(frames) || error("Run preparation before requesting annotations")
    table=CSV.read(frames,DataFrame)
    nrow(table)>0 || return Dict("status"=>"no_frames","frame_count"=>0)
    _workflow_annotation_current(output_dir) && return Dict("status"=>"ready")
    request_path=joinpath(output_dir,"annotation_request.json")
    frames_sha=sha256_file(frames)
    previous=nothing
    if resume && annotator==:herdr && isfile(request_path)
        previous=read_json(request_path)
        if get(previous,"frames_sha256",nothing)==frames_sha && get(previous,"status",nothing)=="launched"
            agent=agent_lookup(previous["agent_name"])
            if agent !== nothing && get(agent,"agent_status","") in ("working","running")
                return previous
            end
            # A saved launch receipt is not liveness. Retain it as recovery
            # evidence, then launch a replacement for unfinished annotation.
            write_json(joinpath(output_dir,"annotation_request.previous.json"),previous)
        end
    end
    prompt=_workflow_annotation_prompt(output_dir,_workflow_python(python))
    prompt_path=joinpath(output_dir,"annotation_prompt.txt")
    write(prompt_path,prompt)
    request=Dict{String,Any}("status"=>"awaiting_external_annotation","frames_sha256"=>frames_sha,
        "frame_count"=>nrow(table),"prompt_path"=>prompt_path,"model"=>"gpt-6.1-terra","effort"=>"high")
    if annotator == :herdr
        device=replace(basename(output_dir),r"[^A-Za-z0-9]"=>"")
        name=_workflow_annotation_name(device)
        # hagent is initialized by the caller/lead once in its Herdr pane.
        agent_spawn(name,prompt)
        request["status"]="launched";request["agent_name"]=name
    elseif annotator != :external
        throw(ArgumentError("annotator must be :herdr or :external"))
    end
    write_json(request_path,request)
    return request
end

function render_review(output_dir::AbstractString; fitted=false,python=nothing)
    _workflow_annotation_tool(output_dir,fitted ? "fit-gallery" : "gallery";python=python)
end

function _workflow_notify_annotation(output_dir,review;notify=:telegram,telegram_sender=nothing)
    notify == :none && return Dict("status"=>"disabled")
    notify == :telegram || throw(ArgumentError("notify must be :telegram or :none"))
    sender=isnothing(telegram_sender) ? get(ENV,"CLICKER_TELEGRAM_SENDER",
        joinpath(homedir(),".agents","skills","telegram-webhook","send.sh")) : telegram_sender
    isfile(sender) || error("Telegram sender is unavailable: $sender")
    receipt=joinpath(output_dir,"annotation_notification.json")
    if isfile(receipt)
        old=read_json(receipt)
        get(old,"annotations_sha256",nothing)==review["annotations_sha256"] && get(old,"status",nothing)=="delivered" && return old
    end
    message="Clicker calibration $(basename(output_dir)): annotations and verification overlays are ready.\n" *
        "Review gallery: $(review["gallery"])\n" *
        "Annotations: $(joinpath(abspath(output_dir),"annotations.csv"))\n" *
        "Counts: $(JSON.json(review["counts"]))\n" *
        "Fitting continues while you check. Tell me the frame ID if a correction is needed."
    run_logged(`$sender --text $message`;operation="Telegram annotation-ready notification",capture=true)
    result=Dict("status"=>"delivered","annotations_sha256"=>review["annotations_sha256"],
        "gallery"=>review["gallery"],"sent_utc"=>string(now(UTC))*"Z")
    write_json(receipt,result)
    return result
end

function run_estimator_review(output_dir;resume=true)
    source=joinpath(@__DIR__,"srp_review.jl")
    isfile(source) || error("Estimator review implementation is missing")
    if !isdefined(@__MODULE__,:review_estimators)
        Base.include(@__MODULE__,source)
    end
    Base.invokelatest() do
        getfield(@__MODULE__,:review_estimators)(output_dir;resume=resume)
    end
end

"""
    run_clicker_calibration(input_root; devices=(:A,:B), output_dir,
        annotator=:herdr, notify=:telegram, resume=true, stage=:all,
        wait_annotations=true, annotation_timeout_s=21600, ...)

Prepare, request independent Terra annotations, notify their review location,
fit each device separately, and compare existing estimators after the first fit.
Stages :prepare, :annotate, :fit, :review, and :all are independently resumable.
External CSV annotation is supported with annotator=:external. The only normal
manual step is visual annotation; :all can wait for atomic completion manifests.
Generated artifacts go under output_dir/A and output_dir/B. Source input is
read-only. No configuration file is rewritten by this function; fitted JSONs
are loaded through the separate explicit device preset setters.
"""
function run_clicker_calibration(input_root::AbstractString; devices=(:A,:B),output_dir::AbstractString,
    annotator=:herdr,notify=:telegram,resume=true,stage=:all,wait_annotations=true,
    annotation_timeout_s=21600,python=nothing,telegram_sender=nothing,
    error_log=get(ENV,"CLICKER_ERROR_LOG",joinpath(pwd(),"docs","error.log")),
    prepare_options=NamedTuple(),fit_options=NamedTuple())
    stage in (:all,:prepare,:annotate,:fit,:review) || throw(ArgumentError("Unknown calibration stage"))
    input_root=abspath(input_root);output_dir=abspath(output_dir)
    (output_dir == input_root || startswith(output_dir,input_root*"/")) && error("Output must be outside read-only source recordings")
    mkpath(output_dir)
    CALIBRATION_ERROR_LOG[]=abspath(error_log)
    sources=Dict(basename(p)=>sha256_file(p) for p in readdir(@__DIR__;join=true) if endswith(p,".jl") || endswith(p,".py"))
    state=Dict{String,Any}("schema_version"=>"clicker-run-v1","status"=>"running","stage"=>string(stage),
        "input_root"=>input_root,"output_dir"=>output_dir,"devices"=>Dict{String,Any}(),
        "source_hashes"=>sources,"julia_version"=>string(VERSION),"started_utc"=>string(now(UTC))*"Z")
    state_path=joinpath(output_dir,"run_manifest.json")
    write_json(state_path,state)
    try
        for device in devices
            name=string(device);name in ("A","B") || error("Device must be A or B")
            dest=joinpath(output_dir,name);mkpath(dest)
            result=Dict{String,Any}();state["devices"][name]=result
            if !(stage in (:all,:prepare)) && !isfile(joinpath(dest,"frames.csv"))
                error("Missing prepared frames for $name; run stage=:prepare first")
            end
            if stage in (:all,:prepare)
                result["prepare"]=prepare_device(joinpath(input_root,name),dest;device=name,resume=resume,prepare_options...)
                write_json(state_path,state)
                if result["prepare"]["status"] != "complete"
                    result["status"]="needs_sync_or_preparation_review"
                    continue
                end
            end
            if stage in (:all,:annotate)
                result["annotations"]=request_annotations(dest;annotator=annotator,python=python,resume=resume)
                if result["annotations"]["status"]=="no_frames"
                    result["status"]="no_frames_to_annotate"
                end
            end
            write_json(state_path,state)
        end
        if stage == :prepare
            state["status"]=all(get(v,"status","")!="needs_sync_or_preparation_review" for v in values(state["devices"])) ? "prepared" : "needs_review"
            write_json(state_path,state);return state
        end
        pending=String[string(d) for d in devices if !(get(state["devices"][string(d)],"status","") in ("needs_sync_or_preparation_review","no_frames_to_annotate"))]
        deadline=time()+annotation_timeout_s
        while !isempty(pending)
            for name in copy(pending)
                dest=joinpath(output_dir,name);result=state["devices"][name]
                if !_workflow_annotation_current(dest)
                    result["status"]="awaiting_annotations"
                    continue
                end
                _workflow_annotation_tool(dest,"validate";python=python)
                result["annotation_review"]=render_review(dest;python=python)
                result["annotation_notification"]=_workflow_notify_annotation(dest,result["annotation_review"];notify=notify,telegram_sender=telegram_sender)
                if stage in (:all,:fit)
                    result["fit"]=fit_device(dest;resume=resume,fit_options...)
                    result["fit_review"]=render_review(dest;fitted=true,python=python)
                end
                if stage in (:all,:review)
                    result["estimator_review"]=run_estimator_review(dest;resume=resume)
                end
                result["status"]=stage==:annotate ? "annotated" : "complete"
                filter!(!=(name),pending)
                write_json(state_path,state)
            end
            isempty(pending) && break
            if !wait_annotations || annotator==:external || time()>deadline
                state["status"]="awaiting_annotations";write_json(state_path,state);return state
            end
            sleep(5)
        end
        state["status"]=all(get(v,"status","") in ("complete","annotated") for v in values(state["devices"])) ? "complete" : "needs_review"
        state["completed_utc"]=string(now(UTC))*"Z"
        write_json(state_path,state)
        return state
    catch exception
        state["status"]="failed";state["error"]=sprint(showerror,exception)
        write_json(state_path,state)
        log_failure("run_clicker_calibration stage=$stage",exception,catch_backtrace())
        rethrow()
    end
end
