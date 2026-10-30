# Scoped postprocessing of frozen preparation. No import-time I/O/computation.
const SYNC_ELIGIBILITY_SCHEMA = "clicker-sync-eligibility-v1"

function _syncelig_identity(path)
    a=stat(path);sha=sha256_file(path);b=stat(path)
    (a.size,a.mtime,a.ctime,a.inode)==(b.size,b.mtime,b.ctime,b.inode) || error("Timing-review input changed while hashing: $path")
    Dict{String,Any}("path"=>abspath(path),"sha256"=>sha,"bytes"=>b.size,
        "mtime"=>b.mtime,"ctime"=>b.ctime,"inode"=>b.inode)
end
_syncelig_digest(bytes)=bytes2hex(SHA.sha256(bytes))
function _syncelig_write_bytes(path,bytes; preserve=false)
    islink(path) && error("Refusing symlink timing-review output: $path")
    if preserve && isfile(path)
        read(path)==bytes || error("Immutable preparation original differs: $path")
        return
    end
    mkpath(dirname(path));tmp=path*".tmp-$(getpid())"
    ispath(tmp) && error("Timing-review temporary already exists: $tmp")
    try
        open(tmp,"w") do io
            write(io,bytes);flush(io);ccall(:fsync,Cint,(Cint,),Base.fd(io))
        end
        # Single integration owner; never replace an existing raw original.
        preserve && ispath(path) && error("Raw original appeared during write: $path")
        mv(tmp,path;force=!preserve)
    finally
        isfile(tmp) && rm(tmp)
    end
end
function _syncelig_write_json(path,value)
    _syncelig_write_bytes(path,Vector{UInt8}(codeunits(JSON.json(json_safe(value),2)*"\n")))
end
function _syncelig_output_entry(manifest,path)
    entries=[q for q in manifest["outputs"] if abspath(q["path"])==abspath(path)]
    length(entries)==1 || error("Expected one preparation output identity: $path")
    only(entries)
end
function _syncelig_check_identity(path,entry)
    ident=_syncelig_identity(path)
    ident["bytes"]==entry["bytes"] && ident["sha256"]==entry["sha256"] || error("Preparation output authentication failed: $path")
    ident
end
function _syncelig_non_event_outputs(manifest,eventpath)
    [deepcopy(q) for q in manifest["outputs"] if abspath(q["path"])!=abspath(eventpath)]
end

# RFC4180 byte spans, including quoted/multiline fields. Used only to replace
# selected status cells; every other byte of the frozen table is retained.
function _syncelig_csv_spans(bytes)
    rows=Vector{Vector{UnitRange{Int}}}();row=UnitRange{Int}[]
    n=length(bytes);first=1;i=1;quoted=false
    while i<=n
        b=bytes[i]
        if b==0x22
            if quoted && i<n && bytes[i+1]==0x22
                i+=2;continue
            elseif quoted
                quoted=false
            elseif i==first
                quoted=true
            else
                error("Malformed quoted CSV field")
            end
        elseif !quoted && (b==0x2c || b==0x0a || b==0x0d)
            push!(row,first:i-1)
            if b!=0x2c
                push!(rows,row);row=UnitRange{Int}[]
                b==0x0d && i<n && bytes[i+1]==0x0a && (i+=1)
            end
            first=i+1
        end
        i+=1
    end
    quoted && error("Unterminated CSV quote")
    if first<=n || !isempty(row)
        push!(row,first:n);push!(rows,row)
    end
    rows
end
function _syncelig_decode_cell(bytes,span)
    raw=bytes[span]
    if length(raw)>=2 && first(raw)==0x22 && last(raw)==0x22
        return replace(String(raw[2:end-1]),"\"\""=>"\"")
    end
    String(raw)
end
function _syncelig_patch_status(bytes,events,additional)
    isempty(additional) && return copy(bytes)
    spans=_syncelig_csv_spans(bytes)
    length(spans)==nrow(events)+1 || error("CSV lexical row count differs from parsed events")
    header=[_syncelig_decode_cell(bytes,r) for r in spans[1]]
    idcol=findfirst(==("event_id"),header);statuscol=findfirst(==("sync_status"),header)
    isnothing(idcol) || isnothing(statuscol) ? error("Missing event/status header") : nothing
    changes=Tuple{UnitRange{Int},Vector{UInt8}}[]
    for (i,row) in enumerate(eachrow(events))
        fields=spans[i+1]
        length(fields)==length(header) || error("CSV column count differs")
        id=String(row.event_id)
        _syncelig_decode_cell(bytes,fields[idcol])==id || error("CSV event identity mismatch")
        if id in additional
            String(row.sync_status)=="mapped" || error("Only original mapped rows may receive an added exclusion")
            r=fields[statuscol];_syncelig_decode_cell(bytes,r)=="mapped" || error("CSV original status differs")
            replacement=bytes[first(r)]==0x22 ? "\"ambiguous\"" : "ambiguous"
            push!(changes,(r,Vector{UInt8}(codeunits(replacement))))
        end
    end
    length(changes)==length(additional) || error("Additional exclusion event missing from raw table")
    io=IOBuffer();cursor=1
    for (span,replacement) in changes
        write(io,bytes[cursor:first(span)-1]);write(io,replacement);cursor=last(span)+1
    end
    write(io,bytes[cursor:end]);take!(io)
end
function _syncelig_check_table_change(raw,patched,additional)
    old=CSV.read(IOBuffer(copy(raw)),DataFrame);new=CSV.read(IOBuffer(copy(patched)),DataFrame)
    names(old)==names(new) && nrow(old)==nrow(new) || error("Timing review changed CSV schema/rows")
    for i in 1:nrow(old),field in propertynames(old)
        expected=field==:sync_status && String(old.event_id[i]) in additional ? "ambiguous" : old[i,field]
        isequal(expected,new[i,field]) || error("Timing review changed forbidden event field: $field row$i")
    end
    new
end

function _syncelig_frame(videos,starts,t,device,video_hashes)
    i=searchsortedlast(starts,t);i==0 && return nothing
    v=videos[i];localtime=t-starts[i]+get(v,"audio_start_s",0.0)
    0<=localtime<v["duration_s"] || return nothing
    pts=v["pts"];k=searchsortedfirst(pts,localtime)
    k=k>length(pts) ? length(pts) : k==1 ? 1 : abs(pts[k-1]-localtime)<=abs(pts[k]-localtime) ? k-1 : k
    Dict{String,Any}("video_file"=>v["path"],"frame_index"=>k-1,"frame_pts_s"=>pts[k],
        "frame_id"=>string(device,"-",video_hashes[v["path"]][1:12],"-f",k-1))
end
_syncelig_frame_key(hit)=isnothing(hit) ? nothing : (hit["video_file"],hit["frame_index"])
function _syncelig_ambiguous(videos,starts,t,error_s,device,hashes)
    center=_syncelig_frame(videos,starts,t,device,hashes)
    isnothing(center) && return true
    any(d->_syncelig_frame_key(_syncelig_frame(videos,starts,t+d,device,hashes))!=_syncelig_frame_key(center),(-error_s,error_s))
end
function _syncelig_uncertainty(fit,t)
    span=max(fit["anchor_max_time_s"]-fit["anchor_min_time_s"],1.0)
    outside=max(fit["anchor_min_time_s"]-t,t-fit["anchor_max_time_s"],0.0)
    fit["max_residual_s"]+1/48000+fit["template_scale_bias_bound_s"]+
        outside*max(fit["clock_scale_uncertainty"],2*fit["max_residual_s"]/span,1/(48000span))
end

function _syncelig_counterfactual(events,manifest,sync,metadata)
    get(sync,"status","")=="valid" || error("Timing eligibility requires saved valid waveform synchronization")
    videos=metadata["videos"];isempty(videos) && error("Missing source chapter metadata")
    stored=Float64[v["global_start_s"] for v in videos];exact=zeros(length(videos));total=0.0
    bases=Any[]
    for (i,v) in enumerate(videos)
        duration=Float64(v["audio_duration_s"])
        isfinite(duration) && duration>0 || error("Invalid declared camera audio duration")
        pts=v["pts"];!isempty(pts) && all(isfinite,pts) && issorted(pts) || error("Invalid original source PTS")
        exact[i]=total;total+=duration
        push!(bases,Dict("video_file"=>v["path"],"stored_start_s"=>stored[i],"exact_declared_start_s"=>exact[i],
            "stored_minus_exact_s"=>stored[i]-exact[i],"declared_audio_duration_s"=>duration))
    end
    issorted(stored) || error("Unsorted saved camera chapter starts")
    path_index=Dict(v["path"]=>i for (i,v) in enumerate(videos))
    anchors=[a for a in sync["anchors"] if get(a,"accepted",false) && get(a,"fit_inlier",false)]
    length(anchors)>=4 || error("Insufficient accepted/inlier waveform anchors")
    x=Float64[a["recording_time_s"] for a in anchors]
    y=Float64[a["video_global_time_s"]-(stored[path_index[a["video_file"]]]-exact[path_index[a["video_file"]]]) for a in anchors]
    all(isfinite,x) && all(isfinite,y) || error("Nonfinite waveform anchors")
    denom=sum(abs2,x.-mean(x));denom>0 || error("No waveform anchor time span")
    scale=sum((x.-mean(x)).*(y.-mean(y)))/denom;offset=mean(y)-scale*mean(x)
    residual=y.-offset.-scale.*x;maxres=maximum(abs,residual);span=max(maximum(x)-minimum(x),1.0)
    applied=Float64[a["template_clock_scale"] for a in anchors]
    # No waveform rerun. Retain unresolved scale relative to every actually used
    # template, rather than letting a common start-anchor bias vanish in a fit.
    scale_error=max(maximum(abs.(scale.-applied)),2*(maxres+1/48000)/span)
    window=maximum(Float64[a["window_s"] for a in anchors])
    fit=Dict{String,Any}("offset_s"=>offset,"scale"=>scale,"drift_ppm"=>(scale-1)*1e6,
        "max_residual_s"=>maxres,"rms_residual_s"=>sqrt(mean(abs2,residual)),
        "clock_scale_uncertainty"=>scale_error,"template_scale_bias_bound_s"=>window*scale_error/2,
        "anchor_min_time_s"=>minimum(x),"anchor_max_time_s"=>maximum(x),"anchor_count"=>length(anchors))
    limit=manifest["request"]["config"]["sync_residual_limit_s"]
    maxres<=limit && abs(scale-1)*1e6<=2000 || error("Exact-base refit fails saved timing bounds; requires independent review")
    hashes=Dict(q["path"]=>q["sha256"] for q in manifest["request"]["inputs"] if haskey(path_index,q["path"]))
    length(hashes)==length(videos) || error("Missing authenticated source video SHA")
    additional=Set{String}();rows=Any[];preserved=String[]
    for e in eachrow(events)
        id=String(e.event_id);oldstatus=String(e.sync_status);t=Float64(e.recording_time_s);device=String(e.device)
        isfinite(t) || error("Nonfinite event recording time")
        oldt=sync["offset_s"]+sync["scale"]*t;newt=offset+scale*t
        oldhit=_syncelig_frame(videos,stored,oldt,device,hashes);newhit=_syncelig_frame(videos,exact,newt,device,hashes)
        stored_error=ismissing(e.sync_error_s) ? nothing : Float64(e.sync_error_s)
        olderror=isnothing(stored_error) ? _syncelig_uncertainty(sync,t) : stored_error
        isfinite(olderror) && olderror>=0 || error("Invalid saved event timing uncertainty")
        computed_error=_syncelig_uncertainty(fit,t)
        # Union uncertainty conservatively: never shrink any saved event bound.
        error_s=max(olderror,computed_error)
        oldkey=(ismissing(e.video_file) || ismissing(e.frame_index)) ? nothing : (String(e.video_file),Int(e.frame_index))
        oldstatus=="mapped" && _syncelig_frame_key(oldhit)!=oldkey && error("Saved mapped event disagrees with original affine/PTS mapping: $id")
        frame_changed=_syncelig_frame_key(newhit)!=oldkey
        newstatus=isnothing(newhit) ? "outside_video" : _syncelig_ambiguous(videos,exact,newt,error_s,device,hashes) ? "ambiguous" : "mapped"
        added=oldstatus=="mapped" && (frame_changed || newstatus!="mapped")
        added && push!(additional,id)
        oldstatus!="mapped" && push!(preserved,id)
        reasons=String[]
        added && frame_changed && push!(reasons,"exact_declared_base_refit_changes_source_frame")
        added && newstatus!="mapped" && push!(reasons,"exact_declared_base_refit_uncertainty_crosses_frame_or_video_boundary")
        push!(rows,Dict("event_id"=>id,"old_sync_status"=>oldstatus,"counterfactual_sync_status"=>newstatus,
            "final_sync_status"=>added ? "ambiguous" : oldstatus,"additional_exclusion"=>added,"reasons"=>reasons,
            "old_recording_time_s"=>t,"old_video_global_time_s"=>oldt,"counterfactual_video_global_time_s"=>newt,
            "old_mapping"=>Dict("video_file"=>e.video_file,"frame_index"=>e.frame_index,"frame_pts_s"=>e.frame_pts_s,"frame_id"=>e.frame_id),
            "old_model_mapping"=>oldhit,"counterfactual_mapping"=>newhit,"central_frame_changed"=>frame_changed,
            "old_stored_error_s"=>stored_error,"old_error_s"=>olderror,"counterfactual_computed_error_s"=>computed_error,"union_error_s"=>error_s))
    end
    Dict{String,Any}("chapter_bases"=>bases,"counterfactual_fit"=>fit,"events"=>rows,
        "additional_excluded_event_ids"=>sort(collect(additional)),"preserved_original_excluded_event_ids"=>preserved),additional
end

function _syncelig_verify_cache(dir,events,manifest)
    path=joinpath(dir,"event_snippets.f64");cache=manifest["snippet_cache"]
    cache["dtype"]=="Float64" && cache["endianness"]=="little" && cache["layout"]=="column-major" || error("Unsupported snippet cache")
    entry=_syncelig_output_entry(manifest,path);entry["sha256"]==cache["sha256"] && entry["bytes"]==cache["bytes"] || error("Inconsistent cache identities")
    ident=_syncelig_check_identity(path,entry);bytes=read(path);ranges=Tuple{Int,Int}[]
    for e in eachrow(events)
        String(e.snippet_file)=="event_snippets.f64" || error("Unexpected snippet file")
        a=Int(e.snippet_byte_offset);n=Int(e.snippet_samples);channels=Int(e.snippet_channels);peak=Int(e.snippet_peak_index)
        a>=0 && a%8==0 && 1<=n<=1001 && channels==3 && 0<=peak<n || error("Invalid snippet range/shape")
        b=a+n*channels*8;b<=length(bytes) || error("Snippet range outside cache")
        _syncelig_digest(bytes[a+1:b])==e.snippet_sha256 || error("Snippet payload SHA mismatch: $(e.event_id)")
        sum(q["samples"] for q in JSON.parse(String(e.window_source_segments)))==n || error("Snippet source-window count differs")
        push!(ranges,(a,b))
    end
    sort!(ranges)
    if isempty(ranges)
        isempty(bytes) || error("Unreferenced snippet cache bytes")
    else
        ranges[1][1]==0 && ranges[end][2]==length(bytes) || error("Snippet cache coverage incomplete")
        all(i->ranges[i][2]==ranges[i+1][1],1:length(ranges)-1) || error("Snippet cache gaps/overlaps")
    end
    sha256_file(path)==ident["sha256"] || error("Snippet cache changed during timing review")
    ident
end

"""
    review_sync_eligibility(output_dir;resume=true) -> Dict

Postprocess frozen preparation using only previously saved waveform anchors and
exact sums of declared camera audio durations. The original waveform fit is not
rerun and this is not independent physical clock truth. All original nonmapped
statuses stay excluded. Original mapped rows receive an additional ambiguity
exclusion when the counterfactual changes their source frame or its conservative
uncertainty interval crosses a frame/video boundary. No labels, image positions,
TDOAs, annotation truth or fit errors enter this decision.

Writes sync_eligibility.json with per-event mappings/errors and source/input SHA.
With zero additions, every canonical producer output (including its manifest)
remains byte-identical. With additions, byte-exact once-only originals are kept
in preparation_raw/, only selected sync_status CSV cells change, and the new
manifest records explicit postprocessing lineage and the new events identity.
Images/cache/frames and original request/input/source binding are never rewritten.
Resume authenticates immutable originals and recognizes an interrupted two-file
commit; arbitrary canonical edits are rejected, never adopted as a new baseline.
"""
function review_sync_eligibility(output_dir;resume=true)
    try
        _review_sync_eligibility(abspath(output_dir);resume=resume)
    catch exception
        log_failure("review_sync_eligibility output=$output_dir",exception,catch_backtrace())
        rethrow()
    end
end
function _review_sync_eligibility(dir;resume)
    eventpath=joinpath(dir,"events.csv");manifestpath=joinpath(dir,"prepare_manifest.json")
    qcpath=joinpath(dir,"sync_eligibility.json");rawdir=joinpath(dir,"preparation_raw")
    rawevent=joinpath(rawdir,"events.csv");rawmanifest=joinpath(rawdir,"prepare_manifest.json")
    islink(rawdir) && error("Raw preparation directory must not be a symlink")
    currenteventbytes=read(eventpath);currentmanifestbytes=read(manifestpath)
    raw_eventbytes=isfile(rawevent) ? read(rawevent) : currenteventbytes
    raw_manifestbytes=isfile(rawmanifest) ? read(rawmanifest) : currentmanifestbytes
    raw=JSON.parse(String(copy(raw_manifestbytes)))
    any(q->get(q,"stage","")=="sync_eligibility",get(raw,"postprocessing",Any[])) && error("Missing immutable original producer manifest")
    get(raw,"status","")=="complete" || error("Timing eligibility requires a complete preparation")
    original_event_sha=_syncelig_digest(raw_eventbytes);original_manifest_sha=_syncelig_digest(raw_manifestbytes)
    original_entry=_syncelig_output_entry(raw,eventpath)
    original_entry["sha256"]==original_event_sha && original_entry["bytes"]==length(raw_eventbytes) || error("Raw events do not match original producer manifest")
    events=CSV.read(IOBuffer(copy(raw_eventbytes)),DataFrame)
    length(unique(events.event_id))==nrow(events) || error("Duplicate event IDs")
    raw["event_count"]==nrow(events) || error("Original event count differs")
    syncpath=joinpath(dir,"synchronization.json");metapath=joinpath(dir,"prepare_metadata.json");framepath=joinpath(dir,"frames.csv")
    identities=Dict{String,Any}()
    for (name,path) in (("synchronization",syncpath),("metadata",metapath),("frames",framepath))
        identities[name]=_syncelig_check_identity(path,_syncelig_output_entry(raw,path))
    end
    identities["snippet_cache"]=_syncelig_verify_cache(dir,events,raw)
    sync=read_json(syncpath);metadata=read_json(metapath)
    review,additional=_syncelig_counterfactual(events,raw,sync,metadata)
    patched=_syncelig_patch_status(raw_eventbytes,events,additional)
    changed_events=_syncelig_check_table_change(raw_eventbytes,patched,additional)
    target_event_sha=_syncelig_digest(patched);source_sha=sha256_file(@__FILE__)
    current_event_sha=_syncelig_digest(currenteventbytes);current_manifest_sha=_syncelig_digest(currentmanifestbytes)
    current_event_sha in (original_event_sha,target_event_sha) || error("Canonical events contain edits outside timing-review changes")
    lineage=nothing
    if current_manifest_sha!=original_manifest_sha
        current=JSON.parse(String(copy(currentmanifestbytes)))
        candidates=[q for q in get(current,"postprocessing",Any[]) if get(q,"stage","")=="sync_eligibility"]
        length(candidates)==1 || error("Canonical manifest has no unique timing-review lineage")
        lineage=only(candidates)
        for (k,v) in (("raw_manifest_sha256",original_manifest_sha),("raw_events_sha256",original_event_sha),
            ("review_source_sha256",source_sha),("events_sha256_after",target_event_sha))
            get(lineage,k,nothing)==v || error("Incompatible timing-review lineage: $k")
        end
        current["request"]==raw["request"] || error("Original request binding changed")
        current["mapped_count"]==count(==("mapped"),changed_events.sync_status) || error("Canonical mapped count differs from reviewed table")
        get(current,"producer_manifest_sha256",nothing)==original_manifest_sha || error("Original producer manifest binding differs")
        get(lineage,"additional_excluded_event_ids",nothing)==sort(collect(additional)) || error("Timing exclusion lineage differs")
        for key in keys(raw)
            key in ("outputs","mapped_count","postprocessing","producer_manifest_sha256") && continue
            get(current,key,nothing)==raw[key] || error("Forbidden producer manifest field changed: $key")
        end
        _syncelig_non_event_outputs(current,eventpath)==_syncelig_non_event_outputs(raw,eventpath) || error("Non-event output inventory changed")
        _syncelig_output_entry(current,eventpath)["sha256"]==target_event_sha || error("Canonical event identity not committed in manifest")
    end
    input_binding=Dict("raw_manifest_sha256"=>original_manifest_sha,"raw_events_sha256"=>original_event_sha,
        "producer_source_sha256"=>raw["request"]["source_sha256"],"review_source_sha256"=>source_sha,
        "saved_inputs"=>identities)
    if resume && isfile(qcpath)
        prior=read_json(qcpath)
        if get(prior,"status","")=="complete" && get(prior,"input_binding",nothing)==input_binding &&
            get(prior,"events_sha256_after",nothing)==current_event_sha && get(prior,"manifest_sha256_after",nothing)==current_manifest_sha
            return prior
        end
    end
    nowutc=string(now(UTC))*"Z"
    qc=Dict{String,Any}("schema"=>SYNC_ELIGIBILITY_SCHEMA,"status"=>"planned","review_source_sha256"=>source_sha,
        "input_binding"=>input_binding,"created_utc"=>nowutc,"algorithm"=>"saved waveform anchors; exact declared camera audio bases; affine refit; union timing exclusions",
        "limitations"=>["No new waveform evidence or independent physical clock truth; declared audio durations retain source metadata precision.",
            "No waveform re-warp is performed; unresolved difference from the actually applied template scales is included in timing uncertainty.",
            "Original error bounds never shrink; all original nonmapped events stay excluded; no label, acoustic prediction or fit-error-driven exclusion.",
            "A changed counterfactual source frame is excluded; original frame/image/cache identities are retained and never re-extracted."],
        "review"=>review,"additional_excluded_event_ids"=>sort(collect(additional)),"additional_exclusion_count"=>length(additional),
        "original_mapped_count"=>count(==("mapped"),events.sync_status),"mapped_count_after"=>count(==("mapped"),changed_events.sync_status),
        "events_sha256_before"=>original_event_sha,"events_sha256_after"=>target_event_sha,"manifest_sha256_before"=>original_manifest_sha,
        "canonical_changes_applied"=>!isempty(additional))
    qc=json_safe(qc)
    if isempty(additional)
        current_event_sha==original_event_sha && current_manifest_sha==original_manifest_sha || error("Zero-addition review cannot adopt modified canonical files")
        qc["status"]="complete";qc["manifest_sha256_after"]=original_manifest_sha
        _syncelig_write_json(qcpath,qc)
        return qc
    end
    # Copy exact raw bytes before any canonical write; existing copies are never
    # overwritten. Partial raw preservation can recover only from authenticated
    # originals, as checked above.
    _syncelig_write_bytes(rawevent,raw_eventbytes;preserve=true)
    _syncelig_write_bytes(rawmanifest,raw_manifestbytes;preserve=true)
    qc["raw_events_path"]=rawevent;qc["raw_manifest_path"]=rawmanifest
    _syncelig_write_json(qcpath,qc)
    if current_event_sha!=target_event_sha
        _syncelig_write_bytes(eventpath,patched)
    end
    sha256_file(eventpath)==target_event_sha || error("Canonical event commit digest differs")
    if isnothing(lineage)
        updated=deepcopy(raw)
        entry=_syncelig_output_entry(updated,eventpath);merge!(entry,_syncelig_identity(eventpath))
        updated["mapped_count"]=qc["mapped_count_after"]
        lineage=Dict{String,Any}("stage"=>"sync_eligibility","schema"=>SYNC_ELIGIBILITY_SCHEMA,"applied_utc"=>nowutc,
            "raw_manifest_path"=>rawmanifest,"raw_manifest_sha256"=>original_manifest_sha,
            "raw_events_path"=>rawevent,"raw_events_sha256"=>original_event_sha,"review_source_sha256"=>source_sha,
            "synchronization_sha256"=>identities["synchronization"]["sha256"],"metadata_sha256"=>identities["metadata"]["sha256"],
            "qc_path"=>qcpath,"events_sha256_after"=>target_event_sha,"additional_excluded_event_ids"=>sort(collect(additional)),
            "exclusion_reason"=>"Conservative union: original mapped event changes source frame or crosses frame/video uncertainty under exact declared camera bases",
            "original_mapped_count"=>qc["original_mapped_count"],"mapped_count_after"=>qc["mapped_count_after"])
        history=get!(updated,"postprocessing",Any[]);push!(history,lineage)
        updated["producer_manifest_sha256"]=original_manifest_sha
        updated["request"]==raw["request"] || error("Original producer request changed")
        _syncelig_non_event_outputs(updated,eventpath)==_syncelig_non_event_outputs(raw,eventpath) || error("Non-event outputs changed")
        _syncelig_write_json(manifestpath,updated)
    end
    qc["status"]="complete";qc["manifest_sha256_after"]=sha256_file(manifestpath)
    _syncelig_write_json(qcpath,qc)
    qc
end
