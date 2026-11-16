# Schema clicker-physical-calibration/2.0: type, loader, lens-aware projection (pinhole), runtime sound speed.
# Included inside ClickerCalibration after workflow.jl/physical_fit.jl. Definitions only; loading does no I/O.
# New code only: the v1 functions (load_calibration, _pf_project_direction, project_tdoas, ...) are untouched.
# See docs/clicker_calibration_v2_contract.md.

# Warm seawater at the site; overridable per deployment (the profile only records the value used for fitting).
const CLICKER_SITE_SOUND_SPEED_M_S = 1530.0   # identical re-definition is a no-op, so repeated includes are safe

const CLICKER_V2_SCHEMA = "clicker-physical-calibration/2.0"
const CLICKER_V2_FAMILIES = ("pinhole", "fisheye", "superview")

"""
    CalibrationV2 <: AbstractDict{String,Any}

A validated `clicker-physical-calibration/2.0` profile. It forwards the dictionary interface to the wrapped
`Dict{String,Any}`, so `project_tdoas`, `project_angles` and every other `AbstractDict` consumer work unchanged.
Build one with [`load_calibration_2`](@ref).
"""
struct CalibrationV2 <: AbstractDict{String,Any}
    data::Dict{String,Any}
end
CalibrationV2(d::AbstractDict) = CalibrationV2(Dict{String,Any}(string(k) => v for (k, v) in d))

Base.getindex(c::CalibrationV2, key) = getindex(c.data, key)
Base.get(c::CalibrationV2, key, default) = get(c.data, key, default)
Base.get(f::Base.Callable, c::CalibrationV2, key) = get(f, c.data, key)
Base.haskey(c::CalibrationV2, key) = haskey(c.data, key)
Base.keys(c::CalibrationV2) = keys(c.data)
Base.values(c::CalibrationV2) = values(c.data)
Base.length(c::CalibrationV2) = length(c.data)
Base.iterate(c::CalibrationV2) = iterate(c.data)
Base.iterate(c::CalibrationV2, state) = iterate(c.data, state)
Base.setindex!(c::CalibrationV2, value, key) = (setindex!(c.data, value, key); c)
Base.delete!(c::CalibrationV2, key) = (delete!(c.data, key); c)
Base.copy(c::CalibrationV2) = CalibrationV2(copy(c.data))
Base.show(io::IO, c::CalibrationV2) = print(io, "CalibrationV2(device=", get(c.data, "device", "?"),
    ", lens=", get(get(c.data, "lens", Dict()), "family", "?"), ", c=", get(c.data, "sound_speed_m_s", "?"), " m/s)")

_v2_fail(msg) = error("Invalid v2 calibration: " * msg)
_v2_require(cond, msg) = cond || _v2_fail(msg)
_v2_isnum(x) = x isa Real && !(x isa Bool) && isfinite(x)

function _v2_key(d, key, what)
    d isa AbstractDict && haskey(d, key) || _v2_fail("missing $(what)$(key)")
    return d[key]
end

function _v2_check_pinholelike(p, what)
    p isa AbstractDict || _v2_fail("$what must be an object")
    for key in ("fx", "fy", "cx", "cy")
        _v2_require(_v2_isnum(_v2_key(p, key, what * ".")), "$what.$key must be finite")
    end
    _v2_require(p["fx"] > 0 && p["fy"] > 0, "$what focal lengths must be positive")
end

function _v2_check_lens(lens)
    lens isa AbstractDict || _v2_fail("lens must be an object")
    family = _v2_key(lens, "family", "lens.")
    _v2_require(family in CLICKER_V2_FAMILIES, "lens.family must be one of $(join(CLICKER_V2_FAMILIES, ", "))")
    params = _v2_key(lens, "params", "lens.")
    _v2_require(params isa AbstractDict, "lens.params must be an object")
    source = _v2_key(lens, "source", "lens.")
    _v2_require(source isa AbstractString && !isempty(source), "lens.source must be a non-empty string")
    if family == "fisheye"
        _v2_check_pinholelike(params, "lens.params")
        for key in ("k1", "k2", "k3", "k4")
            _v2_require(_v2_isnum(_v2_key(params, key, "lens.params.")), "lens.params.$key must be finite")
        end
    elseif family == "superview"
        inner_family = _v2_key(params, "inner_family", "lens.params.")
        _v2_require(inner_family in ("pinhole", "fisheye"), "lens.params.inner_family must be pinhole or fisheye")
        inner = _v2_key(params, "inner", "lens.params.")
        _v2_check_pinholelike(inner, "lens.params.inner")
        inner_family == "fisheye" && for key in ("k1", "k2", "k3", "k4")
            _v2_require(_v2_isnum(_v2_key(inner, key, "lens.params.inner.")), "lens.params.inner.$key must be finite")
        end
        _v2_require(_v2_isnum(_v2_key(params, "sensor_width_px", "lens.params.")) && params["sensor_width_px"] > 0,
            "lens.params.sensor_width_px must be positive")
        coeffs = _v2_key(params, "remap_coeffs", "lens.params.")
        _v2_require(coeffs isa AbstractVector && !isempty(coeffs) && all(_v2_isnum, coeffs),
            "lens.params.remap_coeffs must be a non-empty finite vector")
    end
    return Symbol(family)
end

_v2_date(x, key) = x === nothing ? nothing :
    (x isa AbstractString ? (try Date(x) catch; _v2_fail("$key must be an ISO date or null") end) :
     _v2_fail("$key must be an ISO date or null"))

"""
    validate_calibration_2(raw; device=nothing) -> CalibrationV2

Check every rule of the schema table in docs/clicker_calibration_v2_contract.md; throws `ErrorException`
with a clear message on the first failure. `raw` is a parsed JSON object.
"""
function validate_calibration_2(raw; device=nothing)
    raw isa AbstractDict || _v2_fail("top level must be an object")
    cal = Dict{String,Any}(string(k) => v for (k, v) in raw)
    _v2_require(get(cal, "schema_version", nothing) == CLICKER_V2_SCHEMA, "unsupported schema (need $CLICKER_V2_SCHEMA)")
    dev = _v2_key(cal, "device", "")
    _v2_require(dev in ("A", "B"), "device must be \"A\" or \"B\"")
    isnothing(device) || string(dev) == string(device) || _v2_fail("calibration belongs to device $dev, not $device")
    _v2_require(_v2_key(cal, "method", "") == "agnostic", "method must be \"agnostic\"")
    spec = _v2_key(cal, "estimator_spec", "")
    _v2_require(spec isa AbstractDict && get(spec, "method", nothing) == "agnostic" &&
        _v2_isnum(get(spec, "sample_rate_hz", nothing)) && spec["sample_rate_hz"] > 0,
        "estimator_spec must be {method:\"agnostic\", sample_rate_hz>0}")
    _v2_require(_v2_key(cal, "physical_channels", "") == [3, 4, 5], "unexpected physical channel mapping")
    _v2_require(_v2_key(cal, "reference_channel", "") == 3, "reference_channel must be 3")
    _v2_require(_v2_key(cal, "reference_physical_channel", "") == 5, "reference_physical_channel must be 5")
    geometry = try _pf_matrix(_v2_key(cal, "geometry_m", "")) catch; _v2_fail("geometry_m must be a 3×3 matrix") end
    _v2_require(size(geometry) == (3, 3) && all(isfinite, geometry) && maximum(abs.(geometry .- clicker_geometry())) < 1e-12,
        "geometry_m is not the confirmed D3 array")
    offsets = _v2_key(cal, "relative_delay_offsets_s", "")
    _v2_require(offsets isa AbstractVector && length(offsets) == 3 && all(_v2_isnum, offsets) && all(==(0), offsets),
        "relative_delay_offsets_s must be exactly [0,0,0] in v2")
    speed = _v2_key(cal, "sound_speed_m_s", "")
    _v2_require(_v2_isnum(speed) && 1400 <= speed <= 1600, "sound_speed_m_s outside 1400..1600")
    _v2_require(_v2_key(cal, "sound_speed_role", "") == "runtime_environment", "sound_speed_role must be \"runtime_environment\"")
    rotation = try _pf_matrix(_v2_key(cal, "rotation_array_to_camera", "")) catch; _v2_fail("rotation_array_to_camera must be 3×3") end
    _v2_require(size(rotation) == (3, 3) && all(isfinite, rotation) && norm(rotation' * rotation - I) < 1e-9 && abs(det(rotation) - 1) < 1e-9,
        "rotation must be proper SO(3)")
    euler = _v2_key(cal, "rotation_euler_xyz_rad", "")
    _v2_require(euler isa AbstractVector && length(euler) == 3 && all(_v2_isnum, euler), "rotation_euler_xyz_rad must be 3 finite numbers")
    _v2_require(norm(_pf_rotation(Float64.(euler)) - rotation) < 1e-6, "rotation_euler_xyz_rad disagrees with rotation_array_to_camera (Rz·Ry·Rx)")
    for key in ("width", "height")
        value = _v2_key(cal, key, "")
        _v2_require(value isa Integer && !(value isa Bool) && value > 0, "$key must be a positive integer")
    end
    camera = _v2_key(cal, "camera", "")
    camera isa AbstractDict || _v2_fail("camera must be an object")
    for key in ("fx", "fy", "cx", "cy", "k1", "k2")
        _v2_require(_v2_isnum(_v2_key(camera, key, "camera.")), "camera.$key must be finite")
    end
    _v2_require(camera["fx"] > 0 && camera["fy"] > 0, "camera focal lengths must be positive")
    _v2_check_lens(_v2_key(cal, "lens", ""))
    mode = _v2_key(cal, "lens_mode", "")
    _v2_require(mode === nothing || mode isa AbstractString, "lens_mode must be a string or null")
    zfov = _v2_key(cal, "zfov_deg", "")
    _v2_require(zfov === nothing || (_v2_isnum(zfov) && zfov > 0), "zfov_deg must be a positive number or null")
    identity = _v2_key(cal, "camera_identity", "")
    _v2_require(identity isa AbstractDict && all(k -> get(identity, k, nothing) === nothing || get(identity, k, nothing) isa AbstractString, ("CAME", "CASN", "LINF")),
        "camera_identity must be an object whose CAME, CASN and LINF are strings or null")
    chain = _v2_key(cal, "camera_chain", "")
    _v2_require(chain isa AbstractString && !isempty(chain), "camera_chain must be a non-empty string")
    from = _v2_date(_v2_key(cal, "valid_from", ""), "valid_from")
    until = _v2_date(_v2_key(cal, "valid_until", ""), "valid_until")
    (from === nothing || until === nothing || from <= until) || _v2_fail("valid_from is after valid_until")
    _v2_require(_v2_key(cal, "model_ladder_choice", "") in ("M0", "M1", "M2", "M3"), "model_ladder_choice must be M0..M3")
    _v2_require(_v2_key(cal, "metrics", "") isa AbstractDict, "metrics must be an object")
    _v2_require(_v2_key(cal, "provenance", "") isa AbstractDict, "provenance must be an object")
    return CalibrationV2(cal)
end

"""
    load_calibration_2(path; device=nothing) -> CalibrationV2

Read and validate a schema-2.0 profile. `device` ("A"/"B"/:A/:B) must match the profile when given.
Throws `ErrorException` (clear message) on any schema violation; nothing else is touched.
"""
load_calibration_2(path::AbstractString; device=nothing) = validate_calibration_2(read_json(path); device=device)

"""Copy of `calibration` whose `sound_speed_m_s` is the runtime `c_m_s` (1400..1600); the fitted value is kept as `sound_speed_fit_m_s`."""
function _v2_with_sound_speed(calibration::CalibrationV2, c_m_s)
    _v2_isnum(c_m_s) && 1400 <= c_m_s <= 1600 || throw(ArgumentError("c_m_s must be finite and within 1400..1600 m/s"))
    copied = CalibrationV2(deepcopy(calibration.data))
    copied["sound_speed_fit_m_s"] = calibration["sound_speed_m_s"]
    copied["sound_speed_m_s"] = Float64(c_m_s)
    return copied
end

_v2_family(calibration) = String(get(get(calibration, "lens", Dict{String,Any}()), "family", "pinhole"))

# ---- lens families -------------------------------------------------------------------------------------
# _v2_lens_project(Val(family), params, xn, yn) -> (x_px, y_px) or `nothing` for an invalid projection.
#   xn=X/Z, yn=Y/Z are camera-frame, +y up; pixels are zero-origin, +y down (as v1).
# _v2_lens_unproject(Val(family), params, x_px, y_px) -> (xn, yn), the inverse.
# pinhole: `params` is a camera-like object with fx,fy,cx,cy,k1,k2 (pass calibration["camera"]).
# fisheye/superview methods belong to constrained_fit.jl (T3).
_v2_lens_project(::Val{F}, params, xn, yn) where {F} = throw(ArgumentError("lens family not implemented yet"))
_v2_lens_unproject(::Val{F}, params, x_px, y_px) where {F} = throw(ArgumentError("lens family not implemented yet"))

function _v2_lens_project(::Val{:pinhole}, camera, xn, yn)
    r2 = xn * xn + yn * yn
    k1 = get(camera, "k1", 0.0); k2 = get(camera, "k2", 0.0)
    radial = 1.0 + k1 * r2 + k2 * r2 * r2
    derivative = 1.0 + 3k1 * r2 + 5k2 * r2 * r2
    radial > 0 && derivative > 0 || return nothing
    x = Float64(camera["cx"]) + Float64(camera["fx"]) * xn * radial
    y = Float64(camera["cy"]) - Float64(camera["fy"]) * yn * radial
    return isfinite(x) && isfinite(y) ? (x, y) : nothing
end

function _v2_lens_unproject(::Val{:pinhole}, camera, x_px, y_px)
    xd = (x_px - Float64(camera["cx"])) / Float64(camera["fx"])
    yd = -(y_px - Float64(camera["cy"])) / Float64(camera["fy"])
    k1 = Float64(get(camera, "k1", 0.0)); k2 = Float64(get(camera, "k2", 0.0))
    rd = hypot(xd, yd)
    rd == 0 && return (0.0, 0.0)
    r = rd
    for _ in 1:50   # Newton on r*(1+k1 r²+k2 r⁴) = rd
        r2 = r * r
        f = r * (1 + k1 * r2 + k2 * r2 * r2) - rd
        df = 1 + 3k1 * r2 + 5k2 * r2 * r2
        df == 0 && break
        step = f / df
        r -= step
        abs(step) < 1e-15 * max(1.0, abs(r)) && break
    end
    scale = r / rd
    return (xd * scale, yd * scale)
end

"""
Lens-aware ray projection for `CalibrationV2`. Pinhole reproduces the untyped v1 method exactly (it calls it).
Other families project the camera ray through `_v2_lens_project(Val(family), lens.params, X/Z, Y/Z)`.
"""
function _pf_project_direction(calibration::CalibrationV2, direction)
    family = _v2_family(calibration)
    family == "pinhole" && return invoke(_pf_project_direction, Tuple{Any,Any}, calibration, direction)
    u = Float64.(direction)
    all(isfinite, u) || return _pf_invalid("nonfinite_direction")
    camera_direction = _pf_matrix(calibration["rotation_array_to_camera"]) * u
    X, Y, Z = camera_direction
    Z > 1e-9 || return _pf_invalid("behind_camera"; direction_array=u, direction_camera=camera_direction)
    xy = _v2_lens_project(Val(Symbol(family)), calibration["lens"]["params"], X / Z, Y / Z)
    xy === nothing && return _pf_invalid("invalid_distortion"; direction_array=u, direction_camera=camera_direction)
    x, y = Float64(xy[1]), Float64(xy[2])
    isfinite(x) && isfinite(y) || return _pf_invalid("nonfinite_projection"; direction_array=u, direction_camera=camera_direction)
    inside = 0 <= x < calibration["width"] && 0 <= y < calibration["height"]
    return (x=x, y=y, valid=true, status=inside ? "ok" : "outside_image",
        direction_array=u, direction_camera=camera_direction, in_frame=inside)
end

# ---- per-video checks used by run_contiguous_folders_calibrated -----------------------------------------------

"""
    _v2_gpmf_lens_infos(videos; python=ENV CLICKER_PYTHON or "python3") -> Vector{Union{Dict,Nothing}}

Run `python3 -I tools/gpmf_lens_info.py VIDEO...` once (it prints one JSON line per file: VFOV, ZFOV, CASN, CAME, ...).
Entry `i` is the parsed record of `videos[i]`, or `nothing` (with an `@warn`) when the tool is missing or fails.
"""
function _v2_gpmf_lens_infos(videos::AbstractVector{<:AbstractString}; python=get(ENV, "CLICKER_PYTHON", "python3"))
    none = Union{Dict{String,Any},Nothing}[nothing for _ in videos]
    tool = joinpath(@__DIR__, "..", "..", "tools", "gpmf_lens_info.py")
    isfile(tool) || (@warn "GPMF lens tool not found; lens identity unchecked" tool; return none)
    isempty(videos) && return none
    try
        out = read(pipeline(`$python -I $tool $(collect(String, videos))`; stderr=devnull), String)
        records = [JSON.parse(String(l)) for l in split(out, '\n') if startswith(strip(l), "{")]
        length(records) == length(videos) || (@warn "GPMF lens tool returned $(length(records)) records for $(length(videos)) videos; lens identity unchecked"; return none)
        return Union{Dict{String,Any},Nothing}[Dict{String,Any}(r) for r in records]
    catch err
        @warn "GPMF lens tool failed; lens identity unchecked" exception=(err, catch_backtrace())
        return none
    end
end
_v2_gpmf_lens_info(video::AbstractString; kwargs...) = only(_v2_gpmf_lens_infos([String(video)]; kwargs...))

_v2_str(x) = x === nothing ? nothing : string(x)

"""
    _v2_validate_video(calibration; width, height, lens_infos=nothing, recording_date=nothing) -> NamedTuple

Per-recording checks. Errors (`ErrorException`): resolution differs from the profile; for any chapter, when both
sides are known, GPMF `VFOV` differs from `lens_mode`, `ZFOV` differs from `zfov_deg` by more than 0.01°, or `CASN`
differs from `camera_identity.CASN`. `lens_infos` holds one parsed `gpmf_lens_info.py` record (or `nothing`) per
chapter. Warnings (`@warn`): lens identity unknown, recording date outside `valid_from..valid_until`.
Returns `(warnings, lens, in_validity_window)`; `lens` has the observed `lens_mode`, `zfov_deg`, `casn`
(unique values over the chapters) and the profile's expected ones.
"""
function _v2_validate_video(calibration::AbstractDict; width, height, lens_infos=nothing, recording_date=nothing)
    warnings = String[]
    (Int(width), Int(height)) == (Int(calibration["width"]), Int(calibration["height"])) ||
        error("Video resolution $(width)x$(height) differs from the calibration profile $(calibration["width"])x$(calibration["height"]) (device $(calibration["device"])); refusing to run")
    want_mode = get(calibration, "lens_mode", nothing)
    want_zfov = get(calibration, "zfov_deg", nothing)
    want_casn = get(get(calibration, "camera_identity", Dict()), "CASN", nothing)
    modes = Any[]; zfovs = Any[]; casns = Any[]
    unknown = false
    for (i, info) in enumerate(something(lens_infos, Any[]))
        vfov = info === nothing ? nothing : _v2_str(get(info, "VFOV", nothing))
        zfov = info === nothing ? nothing : get(info, "ZFOV", nothing)
        casn = info === nothing ? nothing : _v2_str(get(info, "CASN", nothing))
        (vfov === nothing || zfov === nothing) && (unknown = true)
        vfov !== nothing && (push!(modes, vfov); want_mode !== nothing && string(want_mode) != vfov &&
            error("Chapter $i GPMF lens mode (VFOV) \"$vfov\" differs from the calibration profile \"$want_mode\"; the intrinsics differ, refusing to run"))
        zfov isa Real && (push!(zfovs, Float64(zfov)); want_zfov !== nothing && abs(Float64(zfov) - Float64(want_zfov)) > 0.01 &&
            error("Chapter $i GPMF ZFOV $(zfov)° differs from the calibration profile zfov_deg $(want_zfov)° (tolerance 0.01°); the intrinsics differ, refusing to run"))
        casn !== nothing && (push!(casns, casn); want_casn !== nothing && string(want_casn) != casn &&
            error("Chapter $i camera serial (CASN) $casn differs from the calibration profile camera_identity.CASN $want_casn; refusing to run"))
    end
    if lens_infos === nothing || unknown
        push!(warnings, "lens identity (GPMF VFOV/ZFOV) unknown for at least one chapter; the profile's lens could not be fully checked")
    end
    in_window = nothing
    if recording_date !== nothing
        from = _v2_date(get(calibration, "valid_from", nothing), "valid_from")
        until = _v2_date(get(calibration, "valid_until", nothing), "valid_until")
        d = recording_date isa DateTime ? Date(recording_date) : recording_date
        in_window = (from === nothing || d >= from) && (until === nothing || d <= until)
        in_window || push!(warnings, "recording date $d is outside the profile validity window ($(from) .. $(until)); the camera mount may have changed")
    end
    foreach(w -> @warn(w), warnings)
    lens = (lens_mode=unique(modes), zfov_deg=unique(zfovs), casn=unique(casns),
            profile_lens_mode=want_mode, profile_zfov_deg=want_zfov, profile_casn=want_casn)
    return (warnings=warnings, lens=lens, in_validity_window=in_window)
end

"""Recording date of `video`: the GPMF tool's `creation_time` (or `creation_time` argument), else the file's mtime. `nothing` if neither works."""
function _v2_recording_date(video::AbstractString; creation_time=nothing)
    if creation_time isa AbstractString && length(creation_time) >= 10
        d = tryparse(Date, creation_time[1:10])
        d === nothing || return d
    end
    try
        return Date(Dates.unix2datetime(mtime(video)))
    catch
        return nothing
    end
end

"""Current git commit of this repository (`git rev-parse HEAD`), or `nothing` when git or the repository is unavailable."""
function _v2_git_commit()
    try
        root = abspath(joinpath(@__DIR__, "..", ".."))
        commit = strip(read(pipeline(Cmd(`git rev-parse HEAD`; dir=root); stderr=devnull), String))
        return isempty(commit) ? nothing : String(commit)
    catch
        return nothing
    end
end

_v2_now_utc() = Dates.now(Dates.UTC)
