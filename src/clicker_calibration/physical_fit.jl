# Native Julia physical clicker calibration. This file is included inside
# ClickerCalibration after support.jl; loading it performs no I/O or computation.

"""D3 sensor positions in metres: rows are physical channels 3,4,5; columns x,y,z.
The equilateral side is 0.18 m (circumradius 0.18/sqrt(3)); angles -30,90,210°.
Array and camera frames are right handed, x right, y up, z forward.
"""
function clicker_geometry()
    angles = deg2rad.([-30.0, 90.0, 210.0])
    radius = 0.18 / sqrt(3.0)
    return hcat(radius .* cos.(angles), radius .* sin.(angles), zeros(3))
end

_pf_matrix(rows) = rows isa AbstractMatrix ? Matrix{Float64}(rows) : reduce(vcat, [permutedims(Float64.(r)) for r in rows])
_pf_rows(matrix) = [collect(row) for row in eachrow(matrix)]
_pf_finite(x) = x !== missing && x !== nothing && x isa Real && isfinite(x)
_pf_text(x) = x === missing || x === nothing ? "" : string(x)
_pf_number(x, default=NaN) = _pf_finite(x) ? Float64(x) : default

function _pf_invalid(status; direction_array=fill(NaN, 3), direction_camera=fill(NaN, 3))
    return (x=NaN, y=NaN, valid=false, status=string(status),
        direction_array=direction_array, direction_camera=direction_camera, in_frame=false)
end

function _pf_project_direction(calibration, direction)
    u = Float64.(direction)
    all(isfinite, u) || return _pf_invalid("nonfinite_direction")
    camera_direction = _pf_matrix(calibration["rotation_array_to_camera"]) * u
    X, Y, Z = camera_direction
    Z > 1e-9 || return _pf_invalid("behind_camera"; direction_array=u, direction_camera=camera_direction)
    camera = calibration["camera"]
    xn, yn = X / Z, Y / Z
    r2 = xn*xn + yn*yn
    radial = 1.0 + get(camera, "k1", 0.0)*r2 + get(camera, "k2", 0.0)*r2*r2
    derivative = 1.0 + 3get(camera, "k1", 0.0)*r2 + 5get(camera, "k2", 0.0)*r2*r2
    radial > 0 && derivative > 0 || return _pf_invalid("invalid_distortion"; direction_array=u, direction_camera=camera_direction)
    x = Float64(camera["cx"]) + Float64(camera["fx"]) * xn * radial
    y = Float64(camera["cy"]) - Float64(camera["fy"]) * yn * radial
    isfinite(x) && isfinite(y) || return _pf_invalid("nonfinite_projection"; direction_array=u, direction_camera=camera_direction)
    inside = 0 <= x < calibration["width"] && 0 <= y < calibration["height"]
    return (x=x, y=y, valid=true, status=inside ? "ok" : "outside_image",
        direction_array=u, direction_camera=camera_direction, in_frame=inside)
end

"""
    project_tdoas(calibration, delays::AbstractVector; fs=nothing, units="samples")
    project_tdoas(calibration, delays::AbstractMatrix; fs=nothing, units="samples")

Project reference-relative delays (length 2 or 3, or N×2/N×3 Float64-compatible
rows) to original-resolution zero-origin pixels, x right/y down. Two columns
are physical channels 3,4 relative to physical 5. Three columns include the
zero reference; nonzero reference is invalid, not silently re-referenced.
`fs` is positive Hz for sample units (scalar or one per matrix row); seconds
need no fs. Measured delays obey dt_i=-(r_i-r_ref)⋅u/c+b_i. Offsets b_i are
seconds. The planar array uses the fixed +z front branch. Impossible disk
TDOAs are invalid and NEVER clamped, including roundoff beyond the unit disk.
Only front-camera rays project. This is far-field direction-only geometry;
range/near-field curvature is unobservable and unmodelled.

Vector input returns `(x,y,valid,status,direction_array,direction_camera,in_frame)`;
matrix input returns a vector of these records. Invalid pixels are NaN.
"""
function project_tdoas(calibration::AbstractDict, delays::AbstractVector; fs=nothing, units="samples")
    length(delays) in (2,3) || throw(ArgumentError("delays must have 2 or 3 selected-channel entries"))
    units in ("samples", "seconds") || throw(ArgumentError("units must be samples or seconds"))
    if units == "samples"
        _pf_finite(fs) && fs > 0 || throw(ArgumentError("positive finite fs in Hz required for sample delays"))
    end
    all(_pf_finite, delays) || return _pf_invalid("nonfinite_tdoa")
    scale = units == "samples" ? 1.0 / Float64(fs) : 1.0
    dt = Float64.(delays) .* scale
    length(dt) == 3 && abs(dt[3]) > 1e-12 && return _pf_invalid("invalid_reference")
    positions = _pf_matrix(calibration["geometry_m"])
    size(positions) == (3,3) || throw(ArgumentError("geometry_m must be three sensor rows x three coordinates"))
    maximum(abs.(positions[:,3] .- positions[3,3])) < 1e-12 || throw(ArgumentError("projection requires a planar D3 array"))
    baseline = positions[1:2,1:2] .- permutedims(positions[3,1:2])
    abs(det(baseline)) > 1e-12 || throw(ArgumentError("degenerate sensor geometry"))
    offsets = Float64.(calibration["relative_delay_offsets_s"])
    c = Float64(calibration["sound_speed_m_s"])
    c > 0 && isfinite(c) || throw(ArgumentError("sound speed must be positive finite m/s"))
    uv = -(baseline \ (dt[1:2] .- offsets[1:2])) .* c
    disk = sum(abs2, uv)
    disk <= 1.0 || return _pf_invalid("impossible_tdoa")
    return _pf_project_direction(calibration, [uv[1], uv[2], sqrt(1.0 - disk)])
end

function project_tdoas(calibration::AbstractDict, delays::AbstractMatrix; fs=nothing, units="samples")
    size(delays,2) in (2,3) || throw(ArgumentError("delays must be N×2 or N×3"))
    fs isa AbstractVector && length(fs) != size(delays,1) && throw(DimensionMismatch("one fs per delay row required"))
    return [project_tdoas(calibration, collect(row); fs=fs isa AbstractVector ? fs[i] : fs, units=units)
        for (i,row) in enumerate(eachrow(delays))]
end

"""
    project_angles(calibration, angles; units="radians")

Project Julia azimuth/elevation pairs (length 2 or N×2) using
u=[sin(az)*cos(el), sin(el), cos(az)*cos(el)]. Both angles are radians unless
units="degrees". +x right, +y up, +z forward; the +z array branch is required.
The output/invalid conventions match project_tdoas. No delay correction is
applied to angles: angles must already represent the corrected array ray.
"""
function project_angles(calibration::AbstractDict, angles::AbstractVector; units="radians")
    length(angles) == 2 || throw(ArgumentError("angles must contain Julia azimuth and elevation"))
    units in ("radians", "degrees") || throw(ArgumentError("units must be radians or degrees"))
    all(_pf_finite, angles) || return _pf_invalid("nonfinite_angle")
    az, el = Float64.(angles) .* (units == "degrees" ? pi / 180 : 1.0)
    abs(el) <= pi/2 || return _pf_invalid("invalid_elevation")
    u = [sin(az)*cos(el), sin(el), cos(az)*cos(el)]
    u[3] >= 0 || return _pf_invalid("back_array_branch"; direction_array=u)
    return _pf_project_direction(calibration, u)
end

function project_angles(calibration::AbstractDict, angles::AbstractMatrix; units="radians")
    size(angles,2) == 2 || throw(ArgumentError("angles must be N×2"))
    return [project_angles(calibration, collect(row); units=units) for row in eachrow(angles)]
end

# Internal normalized parameter box. Rotation is Euler Rz*Ry*Rx, radians;
# sound speed 1500+100q, delay corrections 100μs*q. Focal lengths span 20..150°
# FOV independently. Optional principal point ±10% dimension; k1 ±0.3, k2=0.
function _pf_rotation(angles)
    sx,cx = sincos(angles[1]); sy,cy = sincos(angles[2]); sz,cz = sincos(angles[3])
    return [cz*cy cz*sy*sx-sz*cx cz*sy*cx+sz*sx;
            sz*cy sz*sy*sx+cz*cx sz*sy*cx-cz*sx;
            -sy cy*sx cy*cx]
end

function _pf_decode(q, width, height, model)
    fmin = [width,height] ./ (2*tan(deg2rad(75.0)))
    fmax = [width,height] ./ (2*tan(deg2rad(10.0)))
    focal = exp.(log.(fmin) .+ (q[5:6] .+ 1) ./ 2 .* log.(fmax ./ fmin))
    cx,cy = (width-1)/2, (height-1)/2
    model != "centered_pinhole" && (cx += 0.1*width*q[9]; cy += 0.1*height*q[10])
    return (c=1500.0+100.0*q[1], R=_pf_rotation(pi .* q[2:4]), fx=focal[1], fy=focal[2],
        cx=cx,cy=cy,b=1e-4 .* q[7:8], k1=model == "radial_pinhole" ? 0.3*q[11] : 0.0)
end

_pf_nparams(model) = model == "centered_pinhole" ? 8 : model == "principal_pinhole" ? 10 : model == "radial_pinhole" ? 11 : throw(ArgumentError("unknown physical model $model"))

function _pf_initial(width, height, model)
    q = zeros(_pf_nparams(model))
    # Independent 80° horizontal / 55° vertical starting pinhole.
    for (index,dimension,fov) in ((5,width,80.0),(6,height,55.0))
        flo = dimension / (2*tan(deg2rad(75.0)))
        fhi = dimension / (2*tan(deg2rad(10.0)))
        f = dimension / (2*tan(deg2rad(fov/2)))
        q[index] = 2*log(f/flo)/log(fhi/flo)-1
    end
    return q
end

function _pf_pixels(p, tx, ty, inverse_baseline)
    ux = -p.c*(inverse_baseline[1,1]*(tx-p.b[1])+inverse_baseline[1,2]*(ty-p.b[2]))
    uy = -p.c*(inverse_baseline[2,1]*(tx-p.b[1])+inverse_baseline[2,2]*(ty-p.b[2]))
    disk = ux*ux+uy*uy
    # Optimizer gets an explicit physical constraint penalty; this path never
    # generates a ray/pixel by clamping an impossible observation.
    disk <= 1.0 || return (NaN,NaN,1.0+disk-1.0)
    uz = sqrt(1.0-disk)
    X=p.R[1,1]*ux+p.R[1,2]*uy+p.R[1,3]*uz
    Y=p.R[2,1]*ux+p.R[2,2]*uy+p.R[2,3]*uz
    Z=p.R[3,1]*ux+p.R[3,2]*uy+p.R[3,3]*uz
    Z > 1e-6 || return (NaN,NaN,2.0-Z)
    xn,yn=X/Z,Y/Z
    radial=1+p.k1*(xn*xn+yn*yn)
    derivative=1+3p.k1*(xn*xn+yn*yn)
    radial > 0 && derivative > 0 || return (NaN,NaN,1.0+abs(min(radial,derivative)))
    return (p.cx+p.fx*xn*radial,p.cy-p.fy*yn*radial,0.0)
end

"""True only if no allowed c/delay correction can put this observation on the
physical disk. Minimize ||B⁻¹(tau-b)||² over b∈[-100,100]μs² exactly by checking
the four rectangle edges (and its interior zero); c=1400 minimizes disk radius.
This fixed input-only screen uses no annotations or fitted parameters. Rows
remain in projection/coverage reports; they contribute no fitting influence.
"""
function _pf_definitely_impossible(tx,ty,inverse_baseline)
    isfinite(tx) && isfinite(ty) || return true
    lowx,highx=tx-1e-4,tx+1e-4;lowy,highy=ty-1e-4,ty+1e-4
    lowx <= 0 <= highx && lowy <= 0 <= highy && return false
    a,b=inverse_baseline[1,1],inverse_baseline[1,2]
    c,d=inverse_baseline[2,1],inverse_baseline[2,2]
    g11=a*a+c*c;g12=a*b+c*d;g22=b*b+d*d
    minimum_norm=Inf
    for x in (lowx,highx)
        y=clamp(-g12*x/g22,lowy,highy)
        minimum_norm=min(minimum_norm,(a*x+b*y)^2+(c*x+d*y)^2)
    end
    for y in (lowy,highy)
        x=clamp(-g12*y/g11,lowx,highx)
        minimum_norm=min(minimum_norm,(a*x+b*y)^2+(c*x+d*y)^2)
    end
    # Conservative roundoff guard affects screening only, never ray projection.
    return 1400.0^2*minimum_norm > 1.0+64eps(Float64)
end

function _pf_screen_summary(data)
    screened=[i for i in eachindex(data.x) if _pf_definitely_impossible(data.tau[i,1],data.tau[i,2],data.inverse_baseline)]
    return Dict("rule"=>"input-only exact minimum physical disk radius over c1400..1600m/s and both relative offsets +/-100us; no annotation or fitted-parameter dependence",
        "n_eligible"=>length(data.x),"n_screened_from_objective"=>length(screened),
        "screened_event_rows"=>data.event_row[screened],"retained_in_projection_coverage"=>true)
end

function _pf_loss(q, data, width, height, model, delta)
    excess = sum(abs2, max.(abs.(q) .- 1.0, 0.0))
    excess > 0 && return 1e12*(1+excess)
    p = _pf_decode(q,width,height,model)
    total=0.0; weight=0.0
    # A bounded, smooth redescending transform of pseudo-Huber prevents both
    # impossible rays and arbitrarily large valid pixels from dominating. The
    # fixed invalid-ray contribution is the same cap at every parameter value.
    saturation_residual=max(8delta,0.05hypot(width,height))
    cap=2delta^2*(hypot(1.0,saturation_residual/delta)-1.0)
    for i in eachindex(data.x)
        data.weight[i] > 0 || continue
        _pf_definitely_impossible(data.tau[i,1],data.tau[i,2],data.inverse_baseline) && continue
        x,y,invalid = _pf_pixels(p,data.tau[i,1],data.tau[i,2],data.inverse_baseline)
        contribution=cap
        if invalid == 0
            r=hypot(x-data.x[i],y-data.y[i])
            rho=2delta^2*(hypot(1.0,r/delta)-1.0)
            contribution=cap*(rho/(cap+rho))
        end
        total += data.weight[i]*contribution
        weight += data.weight[i]
    end
    return total / max(weight, eps())
end

function _pf_fit(data,width,height,model; nstarts=6,maxiters=2500,delta=12.0,initial=nothing,fixed=Dict{Int,Float64}())
    isempty(data.x) && throw(ArgumentError("no fit observations"))
    any(i->data.weight[i] > 0 && !_pf_definitely_impossible(data.tau[i,1],data.tau[i,2],data.inverse_baseline),eachindex(data.x)) || throw(ArgumentError("all observations are definitely impossible under allowed sound-speed/delay bounds"))
    nstarts >= 1 || throw(ArgumentError("nstarts must be positive"))
    base = _pf_initial(width,height,model)
    initial !== nothing && (base[1:min(length(base),length(initial))] .= initial[1:min(length(base),length(initial))])
    free = setdiff(collect(eachindex(base)),collect(keys(fixed)))
    for (i,value) in fixed; base[i]=value; end
    rng=MersenneTwister(20261030)
    starts=Vector{Vector{Float64}}([copy(base)])
    for j in 2:nstarts
        trial=copy(base)
        trial[1]=[-0.65,0.65,0.0][mod1(j-1,3)]
        trial[2:4] .+= 0.08 .* randn(rng,3)
        trial[5:6] .+= 0.10 .* randn(rng,2)
        trial[7:8] .+= 0.08 .* randn(rng,2)
        # Cover image/array roll alternatives as well as the local mounting.
        j >= 4 && (trial[4]=[0.5,-0.5,0.95][mod1(j-3,3)])
        trial .= clamp.(trial,-0.95,0.95)
        for (i,value) in fixed; trial[i]=value; end
        push!(starts,trial)
    end
    evidence=Any[]; best_q=copy(base); best_value=Inf
    for (number,start) in enumerate(starts)
        function objective(z)
            q=copy(base); q[free]=z
            return _pf_loss(q,data,width,height,model,delta)
        end
        result=Optim.optimize(objective,start[free],Optim.NelderMead(),
            Optim.Options(iterations=maxiters, x_abstol=1e-9, f_abstol=1e-8, show_trace=false))
        q=copy(base); q[free]=Optim.minimizer(result)
        value=objective(Optim.minimizer(result))
        push!(evidence,Dict("start"=>number,"initial_q"=>start,"objective"=>value,
            "iterations"=>Optim.iterations(result),"converged"=>Optim.converged(result),"q"=>q))
        if isfinite(value) && value < best_value
            best_value=value; best_q=q
        end
    end
    # A local finite-difference BFGS refinement improves the long shallow c/FOV
    # valley after deterministic simplex starts; native Optim only.
    function refine(z)
        q=copy(base); q[free]=z
        return _pf_loss(q,data,width,height,model,delta)
    end
    refinement=Optim.optimize(refine,best_q[free],Optim.BFGS(),
        Optim.Options(iterations=maxiters, g_abstol=1e-6, show_trace=false))
    refined=copy(base); refined[free]=Optim.minimizer(refinement)
    value=refine(Optim.minimizer(refinement))
    push!(evidence,Dict("start"=>"BFGS_refinement","objective"=>value,
        "iterations"=>Optim.iterations(refinement),"converged"=>Optim.converged(refinement),"q"=>refined))
    value < best_value && (best_value=value; best_q=refined)
    return (q=best_q,objective=best_value,starts=evidence)
end

function _pf_calibration(q,width,height,model; device="",method="minmax2")
    p=_pf_decode(q,width,height,model)
    return Dict{String,Any}("schema_version"=>"clicker-physical-calibration/1.0",
        "device"=>device,"width"=>width,"height"=>height,"geometry_m"=>_pf_rows(clicker_geometry()),
        "physical_channels"=>[3,4,5],"reference_channel"=>3,"reference_physical_channel"=>5,
        "sound_speed_m_s"=>p.c,"relative_delay_offsets_s"=>[p.b[1],p.b[2],0.0],
        "rotation_array_to_camera"=>_pf_rows(p.R),"rotation_euler_xyz_rad"=>pi .* q[2:4],
        "camera"=>Dict("fx"=>p.fx,"fy"=>p.fy,"cx"=>p.cx,"cy"=>p.cy,"k1"=>p.k1,"k2"=>0.0),
        "fov_degrees"=>[rad2deg(2atan(width/(2p.fx))),rad2deg(2atan(height/(2p.fy)))],
        "method"=>method,"model"=>model,"numeric_precision"=>"Float64",
        "software"=>Dict("Julia"=>string(VERSION),"Optim"=>string(pkgversion(Optim)),
            "CSV"=>string(pkgversion(CSV)),"DataFrames"=>string(pkgversion(DataFrames)),
            "optimizer"=>"native Julia Optim.NelderMead deterministic starts, finite-difference Optim.BFGS refinement"),
        "robust_objective"=>"bounded smooth transform cap*rho/(cap+rho) of pseudo-Huber rho, cap set at residual max(8delta,0.05 image diagonal); invalid rays contribute fixed cap; definitely-impossible observations screened only from objective and retained in coverage",
        "limits"=>Dict("sound_speed_m_s"=>[1400.0,1600.0],"fov_degrees"=>[20.0,150.0],
            "relative_delay_offsets_s"=>[-1e-4,1e-4],"principal_offset_fraction"=>[-0.1,0.1],"k1"=>[-0.3,0.3]),
        "source_model_spec"=>Dict("geometry_layout"=>"rows selected physical3,4,5; columns array x,y,z metres",
            "tdoa_equation"=>"dt_i=-(r_i-r_ref) dot u / sound_speed + relative_delay_offset_i",
            "angle_equation"=>"u=[sin(az)*cos(el), sin(el), cos(az)*cos(el)]",
            "rotation"=>"SO(3), Rz*Ry*Rx; array to camera; both x right,y up,z forward",
            "projection"=>"x=cx+fx*X/Z*(1+k1*r2+k2*r2^2); y=cy-fy*Y/Z*(1+k1*r2+k2*r2^2)",
            "pixel_origin"=>"zero; centered principal point ((width-1)/2,(height-1)/2)",
            "front_branch"=>"array +z only; camera Z>0; impossible TDOAs invalid without clamping",
            "caveat"=>"Far-field planar direction-only model. Array cannot resolve front/back or range; near-field curvature, refraction and propagation paths are unmodelled."))
end

function _pf_subset(data,indices)
    return (tau=data.tau[indices,:],x=data.x[indices],y=data.y[indices],weight=data.weight[indices],
        block=data.block[indices],event_row=data.event_row[indices],frame=data.frame[indices],
        inverse_baseline=data.inverse_baseline)
end

function _pf_metrics(calibration,data)
    errors=Float64[]; ninside=0; statuses=Dict{String,Int}()
    for i in eachindex(data.x)
        projected=project_tdoas(calibration,collect(data.tau[i,:]);units="seconds")
        statuses[projected.status]=get(statuses,projected.status,0)+1
        if projected.valid
            push!(errors,hypot(projected.x-data.x[i],projected.y-data.y[i]))
            ninside += projected.in_frame
        end
    end
    n=length(data.x); valid=length(errors)
    return Dict{String,Any}("n_eligible"=>n,"n_unique_frames"=>length(unique(data.frame)),
        "n_projected"=>valid,"coverage"=>n == 0 ? nothing : valid/n,
        "in_image_coverage"=>n == 0 ? nothing : ninside/n,"status_counts"=>statuses,
        "rmse_px"=>isempty(errors) ? nothing : sqrt(sum(abs2,errors)/valid),
        "median_px"=>isempty(errors) ? nothing : median(errors),
        "p90_px"=>isempty(errors) ? nothing : quantile(errors,0.9),
        "error_definition"=>"Euclidean pixels on valid projections; coverage denominator includes every eligible event; outside-image front rays count as projected")
end

function _pf_validation_score(calibration,data,width,height)
    total=0.0; weights=0.0
    for i in eachindex(data.x)
        p=project_tdoas(calibration,collect(data.tau[i,:]);units="seconds")
        err = p.valid ? hypot(p.x-data.x[i],p.y-data.y[i]) : 3hypot(width,height)
        total += data.weight[i]*err^2; weights += data.weight[i]
    end
    return sqrt(total/max(weights,eps()))
end

"""Uncalibrated dimension-corrected D3 linear-angle baseline, not a fitted model.
Uses c=1540 m/s, zero delay corrections, identity mounting, original
angle2pxd equation with FOV=(80,54.46086642417594) degrees. Image dimensions
come from the actual frame. The historical 2740 width is never substituted.
Impossible TDOAs remain invalid, avoiding legacy unit-disk clamping.
"""
function _pf_d3_baseline(delays,width,height;fs=nothing,units="seconds")
    length(delays) in (2,3) || throw(ArgumentError("D3 baseline needs 2 or 3 delays"))
    all(_pf_finite,delays) || return _pf_invalid("nonfinite_tdoa")
    units in ("samples","seconds") || throw(ArgumentError("invalid baseline delay units"))
    units == "samples" && (!_pf_finite(fs) || fs <= 0) && return _pf_invalid("invalid_fs")
    tau=Float64.(delays) ./ (units == "samples" ? Float64(fs) : 1.0)
    length(tau) == 3 && abs(tau[3]) > 1e-12 && return _pf_invalid("invalid_reference")
    geometry=clicker_geometry();baseline=geometry[1:2,1:2].-permutedims(geometry[3,1:2])
    uv=-1540.0 .* (baseline \ tau[1:2]);disk=sum(abs2,uv)
    disk <= 1 || return _pf_invalid("impossible_tdoa")
    u=[uv[1],uv[2],sqrt(1-disk)]
    az=atan(u[1],u[3]);el=asin(u[2])
    x=width/2+width*rad2deg(az)/80.0
    y=height/2-height*rad2deg(el)/54.46086642417594
    inside=0 <= x < width && 0 <= y < height
    return (x=x,y=y,valid=true,status=inside ? "ok" : "outside_image",
        direction_array=u,direction_camera=copy(u),in_frame=inside)
end

function _pf_error_summary(errors)
    return Dict("n"=>length(errors),"rmse_px"=>isempty(errors) ? nothing : sqrt(mean(abs2,errors)),
        "median_px"=>isempty(errors) ? nothing : median(errors),
        "p90_px"=>isempty(errors) ? nothing : quantile(errors,0.9))
end

function _pf_paired_baseline(data,calibration,width,height)
    baseline_errors=Float64[];calibrated_errors=Float64[];common_baseline=Float64[];common_calibrated=Float64[]
    for i in eachindex(data.x)
        raw=_pf_d3_baseline(collect(data.tau[i,:]),width,height)
        fitted=project_tdoas(calibration,collect(data.tau[i,:]);units="seconds")
        raw.valid && push!(baseline_errors,hypot(raw.x-data.x[i],raw.y-data.y[i]))
        fitted.valid && push!(calibrated_errors,hypot(fitted.x-data.x[i],fitted.y-data.y[i]))
        if raw.valid && fitted.valid
            push!(common_baseline,last(baseline_errors));push!(common_calibrated,last(calibrated_errors))
        end
    end
    n=length(data.x)
    return Dict("n_eligible"=>n,"baseline_coverage"=>n == 0 ? nothing : length(baseline_errors)/n,
        "calibrated_coverage"=>n == 0 ? nothing : length(calibrated_errors)/n,
        "common_coverage"=>n == 0 ? nothing : length(common_baseline)/n,
        "baseline_all_valid"=>_pf_error_summary(baseline_errors),"calibrated_all_valid"=>_pf_error_summary(calibrated_errors),
        "paired_common_event_baseline"=>_pf_error_summary(common_baseline),
        "paired_common_event_calibrated"=>_pf_error_summary(common_calibrated),
        "paired_median_error_change_px"=>isempty(common_baseline) ? nothing : median(common_calibrated-common_baseline))
end

function _pf_baseline_spec(width,height)
    return Dict("label"=>"dimension-corrected D3 baseline","fitted"=>false,
        "sound_speed_m_s"=>1540.0,"relative_delay_offsets_s"=>[0.0,0.0,0.0],
        "rotation_array_to_camera"=>_pf_rows(Matrix{Float64}(I,3,3)),
        "fov_degrees"=>[80.0,54.46086642417594],"actual_width"=>width,"actual_height"=>height,
        "pixel_equation"=>"x=width/2+width*az_degrees/80; y=height/2-height*el_degrees/54.46086642417594",
        "historical_dimension_bug"=>"Original config declared width2740; calibration uses actual decoded image dimensions (real device2704x1520). This baseline corrects dimensions, preserving the original linear-angle equation.",
        "invalid_rule"=>"Impossible TDOAs invalid, never clamped; +z front array branch")
end

function _pf_select_model(data,width,height; allow_extended=true,nstarts=6,maxiters=2500,delta=12.0)
    blocks=sort(unique(data.block))
    evidence=Dict{String,Any}("uses_heldout"=>false,"rule"=>"predeclared nested blocked validation; extension requires >=5% mean penalized RMSE gain and no fold degradation >2%; invalid rays penalized at 3 image diagonals",
        "candidates"=>Any[],"selected"=>"centered_pinhole")
    # No extra optical freedom without enough independent training blocks.
    if !allow_extended || length(blocks) < 4 || length(unique(data.frame)) < 40
        evidence["reason"]="Centered pinhole retained: extensions disabled or fewer than 4 training blocks / 40 annotated frames"
        return "centered_pinhole",evidence
    end
    models=["centered_pinhole","principal_pinhole","radial_pinhole"]
    scores=Dict{String,Vector{Float64}}()
    for model in models
        folds=Any[]; scores[model]=Float64[]
        for fold in 1:3
            validation_blocks=blocks[fold:3:end]
            vi=findall(in(validation_blocks),data.block); ti=findall(!in(validation_blocks),data.block)
            validation_frames=Set(data.frame[vi])
            # Guard exact frame-label leakage even for irregular/duplicate input
            # timestamps that defeat the usual two-second temporal purge.
            ti=filter(i->!(data.frame[i] in validation_frames),ti)
            training=_pf_subset(data,ti); validation=_pf_subset(data,vi)
            fit=_pf_fit(training,width,height,model;nstarts=nstarts,maxiters=maxiters,delta=delta)
            cal=_pf_calibration(fit.q,width,height,model)
            score=_pf_validation_score(cal,validation,width,height)
            push!(scores[model],score)
            push!(folds,Dict("fold"=>fold,"training_blocks"=>sort(unique(training.block)),
                "validation_blocks"=>validation_blocks,"penalized_rmse_px"=>score,"metrics"=>_pf_metrics(cal,validation),
                "q"=>fit.q,"starts"=>fit.starts))
        end
        push!(evidence["candidates"],Dict("model"=>model,"folds"=>folds,"mean_penalized_rmse_px"=>mean(scores[model])))
    end
    selected="centered_pinhole"
    for candidate in models[2:end]
        if mean(scores[candidate]) < 0.95mean(scores[selected]) && all(scores[candidate] .<= 1.02 .* scores[selected])
            selected=candidate
        end
    end
    evidence["selected"]=selected
    return selected,evidence
end

function _pf_residuals(q,data,width,height,model)
    p=_pf_decode(q,width,height,model); residual=Float64[]; rows=Int[]
    for i in eachindex(data.x)
        x,y,invalid=_pf_pixels(p,data.tau[i,1],data.tau[i,2],data.inverse_baseline)
        if invalid == 0
            push!(residual,x-data.x[i],y-data.y[i]); push!(rows,i)
        end
    end
    return residual,rows
end

function _pf_diagnostics(q,data,width,height,model; profiles=true,maxiters=1200,delta=12.0)
    residual,validrows=_pf_residuals(q,data,width,height,model)
    localdata=_pf_subset(data,validrows); npar=length(q)
    length(residual) >= npar || return Dict("status"=>"too_few_valid_rays_for_jacobian",
        "weak_identifiability"=>true,"parameter_bound_hits"=>findall(abs.(q) .>= 0.98),
        "uncertainty"=>"Not estimable: fewer valid pixel residuals than free parameters")
    J=zeros(length(residual),npar)
    for j in 1:npar
        step=1e-5; plus=copy(q); minus=copy(q); plus[j]+=step; minus[j]-=step
        rp,_=_pf_residuals(plus,localdata,width,height,model)
        rm,_=_pf_residuals(minus,localdata,width,height,model)
        if length(rp) != length(residual) || length(rm) != length(residual)
            return Dict("status"=>"boundary_invalidates_jacobian","weak_identifiability"=>true,
                "parameter_bound_hits"=>findall(abs.(q) .>= 0.98),"uncertainty"=>"Not estimable at physical branch boundary")
        end
        J[:,j]=(rp-rm)/(2step)
    end
    decomposition=svd(J;full=false); s=decomposition.S
    rank=isempty(s) ? 0 : count(s .> max(maximum(s)*1e-6,1e-8))
    condition=length(s) < npar || isempty(s) || minimum(s) <= 1e-12 ? Inf : maximum(s)/minimum(s)
    # Match the fitted frame/uncertainty weighting and pseudo-Huber influence
    # in a robust Gauss-Newton sandwich, rather than treating repeated events
    # or annotated outliers as independent equal-precision measurements.
    robust_weights=Float64[]
    saturation_residual=max(8delta,0.05hypot(width,height))
    cap=2delta^2*(hypot(1.0,saturation_residual/delta)-1.0)
    for i in eachindex(localdata.x)
        r=hypot(residual[2i-1],residual[2i]);rho=2delta^2*(hypot(1.0,r/delta)-1.0)
        w=localdata.weight[i]/hypot(1.0,r/delta)*(cap/(cap+rho))^2
        push!(robust_weights,w,w)
    end
    weighted_J=J .* sqrt.(robust_weights)
    bread=pinv(weighted_J'weighted_J;rtol=1e-12)
    # Cluster sandwich by 30s block: repeated candidate events/frame annotations
    # are not independent. This is a local linear approximation, not a physical
    # confidence guarantee or a substitute for missing directional coverage.
    meat=zeros(npar,npar)
    for block in unique(localdata.block)
        obs=findall(==(block),localdata.block)
        rows=reduce(vcat,([2i-1,2i] for i in obs);init=Int[])
        score=J[rows,:]'*(robust_weights[rows].*residual[rows])
        meat += score*score'
    end
    nblocks=length(unique(localdata.block))
    covariance=nblocks > 1 ? bread*meat*bread*(nblocks/(nblocks-1)) : fill(NaN,npar,npar)
    se=sqrt.(max.(diag(covariance),0.0))
    if rank < npar
        null_influence=vec(sum(abs2,decomposition.V[:,rank+1:end];dims=2))
        se[null_influence .> 1e-6] .= NaN
    end
    nblocks < 4 && (se .= NaN)
    names=["c_normalized","rx_over_pi","ry_over_pi","rz_over_pi","log_fx_normalized","log_fy_normalized","b1_over_100us","b2_over_100us"]
    model != "centered_pinhole" && append!(names,["cx_over_0.1width","cy_over_0.1height"])
    model == "radial_pinhole" && push!(names,"k1_over_0.3")
    weak=rank < npar || condition > 1e5 || nblocks < 4
    profile_results=Any[]
    if profiles
        # Nuisance parameters refit on training data ONLY. Three speed values
        # expose c/FOV/delay tradeoffs without making heldout tuning claims.
        for c in (1400.0,1500.0,1600.0)
            fit=_pf_fit(data,width,height,model;nstarts=2,maxiters=maxiters,delta=delta,initial=q,fixed=Dict(1=>(c-1500)/100))
            cal=_pf_calibration(fit.q,width,height,model)
            push!(profile_results,Dict("fixed_sound_speed_m_s"=>c,"training_objective"=>fit.objective,
                "training_metrics"=>_pf_metrics(cal,data),"fov_degrees"=>cal["fov_degrees"],
                "relative_delay_offsets_s"=>cal["relative_delay_offsets_s"],"q"=>fit.q))
        end
        best=_pf_loss(q,data,width,height,model,delta)
        acceptable=[r["fixed_sound_speed_m_s"] for r in profile_results if r["training_objective"] <= best*1.05+0.1]
        weak |= length(acceptable) >= 2
    end
    return Dict{String,Any}("status"=>"local_linear_diagnostics","parameter_names_normalized"=>names,
        "parameter_q"=>q,"jacobian_singular_values"=>s,"jacobian_rank"=>rank,
        "jacobian_condition_number"=>condition,"weak_identifiability"=>weak,
        "parameter_bound_hits"=>[names[i] for i in findall(abs.(q) .>= 0.98)],
        "weakest_normalized_direction"=>isempty(s) ? nothing : Dict(names[i]=>decomposition.V[i,end] for i in 1:npar),
        "cluster_blocks"=>nblocks,"cluster_sandwich_se_normalized"=>se,
        "sound_speed_local_se_m_s"=>100se[1],"delay_local_se_s"=>1e-4 .* se[7:8],
        "rotation_local_se_rad"=>pi .* se[2:4],
        "training_sound_speed_profiles"=>profile_results,
        "uncertainty"=>"30-second-block cluster sandwich with frame/uncertainty and pseudo-Huber weights; Gauss-Newton linearization, conditional on model and fixed front branch, approximate, unreliable with <4 blocks or rank loss. Null directions are not assigned finite confidence intervals. c, FOV, rotation and relative delays can trade off; speed profiles refit nuisances on training only.")
end

# This opt-in policy describes provisional observations, never revised annotations.
function _pf_policy_inputs(output_dir)
    paths=Dict("events_sha256"=>"events.csv","annotations_sha256"=>"annotations.csv",
        "frames_sha256"=>"frames.csv","annotations_complete_sha256"=>"annotations.complete.json")
    hashes=Dict(key=>sha256_file(joinpath(output_dir,name)) for (key,name) in paths)
    completion=read_json(joinpath(output_dir,"annotations.complete.json"))
    for key in ("frames_sha256","annotations_sha256")
        get(completion,key,nothing)==hashes[key] || throw(ArgumentError("annotation completion hash mismatch: $key"))
    end
    return hashes
end

function _pf_validate_observation_policy(policy,hashes)
    get(policy,"schema_version",nothing)=="clicker-fit-observation-policy/1.0" || throw(ArgumentError("unknown fit observation policy schema"))
    get(policy,"policy",nothing)=="terra_estimated_underwater_coordinates" || throw(ArgumentError("unknown fit observation policy"))
    get(policy,"human_confirmed",nothing) === false || throw(ArgumentError("observation policy must not claim human confirmation"))
    get(policy,"provisional_visual_estimates",nothing) === true || throw(ArgumentError("observation policy must declare provisional visual estimates"))
    get(policy,"allowed_visibility",nothing)==["visible","uncertain"] || throw(ArgumentError("observation policy visibility rule changed"))
    get(policy,"required_medium",nothing)=="underwater" || throw(ArgumentError("observation policy medium rule changed"))
    get(policy,"require_positive_finite_original_radius",nothing) === true || throw(ArgumentError("observation policy radius rule changed"))
    get(policy,"policy_source_sha256",nothing)==sha256_file(@__FILE__) || throw(ArgumentError("observation policy source changed"))
    isempty(strip(_pf_text(get(policy,"reason",nothing)))) && throw(ArgumentError("observation policy needs an explicit reason"))
    for key in ("events_sha256","annotations_sha256","frames_sha256","annotations_complete_sha256")
        get(get(policy,"hashes",Dict()),key,nothing)==hashes[key] || throw(ArgumentError("observation policy input changed: $key"))
    end
    return policy
end

"""
    predeclare_observation_policy(output_dir;
        policy=:terra_estimated_underwater_coordinates, reason) -> Dict

Opt in before any first fit to original visible/uncertain underwater coordinate
estimates with finite positive original radii. Labels and coordinates stay intact;
this is a provisional visual-observation policy with human_confirmed=false.
The sealed annotation/table hashes and this source are bound. Existing policies
resume only when valid; existing fit products prohibit a new declaration.
"""
function predeclare_observation_policy(output_dir;policy=:terra_estimated_underwater_coordinates,reason)
    string(policy)=="terra_estimated_underwater_coordinates" || throw(ArgumentError("unsupported observation policy"))
    isempty(strip(string(reason))) && throw(ArgumentError("explicit observation policy reason required"))
    hashes=_pf_policy_inputs(output_dir)
    path=joinpath(output_dir,"fit_observation_policy.json")
    if isfile(path)
        existing=_pf_validate_observation_policy(read_json(path),hashes)
        existing["reason"]==string(reason) || throw(ArgumentError("existing observation policy reason differs"))
        return existing
    end
    any(isfile(joinpath(output_dir,name)) for name in ("calibration.json","fit_report.json","predictions.csv","heldout_evaluation.json")) &&
        throw(ArgumentError("declare observation policy in a fresh fit directory before fitting"))
    declaration=Dict("schema_version"=>"clicker-fit-observation-policy/1.0","policy"=>string(policy),
        "human_confirmed"=>false,"provisional_visual_estimates"=>true,
        "allowed_visibility"=>["visible","uncertain"],"required_medium"=>"underwater",
        "require_positive_finite_original_radius"=>true,"reason"=>string(reason),"hashes"=>hashes,
        "policy_source_sha256"=>sha256_file(@__FILE__),"declared_at_utc"=>string(Dates.now(Dates.UTC))*"Z",
        "limitations"=>"Original Terra housing-center estimates and radii are provisional visual estimates, not human-confirmed positions or calibrated confidence intervals. No labels/coordinates are changed; sync, geometry, purge and shared-frame exclusions still apply.")
    _pf_validate_observation_policy(declaration,hashes)
    write_json(path,declaration)
    return declaration
end

function _pf_observation_evidence(input)
    return Dict("policy"=>input.observation_policy,
        "policy_sha256"=>get(input.hashes,"fit_observation_policy_sha256",nothing),
        "interpretation"=>input.observation_policy===nothing ? "strict original visible observations" : "provisional original visual estimates; human_confirmed=false",
        "events"=>[Dict("event_id"=>_pf_text(input.events.event_id[i]),"frame_id"=>_pf_text(input.events.frame_id[i]),
            "original_visibility"=>input.original_visibility[i],"original_uncertainty_px"=>input.original_uncertainty[i],
            "observation_admission"=>input.observation_admission[i],"eligibility_status"=>input.reason[i],"split"=>input.split[i]) for i in eachindex(input.reason)])
end

function _pf_load(output_dir,method)
    event_path=joinpath(output_dir,"events.csv"); annotation_path=joinpath(output_dir,"annotations.csv")
    frames_path=joinpath(output_dir,"frames.csv"); complete_path=joinpath(output_dir,"annotations.complete.json")
    for path in (event_path,annotation_path,frames_path,complete_path)
        isfile(path) || throw(ArgumentError("missing required input $path"))
    end
    completion=read_json(complete_path)
    hashes=Dict("events_sha256"=>sha256_file(event_path),"annotations_sha256"=>sha256_file(annotation_path),
        "frames_sha256"=>sha256_file(frames_path),"annotations_complete_sha256"=>sha256_file(complete_path),
        "physical_fit_source_sha256"=>sha256_file(@__FILE__))
    support_path=joinpath(@__DIR__,"support.jl")
    isfile(support_path) && (hashes["support_source_sha256"]=sha256_file(support_path))
    for key in ("frames_sha256","annotations_sha256")
        get(completion,key,nothing) == hashes[key] || throw(ArgumentError("annotation completion hash mismatch: $key"))
    end
    policy_path=joinpath(output_dir,"fit_observation_policy.json")
    observation_policy=isfile(policy_path) ? _pf_validate_observation_policy(read_json(policy_path),hashes) : nothing
    observation_policy===nothing || (hashes["fit_observation_policy_sha256"]=sha256_file(policy_path))
    events=CSV.read(event_path,DataFrame;stringtype=String)
    annotations=CSV.read(annotation_path,DataFrame;stringtype=String)
    frames=CSV.read(frames_path,DataFrame;stringtype=String)
    required=["event_id","device","frame_id","width","height","fs","recording_time_s","sync_status",method*"_d1",method*"_d2",method*"_d3"]
    all(in(names(events)),required) || throw(ArgumentError("events.csv missing required fitting columns"))
    all(in(names(annotations)),["frame_id","x","y","visibility","uncertainty_px","medium"]) || throw(ArgumentError("annotations.csv missing required columns (including medium)"))
    all(value->_pf_text(value) in ("underwater","air","mixed","unknown"),annotations.medium) || throw(ArgumentError("invalid annotation medium; expected underwater/air/mixed/unknown"))
    length(unique(annotations.frame_id)) == nrow(annotations) || throw(ArgumentError("duplicate annotation frame_id"))
    Set(_pf_text.(frames.frame_id)) == Set(_pf_text.(annotations.frame_id)) || throw(ArgumentError("annotation rows must cover exactly frames.csv"))
    annotation_map=Dict(_pf_text(row.frame_id)=>row for row in eachrow(annotations))
    widths=Int[]; heights=Int[]; rows=Int[]; blocks=Int[]; frame_ids=String[]; xs=Float64[]; ys=Float64[]; uncertainty=Float64[]
    tau=Vector{Vector{Float64}}(); reason=fill("",nrow(events)); split=fill("ineligible",nrow(events))
    annotated_x=fill(NaN,nrow(events)); annotated_y=fill(NaN,nrow(events))
    original_visibility=fill("",nrow(events)); original_uncertainty=fill(NaN,nrow(events))
    observation_admission=fill("no_annotation",nrow(events))
    for (i,row) in enumerate(eachrow(events))
        frame=_pf_text(row.frame_id)
        if !haskey(annotation_map,frame)
            reason[i]="no_annotation"; continue
        end
        ann=annotation_map[frame]
        annotated_x[i]=_pf_number(ann.x); annotated_y[i]=_pf_number(ann.y)
        original_visibility[i]=_pf_text(ann.visibility);original_uncertainty[i]=_pf_number(ann.uncertainty_px)
        allowed=original_visibility[i] == "visible" || (observation_policy!==nothing && original_visibility[i] == "uncertain")
        valid_radius=_pf_finite(ann.uncertainty_px) && ann.uncertainty_px > 0
        valid_bounds=_pf_finite(row.width) && _pf_finite(row.height) && row.width >= 2 && row.height >= 2 &&
            _pf_finite(ann.x) && _pf_finite(ann.y) && 0 <= ann.x < row.width && 0 <= ann.y < row.height
        # Annotation admission describes only the declared visual rule. Acoustic
        # synchronization, timing purges and physics remain separate eligibility.
        observation_admission[i]=!allowed ? "visibility_not_admitted" :
            _pf_text(ann.medium)!="underwater" ? "medium_"*_pf_text(ann.medium) :
            (!_pf_finite(ann.x) || !_pf_finite(ann.y)) ? "missing_or_nonfinite_coordinates" :
            observation_policy!==nothing && !valid_radius ? "invalid_original_uncertainty_radius" :
            !valid_bounds ? "invalid_image_or_annotation_bounds" :
            observation_policy===nothing ? "strict_visible" : "provisional_"*original_visibility[i]
        if _pf_text(row.sync_status) != "mapped"
            reason[i]="sync_"*_pf_text(row.sync_status); continue
        elseif _pf_text(ann.medium) != "underwater"
            reason[i]="medium_"*_pf_text(ann.medium); continue
        elseif !allowed || !_pf_finite(ann.x) || !_pf_finite(ann.y)
            reason[i]="annotation_"*_pf_text(ann.visibility); continue
        elseif observation_policy!==nothing && !valid_radius
            reason[i]="invalid_original_uncertainty_radius"; continue
        elseif !valid_bounds
            reason[i]="invalid_image_or_annotation_bounds"; continue
        elseif !_pf_finite(row.fs) || row.fs <= 0 || !_pf_finite(row.recording_time_s) || row.recording_time_s < 0
            reason[i]="invalid_time_or_fs"; continue
        end
        delays=[row[Symbol(method*"_d"*string(j))] for j in 1:3]
        if !all(_pf_finite,delays) || abs(Float64(delays[3])/row.fs) > 1e-12
            reason[i]="invalid_tdoa_or_reference"; continue
        end
        t=Float64(row.recording_time_s); within=mod(t,30.0)
        if within < 2.0 || within >= 28.0
            reason[i]="boundary_purge"; split[i]="purged"; continue
        end
        block=floor(Int,t/30.0)
        push!(rows,i);push!(blocks,block);push!(frame_ids,frame);push!(xs,Float64(ann.x));push!(ys,Float64(ann.y))
        push!(uncertainty,observation_policy===nothing ? max(1.0,_pf_number(ann.uncertainty_px,1.0)) : Float64(ann.uncertainty_px))
        push!(tau,Float64.(delays[1:2])./Float64(row.fs));push!(widths,Int(row.width));push!(heights,Int(row.height))
        split[i]=mod(block,5) == 4 ? "heldout" : "train"; reason[i]="eligible"
    end
    length(unique(widths)) <= 1 && length(unique(heights)) <= 1 || throw(ArgumentError("one image resolution per device required"))
    length(unique(_pf_text.(events.device))) == 1 || throw(ArgumentError("fit_device requires one independent device"))
    # A frame shared by outer train/test is excluded rather than allowing its
    # single annotation into both. Purging normally prevents this already.
    leaked=Set{String}()
    for frame in unique(frame_ids)
        frame_splits=Set(split[rows[findall(==(frame),frame_ids)]])
        length(frame_splits) > 1 && push!(leaked,frame)
    end
    keep=findall(frame->!(frame in leaked),frame_ids)
    for j in setdiff(collect(eachindex(rows)),keep)
        reason[rows[j]]="shared_frame_across_outer_split";split[rows[j]]="ineligible"
    end
    frame_counts=Dict(frame=>count(==(frame),frame_ids[keep]) for frame in unique(frame_ids[keep]))
    weight=[1/(frame_counts[frame_ids[j]]*uncertainty[j]^2) for j in keep]
    geometry=clicker_geometry(); baseline=geometry[1:2,1:2] .- permutedims(geometry[3,1:2])
    tau_matrix=isempty(tau) ? zeros(0,2) : reduce(vcat,permutedims.(tau))
    data=(tau=tau_matrix[keep,:],x=xs[keep],y=ys[keep],weight=weight,block=blocks[keep],
        event_row=rows[keep],frame=frame_ids[keep],inverse_baseline=inv(baseline))
    return (events=events,data=data,reason=reason,split=split,annotated_x=annotated_x,annotated_y=annotated_y,
        width=isempty(widths) ? 0 : first(widths),height=isempty(heights) ? 0 : first(heights),
        device=isempty(events.device) ? "" : _pf_text(first(events.device)),hashes=hashes,
        observation_policy=observation_policy,original_visibility=original_visibility,
        original_uncertainty=original_uncertainty,observation_admission=observation_admission)
end

function _pf_method_comparisons(events,data,calibration)
    comparisons=Dict{String,Any}()
    for method in ("minmax2","minmax","gcc","srp")
        all(in(names(events)),[method*"_d"*string(i) for i in 1:3]) || continue
        dt=fill(NaN,length(data.x),2)
        for (i,rowindex) in enumerate(data.event_row)
            row=events[rowindex,:]
            d=[row[Symbol(method*"_d"*string(j))] for j in 1:3]
            if all(_pf_finite,d) && abs(d[3]/row.fs) <= 1e-12
                dt[i,:]=Float64.(d[1:2])./Float64(row.fs)
            end
        end
        alternate=merge(data,(tau=dt,))
        comparisons[method]=_pf_metrics(calibration,alternate)
    end
    return Dict("interpretation"=>"Post-first-fit descriptive comparison: same baseline-eligible events, same minmax2 physical parameters and offsets, no method refits or heldout method selection. Delays may have method-specific biases; this does not establish tuned heldout performance.",
        "methods"=>comparisons,"recommendation"=>"Inspect same-event coverage/errors for existing minmax and GCC before proposing algorithm changes. Broader algorithms and method-specific tuning require a new training-only selection protocol; optional SRP is compared only when supplied.")
end

function _pf_predictions(input,calibration,evaluation_calibration)
    table=DataFrame(event_id=input.events.event_id,frame_id=input.events.frame_id,
        split=input.split,eligibility_status=input.reason,annotation_x=input.annotated_x,annotation_y=input.annotated_y,
        original_visibility=input.original_visibility,original_uncertainty_px=input.original_uncertainty,
        observation_admission=input.observation_admission)
    for (prefix,cal) in (("final",calibration),("evaluation",evaluation_calibration))
        x=fill(NaN,nrow(input.events));y=copy(x);error=copy(x);valid=falses(length(x));inside=falses(length(x));status=fill("",length(x))
        for (i,row) in enumerate(eachrow(input.events))
            delays=[row[Symbol(cal["method"]*"_d"*string(j))] for j in 1:3]
            result = _pf_finite(row.fs) && row.fs > 0 ? project_tdoas(cal,delays;fs=row.fs) : _pf_invalid("invalid_fs")
            x[i]=result.x;y[i]=result.y;valid[i]=result.valid;inside[i]=result.in_frame;status[i]=result.status
            result.valid && isfinite(input.annotated_x[i]) && isfinite(input.annotated_y[i]) && (error[i]=hypot(result.x-input.annotated_x[i],result.y-input.annotated_y[i]))
        end
        table[!,prefix*"_x"]=x;table[!,prefix*"_y"]=y;table[!,prefix*"_valid"]=valid
        table[!,prefix*"_status"]=status;table[!,prefix*"_in_frame"]=inside;table[!,prefix*"_error_px"]=error
    end
    # Preset/overlay integration aliases. baseline means the UNCALIBRATED D3
    # linear-angle projection; evaluation_* stays the training-only fitted model.
    table.predicted_x=copy(table.final_x);table.predicted_y=copy(table.final_y)
    n=nrow(input.events);bx=fill(NaN,n);by=copy(bx);be=copy(bx);bv=falses(n);bi=falses(n);bs=fill("",n)
    for (i,row) in enumerate(eachrow(input.events))
        delays=[row[Symbol(calibration["method"]*"_d"*string(j))] for j in 1:3]
        raw=_pf_d3_baseline(delays,input.width,input.height;fs=row.fs,units="samples")
        bx[i]=raw.x;by[i]=raw.y;bv[i]=raw.valid;bi[i]=raw.in_frame;bs[i]=raw.status
        raw.valid && isfinite(input.annotated_x[i]) && isfinite(input.annotated_y[i]) && (be[i]=hypot(raw.x-input.annotated_x[i],raw.y-input.annotated_y[i]))
    end
    table.baseline_x=bx;table.baseline_y=by;table.baseline_error_px=be
    table.baseline_valid=bv;table.baseline_in_frame=bi;table.baseline_status=bs
    return table
end

"""
    fit_device(output_dir; method="minmax2", resume=true, nstarts=6,
               maxiters=2500, robust_delta_px=12.0, allow_extended=true,
               diagnostics_profiles=true) -> Dict

Native Julia Optim Float64 fit from events.csv and user/vision annotations.csv.
Absent an explicit fit_observation_policy.json, only medium="underwater" and
visibility="visible" finite annotations qualify. An authenticated opt-in policy
may admit original uncertain coordinates with finite positive original radii;
air/mixed/unknown media cannot be fitted with an underwater sound-speed model.
Requires a hash-valid annotations.complete.json; annotation generation must not
use this model. Baseline first fit is always minmax2. Outer evaluation holds out
block=floor(recording_time_s/30) with block%5==4, purging [0,2) and [28,30)
seconds in EVERY block. Optional optical complexity uses ONLY inner blocked
training validation. Outer evaluation is written/frozen before final deployment
refit on all eligible nonpurged annotations. predictions.csv distinguishes the
training-only evaluation model from the final all-eligible deployment model.
Final errors are descriptive, never heldout accuracy. No empirical pixel warp.
"""
function fit_device(output_dir; method="minmax2",resume=true,nstarts=6,maxiters=2500,
        robust_delta_px=12.0,allow_extended=true,diagnostics_profiles=true)
    try
        method == "minmax2" || throw(ArgumentError("baseline first fit must use minmax2; other methods are descriptive post-fit comparisons"))
        robust_delta_px > 0 || throw(ArgumentError("robust_delta_px must be positive"))
        input=_pf_load(output_dir,method)
        settings=Dict("method"=>method,"nstarts"=>nstarts,"maxiters"=>maxiters,
            "robust_delta_px"=>robust_delta_px,"allow_extended"=>allow_extended,"diagnostics_profiles"=>diagnostics_profiles)
        calibration_path=joinpath(output_dir,"calibration.json"); report_path=joinpath(output_dir,"fit_report.json")
        prediction_path=joinpath(output_dir,"predictions.csv"); frozen_path=joinpath(output_dir,"heldout_evaluation.json")
        products=(calibration_path,report_path,prediction_path,frozen_path)
        if any(isfile,products)
            all(isfile,products) || throw(ArgumentError("partial existing fit products; use a fresh output directory"))
            report=read_json(report_path)
            get(report,"hashes",nothing)==input.hashes || throw(ArgumentError("existing fit input/source/policy hash mismatch; use a fresh output directory"))
            get(report,"settings",nothing)==settings || throw(ArgumentError("existing fit settings mismatch; use a fresh output directory"))
            for (key,path) in (("calibration_sha256",calibration_path),("predictions_sha256",prediction_path),("heldout_evaluation_sha256",frozen_path))
                get(report,key,nothing)==sha256_file(path) || throw(ArgumentError("existing fit product hash mismatch: $key"))
            end
            resume || throw(ArgumentError("existing fit products cannot be overwritten; use a fresh output directory"))
            return report
        end
        data=input.data; train_indices=findall(i->input.split[data.event_row[i]] == "train",eachindex(data.x))
        heldout_indices=findall(i->input.split[data.event_row[i]] == "heldout",eachindex(data.x))
        train=_pf_subset(data,train_indices); heldout=_pf_subset(data,heldout_indices)
        length(train.x) >= 12 || throw(ArgumentError("insufficient eligible training observations (<12); no fit saved"))
        length(unique(train.frame)) >= 8 || throw(ArgumentError("insufficient independent training frames (<8); no fit saved"))
        width,height=input.width,input.height
        # Baseline fits first, before diagnostics/alternative TDOA comparisons.
        baseline=_pf_fit(train,width,height,"centered_pinhole";nstarts=nstarts,maxiters=maxiters,delta=robust_delta_px)
        baseline_cal=_pf_calibration(baseline.q,width,height,"centered_pinhole";device=input.device,method=method)
        comparison=_pf_method_comparisons(input.events,train,baseline_cal)
        model,selection=_pf_select_model(train,width,height;allow_extended=allow_extended,nstarts=nstarts,maxiters=maxiters,delta=robust_delta_px)
        evaluation_fit=model == "centered_pinhole" ? baseline : _pf_fit(train,width,height,model;nstarts=nstarts,maxiters=maxiters,delta=robust_delta_px,initial=baseline.q)
        evaluation_cal=_pf_calibration(evaluation_fit.q,width,height,model;device=input.device,method=method)
        diagnostics=_pf_diagnostics(evaluation_fit.q,train,width,height,model;profiles=diagnostics_profiles,maxiters=min(maxiters,1500),delta=robust_delta_px)
        split_evidence=Dict("time_column"=>"recording_time_s","block_seconds"=>30,"heldout_rule"=>"floor(t/30)%5==4",
            "boundary_purge_seconds"=>2,"purge_rule"=>"mod(t,30)<2 or mod(t,30)>=28",
            "train_event_ids"=>collect(input.events.event_id[train.event_row]),"heldout_event_ids"=>collect(input.events.event_id[heldout.event_row]),
            "train_blocks"=>sort(unique(train.block)),"heldout_blocks"=>sort(unique(heldout.block)),
            "weighting"=>"Each annotation frame total inverse-square uncertainty weight; repeated events do not multiply its training weight",
            "n_all_eligible"=>length(data.x),"n_excluded"=>nrow(input.events)-length(data.x),
            "excluded_reason_counts"=>Dict(reason=>count(==(reason),input.reason) for reason in unique(input.reason) if reason != "eligible"))
        heldout_metrics=_pf_metrics(evaluation_cal,heldout)
        frozen=Dict("schema_version"=>"clicker-heldout-evaluation/1.0","frozen_at_utc"=>string(Dates.now(Dates.UTC))*"Z",
            "hashes"=>input.hashes,"settings"=>settings,"model_selection"=>selection,"evaluation_calibration"=>evaluation_cal,
            "split"=>split_evidence,"training_metrics"=>_pf_metrics(evaluation_cal,train),"heldout_metrics"=>heldout_metrics,
            "evaluation_status"=>isempty(heldout.x) ? "no_eligible_heldout_events" : "frozen_before_final_all_eligible_refit",
            "diagnostics_training_only"=>diagnostics,"optimization_training_only"=>evaluation_fit.starts)
        frozen["uncalibrated_d3_baseline_spec"]=_pf_baseline_spec(width,height)
        frozen["physical_input_screen_training_only"]=_pf_screen_summary(train)
        frozen["paired_baseline_comparison_heldout"]=_pf_paired_baseline(heldout,evaluation_cal,width,height)
        frozen["observation_policy"]=_pf_observation_evidence(input)
        evaluation_cal["observation_policy"]=frozen["observation_policy"]
        write_json(frozen_path,frozen)
        frozen_sha=sha256_file(frozen_path)
        # Only now may heldout observations contribute to the deployment fit.
        final_fit=_pf_fit(data,width,height,model;nstarts=nstarts,maxiters=maxiters,delta=robust_delta_px,initial=evaluation_fit.q)
        calibration=_pf_calibration(final_fit.q,width,height,model;device=input.device,method=method)
        calibration["hashes"]=input.hashes;calibration["settings"]=settings
        calibration["observation_policy"]=frozen["observation_policy"]
        calibration["heldout_evaluation_sha256"]=frozen_sha
        calibration["fit_status"]=get(diagnostics,"weak_identifiability",true) ? "fitted_weak_identifiability" : "fitted"
        calibration["evaluation_status"]=frozen["evaluation_status"]
        calibration["deployment_fit_scope"]="all eligible nonpurged annotations, after frozen outer evaluation; errors are in-sample"
        calibration["parameter_bound_hits_normalized"]=findall(abs.(final_fit.q) .>= 0.98)
        write_json(calibration_path,calibration)
        predictions=_pf_predictions(input,calibration,evaluation_cal)
        tmp=prediction_path*".tmp."*string(getpid());CSV.write(tmp,predictions);mv(tmp,prediction_path;force=true)
        report=Dict{String,Any}("schema_version"=>"clicker-fit-report/1.0","device"=>input.device,"fit_status"=>calibration["fit_status"],
            "hashes"=>input.hashes,"settings"=>settings,"selected_model"=>model,"split"=>split_evidence,
            "baseline_training_metrics"=>_pf_metrics(baseline_cal,train),"heldout_metrics"=>heldout_metrics,
            "final_all_eligible_metrics_descriptive"=>_pf_metrics(calibration,data),"model_selection"=>selection,
            "diagnostics_training_only"=>diagnostics,"post_first_fit_method_comparisons_training_only"=>comparison,
            "physical_input_screen_training_only"=>_pf_screen_summary(train),
            "physical_input_screen_final_all_eligible"=>_pf_screen_summary(data),
            "uncalibrated_d3_baseline_spec"=>_pf_baseline_spec(width,height),
            "paired_baseline_comparison_training"=>_pf_paired_baseline(train,evaluation_cal,width,height),
            "paired_baseline_comparison_heldout"=>frozen["paired_baseline_comparison_heldout"],
            "paired_baseline_comparison_final_descriptive"=>_pf_paired_baseline(data,calibration,width,height),
            "optimization_final_all_eligible"=>final_fit.starts,"heldout_evaluation_sha256"=>frozen_sha,
            "calibration_sha256"=>sha256_file(calibration_path),"predictions_sha256"=>sha256_file(prediction_path),
            "limits"=>["Far-field and fixed array-front branch assumptions; planar array cannot determine range or front/back.",
                "Pixel residuals and local uncertainty do not establish sound-speed truth; c/FOV/delay tradeoffs require coverage and profile inspection.",
                "Final deployment model uses heldout labels after evaluation freeze; only evaluation_* predictions from training-only model provide heldout evidence.",
                "No acoustic predictions are used to generate or choose annotations; annotation SHA changes invalidate dependent fits."])
        report["observation_policy"]=frozen["observation_policy"]
        write_json(report_path,report)
        return report
    catch exception
        log_failure("fit_device $(output_dir)",exception,catch_backtrace())
        rethrow()
    end
end
