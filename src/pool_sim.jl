# Image-source simulation of dolphin clicks in a rectangular pool as recorded by the D3 array, to benchmark the
# multi-arrival methods (src/multiarrival.jl) against ground truth. Not used by the pipeline.
# plan: reviews/plans/multiarrival-round3-pool-simulation_2026-10-23.md
# pool frame: x along the length, y across, z depth (0 at the surface, positive down)
# array frame as in cost_tdoa2ang / tdoa2uv: x right, y up, z forward; the array is level and faces +x of the pool
using LinearAlgebra, Random, Statistics, DataFrames, FFTW

Base.@kwdef struct Pool
    L::Float64 = 20.0             # length (x, m)
    W::Float64 = 10.0             # width (y, m)
    D::Float64 = 4.0              # water depth (z, m)
    refl_surface::Float64 = -0.95 # pressure-release surface: polarity inverted
    refl_bottom::Float64 = 0.6
    refl_wall::Float64 = 0.5
end

_pool_planes(p::Pool) = [(3, 0.0, p.refl_surface), (3, p.D, p.refl_bottom), (1, 0.0, p.refl_wall), (1, p.L, p.refl_wall),
    (2, 0.0, p.refl_wall), (2, p.W, p.refl_wall)]

"""
    pool_images(src, heading, pool; order=3)

Image sources of a source at `src` with beam axis `heading` (pool frame), up to `order` reflections: a Vector of
(; pos, heading, refl, order). Position and heading are mirrored together; `refl` is the product of the reflection
coefficients. An image reached by reflections in a different order is kept once.
"""
function pool_images(src, heading, pool::Pool; order=3)
    planes = _pool_planes(pool)
    first_img = (; pos=Float64.(collect(src)), heading=Float64.(collect(heading)), refl=1.0, order=0, last=0)
    imgs = [first_img]; front = [first_img]
    seen = Set([round.(first_img.pos; digits=6)])
    for o in 1:order
        nxt = typeof(first_img)[]
        for im in front, (k, (ax, v, r)) in enumerate(planes)
            k == im.last && continue
            p = copy(im.pos); p[ax] = 2v - p[ax]
            key = round.(p; digits=6)
            key in seen && continue
            push!(seen, key)
            h = copy(im.heading); h[ax] = -h[ax]
            push!(nxt, (; pos=p, heading=h, refl=im.refl * r, order=o, last=k))
        end
        append!(imgs, nxt); front = nxt
    end
    return [(; pos=im.pos, heading=im.heading, refl=im.refl, order=im.order) for im in imgs]
end

# click beam: Gaussian main lobe (θ0 = 20°: -3 dB at about 12° off axis) above a floor (-30 dB)
beam_gain(θdeg; θ0=20.0, floor=0.03) = max(exp(-(θdeg / θ0)^2), floor)

# unit vector from the array towards point p, in the array frame (x right = -y pool, y up = -z pool, z forward = +x pool)
function array_dir(p, array_pos)
    d = p .- array_pos; d = d ./ norm(d)
    return (-d[2], -d[3], d[1])
end

# az, el (deg) as the planar array perceives a direction (front hemisphere, as tdoa2dir_ls)
perceived_azel(dir) = (el = asind(clamp(dir[2], -1, 1)); (cosd(el) > 1e-12 ? asind(clamp(dir[1] / cosd(el), -1, 1)) : 0.0, el))

"""
    click_arrivals(src, heading, array_pos, pool; order=3, c=default_c, θ0=20.0, floor=0.03, rel_min=0.01)

Arrivals of one click at the array centre: a Vector of (; delay (s), gain, dir (array frame), order), with
gain = reflection product x beam gain toward the array / path length. Arrivals weaker than `rel_min` x the strongest
are dropped.
"""
function click_arrivals(src, heading, array_pos, pool::Pool; order=3, c=default_c, θ0=20.0, floor=0.03, rel_min=0.01)
    out = NamedTuple{(:delay, :gain, :dir, :order),Tuple{Float64,Float64,NTuple{3,Float64},Int}}[]
    for im in pool_images(src, heading, pool; order=order)
        d = array_pos .- im.pos; r = norm(d)
        θ = acosd(clamp(dot(im.heading, d) / (r * norm(im.heading)), -1, 1))
        push!(out, (; delay=r / c, gain=im.refl * beam_gain(θ; θ0=θ0, floor=floor) / r, dir=array_dir(im.pos, array_pos), order=im.order))
    end
    gmax = maximum(a -> abs(a.gain), out)
    return filter(a -> abs(a.gain) >= rel_min * gmax, out)
end

# position folded into [lo, hi] (bouncing off the ends) and the sign of the folded velocity
function _fold(x, lo, hi)
    span = hi - lo
    y = mod(x - lo, 2span)
    return y <= span ? (lo + y, 1.0) : (hi - (y - span), -1.0)
end

"""
    random_source(rng, pool, array_pos, dur; speed=(0.5, 2.0), margin=0.5, p_burst=0.5, p_face=0.3, level_db=4.0, ntemplates=1)

A source moving at constant speed (bouncing off the pool walls, surface and bottom) and its click schedule:
trains of burst pulses (ICI 1.8-2.5 ms, 0.1-0.6 s) or echolocation clicks (ICI 8-40 ms, drifting, 0.3-1.5 s) separated
by 0.05-0.5 s. Per train: a level (±`level_db` dB), a template index and a beam direction, either along the
movement turned by random yaw (σ 30°) and pitch (σ 10°) offsets or, with probability `p_face`, towards the array.
Returns (; pos(t), heading(t, click), clicks) where `clicks` holds (; t, level, tmpl, yaw, pitch, face, burst).
"""
function random_source(rng, pool::Pool, array_pos, dur; speed=(0.5, 2.0), margin=0.5, p_burst=0.5, p_face=0.3, level_db=4.0, ntemplates=1)
    U(a, b) = a + (b - a) * rand(rng)
    lo = [margin, margin, margin]; hi = [pool.L - margin, pool.W - margin, pool.D - margin]
    p0 = [U(min(array_pos[1] + 1.5, hi[1]), min(array_pos[1] + 12, hi[1])), U(lo[2], hi[2]), U(lo[3], hi[3])]
    φ = U(0, 360); s = U(speed...)
    v = [s * cosd(φ), s * sind(φ), 0.2s * randn(rng)]
    pos(t) = [first(_fold(p0[a] + v[a] * t, lo[a], hi[a])) for a in 1:3]
    vel(t) = [v[a] * last(_fold(p0[a] + v[a] * t, lo[a], hi[a])) for a in 1:3]
    clicks = NamedTuple{(:t, :level, :tmpl, :yaw, :pitch, :face, :burst),Tuple{Float64,Float64,Int,Float64,Float64,Bool,Bool}}[]
    t = U(0, 0.3)
    while t < dur
        burst = rand(rng) < p_burst
        t_end = min(t + (burst ? U(0.1, 0.6) : U(0.3, 1.5)), dur)
        ici = burst ? U(1.8e-3, 2.5e-3) : U(8e-3, 40e-3)
        level = 10^(level_db * randn(rng) / 20); tmpl = rand(rng, 1:ntemplates)
        yaw, pitch = 30randn(rng), 10randn(rng); face = rand(rng) < p_face
        while t < t_end
            push!(clicks, (; t, level, tmpl, yaw, pitch, face, burst))
            ici = burst ? clamp(ici * (1 + 0.03randn(rng)), 1.5e-3, 3e-3) : clamp(ici * (1 + 0.02randn(rng)), 5e-3, 60e-3)
            t += ici
        end
        t += U(0.05, 0.5)
    end
    function heading(t, cl)
        base = cl.face ? array_pos .- pos(t) : vel(t)
        φh = atand(base[2], base[1]); ψh = asind(clamp(-base[3] / norm(base), -1, 1))   # horizontal angle, elevation (up)
        φh += cl.face ? 5randn(rng) : cl.yaw; ψh += cl.face ? 5randn(rng) : cl.pitch
        return [cosd(ψh) * cosd(φh), cosd(ψh) * sind(φh), -sind(ψh)]
    end
    return (; pos, heading, clicks)
end

# fractional delay (samples) of a short zero-padded segment via the FFT
function _sim_frac_delay(x, d)
    n = length(x); N = nextprod([2, 3, 5], 2n)
    X = rfft([x; zeros(N - n)])
    return irfft(X .* cis.(-2π .* (0:length(X)-1) ./ N .* d), N)[1:n]
end

"""
    simulate_scenario(rng, templates, noise; fs=192000, dur=4.0, nsrc=2, pool=Pool(), array_pos=[1.5, 5.0, 1.5], A0=1.0,
        rx, c=default_c, order=3, θ0=20.0, beam_floor=0.03, rel_min=0.01, source_kwargs=(;))

Simulates `nsrc` moving clicking sources (random_source) in `pool`, recorded by the array at `array_pos` (geometry
`rx`, 3 x 3, as rx_vect's hydrophone columns). Every click is the template `templates[tmpl]` (peak at index 41) times
A0 x level x gain of each arrival (click_arrivals), placed on each hydrophone at its plane-wave delay; `noise` (a
matrix with >= dur x fs rows, 3 columns, e.g. real quiet recording) is added from a random offset.
Returns (; data, arrivals, clicks): `data` has 5 columns (1-2 zero, 3-5 the hydrophones, as the D3 WAV files);
`arrivals` one row per arrival (src, click, order, sample = peak position on the reference hydrophone (1-based),
amp, az, el as the planar array perceives them, behind); `clicks` one row per click (src, click, t, x, y, z, beam
angle to the array, level, burst, the direct arrival's az/el/amp).
"""
function simulate_scenario(rng, templates, noise; fs=192000, dur=4.0, nsrc=2, pool=Pool(), array_pos=[1.5, 5.0, 1.5], A0=1.0,
    rx, c=default_c, order=3, θ0=20.0, beam_floor=0.03, rel_min=0.01, ref=3, source_kwargs=(;))

    n = round(Int, dur * fs) + 2000
    pad = 8; pk0 = 41
    banks = [[_sim_frac_delay([zeros(pad); tm; zeros(pad)], f / 32) for f in 0:31] for tm in templates]
    hyd = zeros(n, 3)
    arrivals = DataFrame(src=Int[], click=Int[], order=Int[], sample=Float64[], amp=Float64[], az=Float64[], el=Float64[], behind=Bool[])
    clicks = DataFrame(src=Int[], click=Int[], t=Float64[], x=Float64[], y=Float64[], z=Float64[], beam_angle=Float64[], level=Float64[],
        burst=Bool[], direct_az=Float64[], direct_el=Float64[], direct_amp=Float64[])
    for s in 1:nsrc
        src = random_source(rng, pool, array_pos, dur; ntemplates=length(templates), source_kwargs...)
        for (ci, cl) in enumerate(src.clicks)
            p = src.pos(cl.t); h = src.heading(cl.t, cl)
            arr = click_arrivals(p, h, array_pos, pool; order=order, c=c, θ0=θ0, floor=beam_floor, rel_min=rel_min)
            dd = array_pos .- p
            for a in arr
                τ = -(rx' * collect(a.dir)) ./ c .* fs               # per-hydrophone delay (samples), as get_tdoa_srp
                amp = A0 * cl.level * a.gain
                for ch in 1:3
                    x = (cl.t + a.delay) * fs + τ[ch] + 1             # position of the template peak (1-based)
                    ip = floor(Int, x); f = round(Int, (x - ip) * 32)
                    f == 32 && (ip += 1; f = 0)
                    seg = banks[cl.tmpl][f+1]
                    i0 = ip - pk0 - pad + 1
                    (i0 < 1 || i0 + length(seg) - 1 > n) && continue
                    @views hyd[i0:i0+length(seg)-1, ch] .+= amp .* seg
                end
                az, el = perceived_azel(a.dir)
                push!(arrivals, (s, ci, a.order, (cl.t + a.delay) * fs + τ[ref] + 1, amp, az, el, a.dir[3] < 0))
            end
            d0 = arr[findfirst(a -> a.order == 0, arr)]
            az0, el0 = perceived_azel(d0.dir)
            push!(clicks, (s, ci, cl.t, p..., acosd(clamp(dot(h, dd) / norm(dd), -1, 1)), cl.level, cl.burst, az0, el0, A0 * cl.level * d0.gain))
        end
    end
    off = rand(rng, 0:size(noise, 1)-n)
    hyd .+= @view noise[off+1:off+n, :]
    return (; data=[zeros(n, 2) hyd], arrivals, clicks)
end

"""
    random_pool_setup(rng; L=(15, 30), W=(8, 15), D=(3, 7)) -> (pool, array_pos)

A random rectangular pool and an array position near the x = 0 wall, facing into the pool (the real array is
portable, so geometry varies between recordings).
"""
function random_pool_setup(rng; L=(15.0, 30.0), W=(8.0, 15.0), D=(3.0, 7.0))
    U(a, b) = a + (b - a) * rand(rng)
    pool = Pool(L=U(L...), W=U(W...), D=U(D...))
    array_pos = [U(0.5, 3.0), U(0.3, 0.7) * pool.W, U(0.5, min(3.0, pool.D - 0.5))]
    return pool, array_pos
end
