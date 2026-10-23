# Separating the clicks of n dolphins (sources) and deciding which follower detections are echoes, from one selected
# direction per detection (e.g. the longest-track rule of src/multiarrival.jl). Works for any number of sources: `n`
# given caps the simultaneously active sources (labels 1..n, reused like track colours), `n = nothing` estimates them.
# plan: reviews/plans/multiarrival-round4-realism-sources_2026-10-24.md
using Statistics, LinearAlgebra, DataFrames, Clustering, Random

_uv(az, el) = [sind(az) * cosd(el), sind(el)]
_angle2(a, b) = (w(x) = sqrt(max(1 - x[1]^2 - x[2]^2, 0.0)); acosd(clamp(a[1] * b[1] + a[2] * b[2] + w(a) * w(b), -1, 1)))

# separate_sources parameters for D3-like directions (selected-direction noise of about 15°), tuned on the calibrated
# realistic pool simulation (improvements/tdoa_localization_2026-10-22/scripts/tune_sources.jl, round 4): pairwise F1
# 0.84 with 2 dolphins and 0.83 with 3 (held out) vs 0.78 / 0.77 with the defaults, which suit low-noise directions
const SEPARATION_PARAMS_D3 = (; sigma_meas_deg=8.0, gate_deg=25.0, force_gate_deg=40.0, t_drop=1.0, follow_gate_deg=12.0, merge_deg=8.0)

mutable struct _Source
    id::Int
    slot::Int
    x::Vector{Float64}     # direction cosines (u, v)
    p::Float64             # isotropic variance
    t_last::Float64
    times::Vector{Float64}
    confirmed::Bool
end

"""
    separate_sources(t, az, el; follower=nothing, echo_window_ms=1.6, n=nothing, gate_deg=8.0, sigma_meas_deg=2.0,
        sigma_q_deg=20.0, gate_chi2=9.21, t_drop=1.0, n_confirm=5, t_confirm=0.5, force_gate_deg=12.0,
        follow_gate_deg=5.0, leader_diff_deg=10.0, repeat_deg=3.0, echo_ref_ms=3.0, merge_deg=4.0)

Assigns each detection (time `t`, s, sorted; selected direction `az`, `el`, deg) to a source.
- A detection within `echo_window_ms` of the previous one is a *follower* (or pass `follower`); the others start
  *events*. Each event goes to the active source whose Kalman-predicted direction (random walk `sigma_q_deg` °/√s,
  measurement noise `sigma_meas_deg`) fits best within `gate_deg` and the chi-square gate, else to a new source.
  With `n` given, at most n sources are active; when all n are busy, an unconfirmed source is replaced, or the event
  joins the nearest source within `force_gate_deg`, or stays unassigned. Sources end after `t_drop` s without events
  and are confirmed by `n_confirm` events within `t_confirm` s; detections of never-confirmed sources become label 0.
  Sources that converge within `merge_deg` are merged into the one with more events.
- A follower is compared with its *reference*: the genuine click (leader, own or repeat) of the last `echo_ref_ms`
  with the closest azimuth (vertical images keep the azimuth), so that each dolphin's echoes go with its own click.
  It is an *echo* (taking the reference's direction and source) unless it repeats the reference's direction (within
  `repeat_deg`: another click of the same source) or it fits another confirmed source within `follow_gate_deg` while
  differing from the reference by more than `leader_diff_deg` (*own*: the other dolphin's click, keeping its direction).
Returns (; label, az, el, role, sources): `label[i]` the source (slot 1..n when `n` is given, else a source id; 0 =
unassigned), `az`/`el` the output direction (echoes take the leader's), `role[i]` "leader", "echo", "repeat", "own" or
"unassigned", and `sources` a DataFrame (id, slot, t_start, t_end, events, confirmed, az_med, el_med, ici_med_ms).
"""
function separate_sources(t, az, el; follower=nothing, echo_window_ms=1.6, n=nothing, gate_deg=8.0, sigma_meas_deg=2.0,
    sigma_q_deg=20.0, gate_chi2=9.21, t_drop=1.0, n_confirm=5, t_confirm=0.5, force_gate_deg=12.0,
    follow_gate_deg=5.0, leader_diff_deg=10.0, repeat_deg=3.0, echo_ref_ms=3.0, merge_deg=4.0)

    N = length(t)
    fol = isnothing(follower) ? [i > 1 && (t[i] - t[i-1]) * 1000 < echo_window_ms for i in 1:N] : collect(follower)
    r = deg2rad(sigma_meas_deg)^2; q = deg2rad(sigma_q_deg)^2
    src_of = zeros(Int, N); oaz = Float64.(copy(az)); oel = Float64.(copy(el)); role = fill("leader", N)
    active = _Source[]; allsrc = _Source[]; next_id = 1
    genuine = Int[]                                   # recent leader / own / repeat detections (echo references)
    alias = Dict{Int,Int}(); root(id) = (while haskey(alias, id); id = alias[id]; end; id)
    cost(s, z, ti) = sum(abs2, z .- s.x) / (s.p + q * (ti - s.t_last) + r)
    function update!(s, z, ti)
        pp = s.p + q * (ti - s.t_last); g = pp / (pp + r)
        s.x .+= g .* (z .- s.x); s.p = (1 - g) * pp; s.t_last = ti; push!(s.times, ti)
        s.confirmed |= count(>=(ti - t_confirm), s.times) >= n_confirm
    end
    function newsource!(z, ti)
        used = Set(s.slot for s in active)
        slot = isnothing(n) ? next_id : something(findfirst(k -> !(k in used), 1:n), 0)
        s = _Source(next_id, slot, copy(z), r, ti, [ti], n_confirm <= 1)
        next_id += 1; push!(active, s); push!(allsrc, s)
        return s
    end
    function merge_converged!()
        a = 1
        while a <= length(active)
            b = findfirst(k -> k != a && _angle2(active[a].x, active[k].x) < merge_deg, eachindex(active))
            if isnothing(b)
                a += 1; continue
            end
            keep, drop = length(active[a].times) >= length(active[b].times) ? (a, b) : (b, a)
            K_, D_ = active[keep], active[drop]
            alias[D_.id] = K_.id
            K_.times = sort([K_.times; D_.times]); K_.confirmed |= D_.confirmed; K_.t_last = max(K_.t_last, D_.t_last)
            deleteat!(active, drop); a = 1
        end
    end
    for i in 1:N
        filter!(s -> t[i] - s.t_last <= t_drop, active)
        filter!(j -> (t[i] - t[j]) * 1000 <= echo_ref_ms, genuine)
        z = _uv(az[i], el[i])
        if fol[i] && !isempty(genuine)
            ref = genuine[argmin([(abs(az[i] - oaz[j]), -j) for j in genuine])]   # closest azimuth, then the latest
            sL = root(src_of[ref])
            dl = _angle2(z, _uv(oaz[ref], oel[ref]))
            if dl <= repeat_deg
                src_of[i] = sL; role[i] = "repeat"; push!(genuine, i)
                s = findfirst(s -> s.id == sL, active); isnothing(s) || update!(active[s], z, t[i])
            else
                cands = [s for s in active if s.confirmed && s.id != sL && _angle2(z, s.x) <= follow_gate_deg]
                if !isempty(cands) && dl > leader_diff_deg
                    s = argmin(s -> cost(s, z, t[i]), cands); update!(s, z, t[i]); src_of[i] = s.id; role[i] = "own"; push!(genuine, i)
                else
                    src_of[i] = sL; oaz[i] = oaz[ref]; oel[i] = oel[ref]; role[i] = "echo"
                end
            end
            continue
        end
        push!(genuine, i)
        best = nothing; bc = Inf
        for s in active
            _angle2(z, s.x) <= gate_deg || continue
            c = cost(s, z, t[i])
            c <= gate_chi2 && c < bc && (bc = c; best = s)
        end
        if !isnothing(best)
            update!(best, z, t[i])
        elseif isnothing(n) || length(active) < n
            best = newsource!(z, t[i])
        else
            unconf = findall(s -> !s.confirmed, active)
            if !isempty(unconf)
                deleteat!(active, unconf[argmin([active[k].t_last for k in unconf])])
                best = newsource!(z, t[i])
            else
                near = [s for s in active if _angle2(z, s.x) <= force_gate_deg]
                isempty(near) || (best = argmin(s -> _angle2(z, s.x), near); update!(best, z, t[i]))
            end
        end
        if isnothing(best)
            role[i] = "unassigned"
        else
            src_of[i] = best.id
            merge_converged!()
        end
    end
    conf = Dict(s.id => s.confirmed for s in allsrc); slot = Dict(s.id => s.slot for s in allsrc)
    src_of .= [id > 0 ? root(id) : 0 for id in src_of]
    allsrc = filter(s -> !haskey(alias, s.id), allsrc)
    label = zeros(Int, N)
    for i in 1:N
        id = src_of[i]
        if id > 0 && conf[id]
            label[i] = isnothing(n) ? id : slot[id]
        else
            role[i] = "unassigned"
        end
    end
    sources = DataFrame(id=[s.id for s in allsrc], slot=[s.slot for s in allsrc], t_start=[s.times[1] for s in allsrc],
        t_end=[s.times[end] for s in allsrc], events=[length(s.times) for s in allsrc], confirmed=[s.confirmed for s in allsrc])
    sources.az_med = [(ix = findall(src_of .== s.id); isempty(ix) ? NaN : median(oaz[ix])) for s in allsrc]
    sources.el_med = [(ix = findall(src_of .== s.id); isempty(ix) ? NaN : median(oel[ix])) for s in allsrc]
    sources.ici_med_ms = [length(s.times) < 2 ? NaN : 1000median(diff(s.times)) for s in allsrc]
    return (; label, az=oaz, el=oel, role, sources)
end

"""
    pairwise_scores(true_label, pred_label, t; max_dt=0.2) -> (; precision, recall, f1, pairs)

"Same source" agreement over pairs of detections at most `max_dt` s apart (robust to label permutations and to
identity resets between bouts): precision = share of predicted same-source pairs that are truly the same source,
recall = share of true same-source pairs predicted as such. Detections with true label 0 are skipped; predicted
label 0 never pairs.
"""
function pairwise_scores(true_label, pred_label, t; max_dt=0.2)
    tp = fp = fn = 0
    for i in eachindex(t)
        true_label[i] > 0 || continue
        j = i + 1
        while j <= lastindex(t) && t[j] - t[i] <= max_dt
            if true_label[j] > 0
                st = true_label[i] == true_label[j]; sp = pred_label[i] > 0 && pred_label[i] == pred_label[j]
                tp += st & sp; fp += !st & sp; fn += st & !sp
            end
            j += 1
        end
    end
    pr = tp / max(tp + fp, 1); rc = tp / max(tp + fn, 1)
    return (; precision=pr, recall=rc, f1=2pr * rc / max(pr + rc, eps()), pairs=tp + fp + fn)
end

# adjusted Rand index over detections with a true label (predicted 0 counts as one more class)
adjusted_rand(true_label, pred_label) = (ok = true_label .> 0; randindex(true_label[ok], pred_label[ok] .+ 1)[1])

"""
    kmeans_sources(t, az, el, k; follower=nothing, echo_window_ms=1.6, seed=1)

Baseline: k-means on the direction cosines of the event leaders (time ignored); followers take their leader's label.
"""
function kmeans_sources(t, az, el, k; follower=nothing, echo_window_ms=1.6, seed=1)
    N = length(t)
    fol = isnothing(follower) ? [i > 1 && (t[i] - t[i-1]) * 1000 < echo_window_ms for i in 1:N] : collect(follower)
    lead = findall(.!fol)
    X = reduce(hcat, [_uv(az[i], el[i]) for i in lead])
    lab = zeros(Int, N)
    length(lead) >= k || return lab
    a = assignments(kmeans(X, k; rng=Random.Xoshiro(seed)))
    lab[lead] = a
    for i in 2:N; fol[i] && (lab[i] = lab[i-1]); end
    return lab
end
