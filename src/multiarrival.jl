# Multiple arrivals per click: SRP-PHAT direction candidates (top-K local maxima instead of the global one),
# per-click quality features, follower labels and a light multi-target tracker over the candidates.
# plan: reviews/plans/multiarrival-candidates-tracking_2026-10-22.md, notes: improvements/tdoa_localization_2026-10-22/
# usage: include("run_example.jl"); set_device__ophk_acoustic_D3(); include("run_multitdoa.jl"); include("multiarrival.jl")
# existing functions are untouched; this file only adds new ones
using FFTW, LinearAlgebra, Statistics, Printf, CSV, DataFrames, JLD2, WAV, DSP, Peaks

# ---- SRP candidates ------------------------------------------------------------------------------------------
# (u, v) grid of get_tdoa_srp plus its lattice position, for 8-neighbour local maxima
function _srp_lattice(grid_deg)
    h = sind(grid_deg)
    us = -1:h:1; vs = -1:h:1
    grid = Tuple{Float64,Float64}[]; ij = Tuple{Int,Int}[]
    lat = zeros(Int, length(us), length(vs))
    for (i, u) in enumerate(us), (j, v) in enumerate(vs)
        u^2 + v^2 <= 1 || continue
        push!(grid, (u, v)); push!(ij, (i, j)); lat[i, j] = length(grid)
    end
    return h, grid, ij, lat
end

_uvw(u, v) = (u, v, sqrt(max(1 - u^2 - v^2, 0.0)))
_angle_uv(a, b) = acosd(clamp(sum(_uvw(a...) .* _uvw(b...)), -1, 1))

# delay x by d samples (fractional, circular) via the FFT
function _frac_delay(x, d)
    n = length(x); N = nextfastfft(2n)
    X = rfft([x; zeros(N - n)])
    f = (0:length(X)-1) ./ N
    return irfft(X .* cis.(-2π .* f .* d), N)[1:n]
end

"""
    get_srp_candidates(data, peaks; window, ref_channel, fs, rx_vect, c, weighting=:phat, band=nothing, interp=8,
        grid_deg=1.0, refine=10, K=3, min_sep_deg=8.0, min_rel=0.5, taper=0.0)

SRP-PHAT as `get_tdoa_srp` (same GCC curves, grid and refinement, so candidate 1 is its answer), but returning up
to `K` local maxima of the SRP map per click, at least `min_sep_deg` apart and scoring at least `min_rel` x the best.
Returns a NamedTuple with, per click i and candidate k (NaN / 0 where fewer were found):
- `uv[i, k, :]` direction cosines (u, v) = (sin(az)cos(el), sin(el)); `score[i, k]` SRP (mean pair GCC, 1 = all in phase);
- `tdoa[i, k, :]` TDOAs predicted there (samples, relative to `ref_channel`, as get_tdoa_srp);
- `tarr[i, k]` beam arrival time: channels aligned to the candidate and summed, envelope peak, samples from the click;
- `nfound[i]`; per click quality: `closure` |τ12+τ23-τ13| of the per-pair GCC peaks (samples), `ncc_min` the smallest
  time-domain normalized cross-correlation over the pairs at their GCC lag, `gcc_peak_min` the smallest pair GCC peak,
  `gcc_npeaks` the largest number of GCC local maxima >= half the peak in one pair.
"""
function get_srp_candidates(data, peaks; window=default_window, ref_channel=ref_channel, fs, rx_vect=rx_vect, c=default_c,
    weighting=:phat, band=nothing, interp=8, grid_deg=1.0, refine=10, K=3, min_sep_deg=8.0, min_rel=0.5, taper=0.0)

    m = size(data, 2)
    rx = _rx_for(data, rx_vect)
    τmax = tau_max_samples(rx, fs; c=c)
    band = default_gcc_band(band, fs)
    pairs = [(j, k) for j in 1:m for k in j+1:m]
    arrival(u, v) = -(rx' * [u, v, sqrt(max(1 - u^2 - v^2, 0.0))]) ./ c .* fs   # as in get_tdoa_srp
    pairdelays(u, v) = (t = arrival(u, v); [t[j] - t[k] for (j, k) in pairs])
    h, grid, ij, lat = _srp_lattice(grid_deg)
    D = reduce(hcat, [pairdelays(u, v) for (u, v) in grid])
    fine = range(-h, h; length=2refine + 1)
    nb = [(di, dj) for di in -1:1 for dj in -1:1 if (di, dj) != (0, 0)]

    n = length(peaks)
    uv = fill(NaN, n, K, 2); score = fill(NaN, n, K); tdoa = fill(NaN, n, K, m); tarr = fill(NaN, n, K)
    nfound = zeros(Int, n); closure = fill(NaN, n); ncc_min = fill(NaN, n); gcc_peak_min = fill(NaN, n); gcc_npeaks = zeros(Int, n)
    Threads.@threads for i in eachindex(peaks)
        win = _clip_window(window .+ peaks[i], size(data, 1))
        snip = Float64.(data[win, :])
        taper > 0 && (snip .*= tukey(length(win), taper))
        curves = [gcc_curve(snip[:, j], snip[:, k]; fs=fs, maxlag=τmax[j, k], weighting=weighting, band=band, interp=interp) for (j, k) in pairs]
        sc(d) = sum(p -> curve_at(curves[p], d[p]), eachindex(pairs))
        S = [sc(view(D, :, g)) for g in eachindex(grid)]

        # local maxima over the 8 lattice neighbours, best first (ties keep grid order, like argmax)
        ismax(g) = all(((di, dj),) -> (a = ij[g][1] + di; b = ij[g][2] + dj;
            !(checkbounds(Bool, lat, a, b)) || lat[a, b] == 0 || S[lat[a, b]] <= S[g]), nb)
        cand = sort(filter(ismax, eachindex(grid)); by=g -> -S[g])
        isempty(cand) && (cand = [argmax(S)])
        chosen = Int[]
        for g in cand
            length(chosen) >= K && break
            S[g] >= min_rel * S[cand[1]] || break
            all(q -> _angle_uv(grid[g], grid[q]) >= min_sep_deg, chosen) && push!(chosen, g)
        end
        for (k, g) in enumerate(chosen)
            best, (bu, bv) = S[g], grid[g]
            u0, v0 = bu, bv
            for du in fine, dv in fine
                u, v = u0 + du, v0 + dv
                u^2 + v^2 <= 1 || continue
                s = sc(pairdelays(u, v))
                s > best && ((best, bu, bv) = (s, u, v))
            end
            t = arrival(bu, bv)
            uv[i, k, :] = [bu, bv]; score[i, k] = best / length(pairs); tdoa[i, k, :] = t .- t[ref_channel]
            beam = sum(ch -> _frac_delay(snip[:, ch], -(t[ch] - t[ref_channel])), 1:m)
            tarr[i, k] = win[argmax(abs.(hilbert(beam)))] - peaks[i]
        end
        nfound[i] = length(chosen)

        # per-pair GCC peaks: closure, peak heights, peak count, waveform similarity at the GCC lag
        τ = zeros(length(pairs)); pk = zeros(length(pairs)); npk = zeros(Int, length(pairs)); ncc = zeros(length(pairs))
        for (p, (j, k)) in enumerate(pairs)
            lags, r = curves[p]
            kf, pk[p] = subsample_peak(r, argmax(r))
            τ[p] = first(lags) + (kf - 1) * step(lags)
            npk[p] = count(q -> r[q] >= r[q-1] && r[q] >= r[q+1] && r[q] >= 0.5pk[p], 2:length(r)-1)
            x, y = snip[:, j], snip[:, k]   # x[n] ≈ y[n - τ]: compare x with y delayed by τ (fractional)
            ncc[p] = abs(sum(x .* _frac_delay(y, τ[p]))) / sqrt(max(sum(abs2, x) * sum(abs2, y), floatmin(Float64)))
        end
        closure[i] = m == 3 ? abs(τ[1] + τ[3] - τ[2]) : NaN    # pairs (1,2), (1,3), (2,3)
        gcc_peak_min[i] = minimum(pk); gcc_npeaks[i] = maximum(npk); ncc_min[i] = minimum(ncc)
    end
    return (; uv, score, tdoa, tarr, nfound, closure, ncc_min, gcc_peak_min, gcc_npeaks)
end

# envelope peaks (>= rel x the channel's highest one, min_sep samples apart) inside the window, per click and channel
function envelope_arrivals(data, peaks; window=default_window, rel=0.3, min_sep=8)
    out = zeros(Int, length(peaks), size(data, 2))
    Threads.@threads for i in eachindex(peaks)
        win = _clip_window(window .+ peaks[i], size(data, 1))
        for ch in axes(data, 2)
            pk = findmaxima(abs.(hilbert(Float64.(data[win, ch]))), min_sep)
            out[i, ch] = isempty(pk.heights) ? 0 : count(pk.heights .>= rel * maximum(pk.heights))
        end
    end
    return out
end

# gap to the previous / next detection (ms) and follower = within echo_window_ms after the previous one
function follower_flags(samples, fs; echo_window_ms=1.6)
    n = length(samples)
    gap_prev = [i == 1 ? Inf : (samples[i] - samples[i-1]) / fs * 1000 for i in 1:n]
    gap_next = [i == n ? Inf : (samples[i+1] - samples[i]) / fs * 1000 for i in 1:n]
    return (; gap_prev, gap_next, follower=gap_prev .< echo_window_ms)
end

# ---- tracker --------------------------------------------------------------------------------------------------
mutable struct _Track
    id::Int
    x::Vector{Float64}     # (u, v)
    p::Float64             # isotropic variance of x
    t_last::Float64
    hits::Vector{Float64}  # hit times
    confirmed::Bool
end

"""
    track_candidates(t, uv, score; gate_deg=6, sigma_meas_deg=1.5, sigma_q_deg=20, n_confirm=3, t_confirm=0.05,
        t_delete=0.25, t_delete_tentative=0.03, birth_min_score=0.2, gate_chi2=9.21)

Links direction candidates over time. `t[i]` click time (s, sorted), `uv[i, k, :]`, `score[i, k]` as from
get_srp_candidates. Each track is a Kalman filter on (u, v) with a random walk (`sigma_q_deg` deg/√s) and
measurement noise `sigma_meas_deg`. Per click, (track, candidate) pairs within `gate_deg` and the chi-square gate
are assigned greedily by distance (one candidate per track; one click's candidates may feed different tracks,
e.g. a direct path and its echo). Unassigned candidates scoring >= `birth_min_score` start tentative tracks; a track
is confirmed once it has `n_confirm` hits within `t_confirm` s, deleted after `t_delete` s (tentative:
`t_delete_tentative`) without hits. Two tracks closer than `merge_deg` are merged into the one with more hits (a
clutter candidate can steal a track's update, so the true path starts a duplicate). Returns (; track, filt, tracks): `track[i, k]` (0 = none), `filt[i, k, :]` the
filtered (u, v) after that update, `tracks` a DataFrame (id, t_start, t_end, hits, confirmed).
"""
function track_candidates(t, uv, score; gate_deg=6.0, sigma_meas_deg=1.5, sigma_q_deg=20.0, n_confirm=3, t_confirm=0.05,
    t_delete=0.25, t_delete_tentative=0.03, birth_min_score=0.2, gate_chi2=9.21, merge_deg=3.0, first_id=1)

    n, K = size(score)
    r = deg2rad(sigma_meas_deg)^2; q = deg2rad(sigma_q_deg)^2
    track = zeros(Int, n, K); filt = fill(NaN, n, K, 2)
    active = _Track[]; done = _Track[]; next_id = first_id
    alias = Dict{Int,Int}()
    root(id) = (while haskey(alias, id); id = alias[id]; end; id)
    for i in 1:n
        # drop stale tracks
        keep = _Track[]
        for tr in active
            (t[i] - tr.t_last > (tr.confirmed ? t_delete : t_delete_tentative)) ? push!(done, tr) : push!(keep, tr)
        end
        active = keep
        ks = [k for k in 1:K if !isnan(score[i, k])]
        pairs_ = Tuple{Float64,Int,Int}[]
        for (a, tr) in enumerate(active), k in ks
            z = (uv[i, k, 1], uv[i, k, 2])
            _angle_uv(z, (tr.x[1], tr.x[2])) <= gate_deg || continue
            d2 = ((z[1] - tr.x[1])^2 + (z[2] - tr.x[2])^2) / (tr.p + q * (t[i] - tr.t_last) + r)
            d2 <= gate_chi2 && push!(pairs_, (d2, a, k))
        end
        used_t = falses(length(active)); used_k = falses(K)
        for (_, a, k) in sort(pairs_)
            (used_t[a] || used_k[k]) && continue
            used_t[a] = used_k[k] = true
            tr = active[a]
            pp = tr.p + q * (t[i] - tr.t_last)
            g = pp / (pp + r)
            tr.x .+= g .* ([uv[i, k, 1], uv[i, k, 2]] .- tr.x); tr.p = (1 - g) * pp
            tr.t_last = t[i]; push!(tr.hits, t[i])
            tr.confirmed |= count(>=(t[i] - t_confirm), tr.hits) >= n_confirm
            track[i, k] = tr.id; filt[i, k, :] = tr.x
        end
        for k in ks
            (used_k[k] || score[i, k] < birth_min_score) && continue
            tr = _Track(next_id, [uv[i, k, 1], uv[i, k, 2]], r, t[i], [t[i]], n_confirm <= 1)
            push!(active, tr); next_id += 1
            track[i, k] = tr.id; filt[i, k, :] = tr.x
        end
        # merge converged tracks into the one with more hits
        a = 1
        while a <= length(active)
            b = findfirst(bb -> bb != a && _angle_uv(Tuple(active[a].x), Tuple(active[bb].x)) < merge_deg, eachindex(active))
            if isnothing(b)
                a += 1; continue
            end
            keep, drop = length(active[a].hits) >= length(active[b].hits) ? (a, b) : (b, a)
            K_, D_ = active[keep], active[drop]
            alias[D_.id] = K_.id
            K_.hits = sort([K_.hits; D_.hits]); K_.confirmed |= D_.confirmed; K_.t_last = max(K_.t_last, D_.t_last)
            deleteat!(active, drop)
            a = 1
        end
    end
    track .= root.(track) .* (track .> 0)
    all_tr = sort(filter(tr -> !haskey(alias, tr.id), [done; active]); by=tr -> tr.id)
    tracks = DataFrame(id=[tr.id for tr in all_tr], t_start=[tr.hits[1] for tr in all_tr], t_end=[tr.hits[end] for tr in all_tr],
        hits=[length(tr.hits) for tr in all_tr], confirmed=[tr.confirmed for tr in all_tr])
    return (; track, filt, tracks)
end

# Okabe-Ito hues plus white: distinct on video and for colour-blind viewers
const TRACK_PALETTE = [(0.902, 0.624, 0.0), (0.337, 0.706, 0.914), (0.0, 0.620, 0.451), (0.941, 0.894, 0.259),
    (0.0, 0.447, 0.698), (0.835, 0.369, 0.0), (0.800, 0.475, 0.655), (1.0, 1.0, 1.0)]

# colour index per track (rows of `tracks`, needs t_start, t_end): the first palette entry no overlapping track uses
function track_colours(tracks; npalette=length(TRACK_PALETTE))
    col = zeros(Int, nrow(tracks))
    for i in sortperm(tracks.t_start)
        busy = Set(col[j] for j in 1:nrow(tracks) if col[j] > 0 && tracks.t_start[j] <= tracks.t_end[i] && tracks.t_end[j] >= tracks.t_start[i])
        free = findfirst(c -> !(c in busy), 1:npalette)
        col[i] = something(free, mod1(i, npalette))   # more concurrent tracks than colours: reuse
    end
    return col
end

uv2azel(u, v) = (el = asin(clamp(v, -1, 1)); (cos(el) > 1e-12 ? asin(clamp(u / cos(el), -1, 1)) : 0.0, el))   # as tdoa2dir_ls, radians

# ---- post-processing driver ----------------------------------------------------------------------------------
# filtered relevant channels of samples rng of a WAV file, as detect_impulse filters the whole file
function read_filtered_block(aufname, rng; fs)
    x, _ = wavread(aufname; subrange=rng, format="native")
    return filter_simple(Float64.(x[:, get_relevant_channels(rx_vect)]), impulsive_band_pass; fs=fs)
end

"""
    run_multiarrival(res_folder; out_dir, viewer_dir=res_folder, K=3, block_s=10.0, pad_s=0.05, echo_window_ms=1.6,
        closure_max=1.0, candidate_kwargs=(;), tracker_kwargs=(;))

Post-processes a run_contiguous_folders_multitdoa result folder (no re-detection, no video work): direction candidates,
features and tracks per click. Writes `out_dir/data/click_features_<folder>.csv`, `out_dir/data/tracks_<folder>.csv`,
`out_dir/multiarrival_summary_<folder>.md/.csv`, and in `viewer_dir` the webui files Impulse7_srpTracks (one row per
click and confirmed track, coloured by track) and Impulse8_srpCandidates (all candidates; rank > 1 as triangles).
"""
function run_multiarrival(res_folder; out_dir, viewer_dir=res_folder, K=3, block_s=10.0, pad_s=0.05, echo_window_ms=1.6,
    closure_max=1.0, candidate_kwargs=(;), tracker_kwargs=(;))

    folder = basename(rstrip(res_folder, '/'))
    mkpath(joinpath(out_dir, "data"))
    files = joinpath.(res_folder, sort(filter(f -> occursin("__cps", f) && endswith(f, ".jld2") && !endswith(f, "_angles.jld2"), readdir(res_folder))))
    w6 = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__Impulse6_srpPhat.csv"), DataFrame)
    w5 = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__Impulse5_gccPhatGuided.csv"), DataFrame)
    rx = rx_vect[:, get_relevant_channels(rx_vect)]
    parts = DataFrame[]; trk_parts = DataFrame[]; cands = []; g0 = 0; next_id = 1
    for (fi, f) in enumerate(files)
        j = jldopen(f); tr = j["res_impulsetrain"]; fs = j["fs"]; ri = j["res_impulse"]; close(j)
        aufname, len_data = ri.aufname, ri.len_data
        P = tr.pind_good; n = length(P)
        @info "multiarrival: $(basename(aufname)), $n clicks"
        # blocks: consecutive clicks of a train, at most block_s long, padded for the filter
        st = [tr.train_start_ind; n + 1]; blocks = UnitRange{Int}[]
        for t in 1:length(tr.train_start_ind)
            a = st[t]
            while a < st[t+1]
                b = a
                while b + 1 < st[t+1] && P[b+1] - P[a] <= block_s * fs; b += 1; end
                push!(blocks, a:b); a = b + 1
            end
        end
        C = (; uv=fill(NaN, n, K, 2), score=fill(NaN, n, K), tdoa=fill(NaN, n, K, 3), tarr=fill(NaN, n, K), nfound=zeros(Int, n),
            closure=fill(NaN, n), ncc_min=fill(NaN, n), gcc_peak_min=fill(NaN, n), gcc_npeaks=zeros(Int, n), env=zeros(Int, n, 3))
        pad = round(Int, pad_s * fs)
        for bl in blocks
            lo = max(P[bl[1]] - pad, 1)
            d = read_filtered_block(aufname, lo:min(P[bl[end]] + pad, len_data); fs=fs)
            pk = P[bl] .- lo .+ 1
            c = get_srp_candidates(d, pk; window=window_impulsive, ref_channel=ref_channel, fs=fs, rx_vect=rx, K=K, candidate_kwargs...)
            for nm in (:uv, :tdoa); getfield(C, nm)[bl, :, :] = getfield(c, nm); end
            for nm in (:score, :tarr); getfield(C, nm)[bl, :] = getfield(c, nm); end
            for nm in (:nfound, :closure, :ncc_min, :gcc_peak_min, :gcc_npeaks); getfield(C, nm)[bl] = getfield(c, nm); end
            C.env[bl, :] = envelope_arrivals(d, pk; window=window_impulsive)
        end
        ff = follower_flags(P, fs; echo_window_ms=echo_window_ms)
        tsec = (P .- 1) ./ fs
        T = track_candidates(tsec, C.uv, C.score; first_id=next_id, tracker_kwargs...)
        next_id = isempty(T.tracks.id) ? next_id : maximum(T.tracks.id) + 1
        T.tracks.file .= fi
        push!(trk_parts, T.tracks)
        rows = DataFrame(file=basename(aufname), click=1:n, sample=P, t_s=tsec, frame=w6.frame[g0 .+ (1:n)],
            train=[searchsortedlast(tr.train_start_ind, i) for i in 1:n], ppeak=tr.ppeak,
            gap_prev_ms=ff.gap_prev, gap_next_ms=ff.gap_next, follower=ff.follower,
            closure_phat=C.closure, ncc_min=C.ncc_min, gcc_peak_min=C.gcc_peak_min, gcc_npeaks=C.gcc_npeaks,
            env_peaks=vec(maximum(C.env; dims=2)), impulse5_flagged=w5.shape[g0 .+ (1:n)] .== "triangle",
            nfound=C.nfound)
        for k in 1:K
            ae = [isnan(C.uv[i, k, 1]) ? (NaN, NaN) : uv2azel(C.uv[i, k, 1], C.uv[i, k, 2]) for i in 1:n]
            rows[!, "az_$k"] = rad2deg.(first.(ae)); rows[!, "el_$k"] = rad2deg.(last.(ae))
            rows[!, "score_$k"] = C.score[:, k]; rows[!, "tarr_$k"] = C.tarr[:, k]; rows[!, "track_$k"] = T.track[:, k]
        end
        push!(parts, rows); push!(cands, (; g0, n, filt=T.filt, uv=C.uv))
        g0 += n
    end
    feat = reduce(vcat, parts); tracks = reduce(vcat, trk_parts)
    N = nrow(feat)
    N <= nrow(w6) || error("more clicks ($N) than rows in the Impulse6 file")

    # per-track statistics from the associations
    tracks.colour = track_colours(tracks)
    st = Dict(id => (az=Float64[], el=Float64[], s=Float64[]) for id in tracks.id)
    for k in 1:K, i in 1:N
        id = feat[i, "track_$k"]; id == 0 && continue
        push!(st[id].az, feat[i, "az_$k"]); push!(st[id].el, feat[i, "el_$k"]); push!(st[id].s, feat[i, "score_$k"])
    end
    tv(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
    tracks.az_mean = [median(st[id].az) for id in tracks.id]; tracks.el_mean = [median(st[id].el) for id in tracks.id]
    tracks.jitter_deg = [length(st[id].az) < 2 ? NaN : median([acosd(clamp(sum(tv(st[id].az[q], st[id].el[q]) .* tv(st[id].az[q+1], st[id].el[q+1])), -1, 1)) for q in 1:length(st[id].az)-1]) for id in tracks.id]
    tracks.score_med = [median(st[id].s) for id in tracks.id]
    tracks.duration_s = tracks.t_end .- tracks.t_start
    conf = Dict(zip(tracks.id, tracks.confirmed)); colour = Dict(zip(tracks.id, tracks.colour))
    feat.track_confirmed = [any(k -> (id = feat[i, "track_$k"]; id > 0 && conf[id]), 1:K) for i in 1:N]
    CSV.write(joinpath(out_dir, "data", "click_features_$folder.csv"), feat)
    CSV.write(joinpath(out_dir, "data", "tracks_$folder.csv"), tracks)

    # webui files: impulse rows built here, tonal rows copied from the Impulse6 file
    tonal = w6[N+1:end, :]
    imp6 = w6[1:N, :]
    px_of(u, v) = angle2px(reshape(collect(uv2azel(u, v)), 1, 2), fov_angle)
    r7 = DataFrame(frame=Float64[], px=Float64[], py=Float64[], radius=Float64[], r=Float64[], g=Float64[], b=Float64[], a=Float64[], shape=String[])
    r8 = similar(r7, 0)
    for cd in cands, i in 1:cd.n, k in 1:K
        g = cd.g0 + i
        isnan(cd.uv[i, k, 1]) && continue
        p = px_of(cd.uv[i, k, 1], cd.uv[i, k, 2])
        push!(r8, (imp6.frame[g], p[1], p[2], imp6.radius[g], imp6.r[g], imp6.g[g], imp6.b[g],
            imp6.a[g] * feat[g, "score_$k"] / feat[g, "score_1"], k == 1 ? "circle" : "triangle"))
        id = feat[g, "track_$k"]
        if id > 0 && conf[id]
            q = px_of(cd.filt[i, k, 1], cd.filt[i, k, 2]); col = TRACK_PALETTE[colour[id]]
            push!(r7, (imp6.frame[g], q[1], q[2], imp6.radius[g], 255col[1], 255col[2], 255col[3], imp6.a[g], "circle"))
        end
    end
    for (key, df) in [("Impulse7_srpTracks", r7), ("Impulse8_srpCandidates", r8)]
        df.type = df.shape
        out = vcat(df, tonal[:, names(df)])
        CSV.write(joinpath(viewer_dir, "detection_pixels__color255_webui__$key.csv"), out)
    end

    # summary
    inside(df) = mean((0 .<= df.px .<= imsize[2]) .& (0 .<= df.py .<= imsize[1]))
    flag6 = feat.impulse5_flagged
    match1 = [acosd(clamp(sum(tv(feat.az_1[g], feat.el_1[g]) .* tv(px2angle(imp6.px[g], imp6.py[g])...)), -1, 1)) for g in 1:N]
    ct = tracks[tracks.confirmed, :]
    summ = DataFrame(metric=String[], value=Float64[])
    add(nm, v) = push!(summ, (nm, v))
    add("clicks", N); add("followers_pct", 100mean(feat.follower))
    add("candidate1_vs_Impulse6_within_0.5deg_pct", 100mean(match1 .<= 0.5))
    for k in 1:K; add("clicks_with_ge$(k)_candidates_pct", 100mean(feat.nfound .>= k)); end
    add("score2_over_score1_median", median(filter(!isnan, feat.score_2 ./ feat.score_1)))
    add("tracks_total", nrow(tracks)); add("tracks_confirmed", nrow(ct))
    add("confirmed_duration_s_median", median(ct.duration_s)); add("confirmed_hits_median", median(ct.hits))
    add("confirmed_jitter_deg_median", median(filter(!isnan, ct.jitter_deg)))
    add("clicks_on_confirmed_track_pct", 100mean(feat.track_confirmed))
    add("clicks_on_confirmed_track_pct_isolated", 100mean(feat.track_confirmed[.!feat.follower]))
    add("clicks_on_confirmed_track_pct_followers", 100mean(feat.track_confirmed[feat.follower]))
    add("closure_le_max_pct", 100mean(feat.closure_phat .<= closure_max))
    add("inside_frame_pct_Impulse6", 100inside(imp6)); add("inside_frame_pct_Impulse6g", 100inside(imp6[.!flag6, :]))
    add("inside_frame_pct_Impulse7", 100inside(r7)); add("rows_Impulse7", nrow(r7)); add("rows_Impulse8", nrow(r8))
    CSV.write(joinpath(out_dir, "multiarrival_summary_$folder.csv"), summ)
    open(joinpath(out_dir, "multiarrival_summary_$folder.md"), "w") do io
        println(io, "# Multi-arrival summary: $folder\n\nSource: `$res_folder`\n\n| metric | value |\n|---|---|")
        for r in eachrow(summ); println(io, "| ", r.metric, " | ", @sprintf("%.4g", r.value), " |"); end
    end
    show(stdout, MIME("text/plain"), summ); println()
    return (; feat, tracks, summary=summ)
end

# az, el (deg) of a pixel, inverse of angle2pxd
px2angle(px, py) = (px - imsize[2] / 2) / (imsize[2] / 2) * fov_angle[1] / 2, -(py - imsize[1] / 2) / (imsize[1] / 2) * fov_angle[2] / 2


# ==== round 2: CLEAN candidates, track-level selection ====================================================
# plan: reviews/plans/multiarrival-round2-clean-selection_2026-10-23.md

# With PHAT the GCC curve of one arrival at lag δ is h(τ - δ), h(τ) = mean over the band bins k of cos(2πkτ/N)
# (gcc_curve: N = nextfastfft(2L), bins outside `band` zeroed, peak normalized to 1). Tabulated on a fine lag grid.
function _phat_kernel_table(L, fs, band; span=64.0, step=1 / 64)
    N = nextfastfft(2L)
    f = rfftfreq(N, fs)
    ks = findall((f .>= band[1]) .& (f .<= band[2])) .- 1
    τ = -span:step:span
    return (τ, [mean(k -> cos(2π * k * t / N), ks) for t in τ])
end

# linear interpolation of curve values r on lags at positions x (same formula as curve_at), added into S
function _add_curve!(S, (lags, r), x)
    l1, st, n = first(lags), step(lags), length(r)
    @inbounds for g in eachindex(S)
        y = clamp((x[g] - l1) / st + 1, 1, n)
        k = min(floor(Int, y), n - 1)
        S[g] += r[k] + (y - k) * (r[k+1] - r[k])
    end
    return S
end

# per-pair GCC peak features of round 1 (closure, ncc_min, gcc_peak_min, gcc_npeaks), from the curves of one click
function _gcc_quality(curves, snip, pairs, m)
    τ = zeros(length(pairs)); pk = zeros(length(pairs)); npk = zeros(Int, length(pairs)); ncc = zeros(length(pairs))
    for (p, (j, k)) in enumerate(pairs)
        lags, r = curves[p]
        kf, pk[p] = subsample_peak(r, argmax(r))
        τ[p] = first(lags) + (kf - 1) * step(lags)
        npk[p] = count(q -> r[q] >= r[q-1] && r[q] >= r[q+1] && r[q] >= 0.5pk[p], 2:length(r)-1)
        x, y = snip[:, j], snip[:, k]
        ncc[p] = abs(sum(x .* _frac_delay(y, τ[p]))) / sqrt(max(sum(abs2, x) * sum(abs2, y), floatmin(Float64)))
    end
    return (closure=m == 3 ? abs(τ[1] + τ[3] - τ[2]) : NaN, ncc_min=minimum(ncc), gcc_peak_min=minimum(pk), gcc_npeaks=maximum(npk))
end

"""
    get_srp_candidates_clean(data, peaks; window, ref_channel, fs, rx_vect, c, band=nothing, interp=8, grid_deg=1.0,
        refine=10, K=3, min_sep_deg=8.0, min_rel=0.0, min_pair=0.0, threaded=true)

SRP-PHAT candidates with CLEAN sidelobe suppression. Candidate 1 is the SRP maximum (= `get_tdoa_srp`). After each
accepted candidate, its PHAT kernel (the GCC shape of a single arrival) is subtracted from every pair curve at the
candidate's lag, scaled by the curve value there; the next candidate is the maximum of the residual SRP map at least
`min_sep_deg` away. So sidelobes of an accepted arrival vanish with it, and a second candidate needs evidence of its own.
The search stops when the residual score drops below `min_rel` x candidate 1's score or the weakest pair support is
below `min_pair`. Per click i and candidate k: `uv`, `tdoa`, `tarr` (beam arrival time, samples from the click) as
get_srp_candidates; `score` the SRP of the original curves there; `rscore` the residual SRP when it was picked;
`pmin` the smallest pair-curve value (support) at its lags. Also `nfound` and the round-1 click features
(`closure`, `ncc_min`, `gcc_peak_min`, `gcc_npeaks`). `threaded=false` runs the clicks in sequence (for callers
that parallelize over blocks).
"""
function get_srp_candidates_clean(data, peaks; window=default_window, ref_channel=ref_channel, fs, rx_vect=rx_vect, c=default_c,
    band=nothing, interp=8, grid_deg=1.0, refine=10, K=3, min_sep_deg=8.0, min_rel=0.0, min_pair=0.0, threaded=true)

    m = size(data, 2)
    rx = _rx_for(data, rx_vect)
    τmax = tau_max_samples(rx, fs; c=c)
    band = default_gcc_band(band, fs)
    pairs = [(j, k) for j in 1:m for k in j+1:m]
    np = length(pairs)
    arrival(u, v) = -(rx' * [u, v, sqrt(max(1 - u^2 - v^2, 0.0))]) ./ c .* fs
    pairdelays(u, v) = (t = arrival(u, v); [t[j] - t[k] for (j, k) in pairs])
    h, grid, _, _ = _srp_lattice(grid_deg)
    G = length(grid)
    D = reduce(hcat, [pairdelays(u, v) for (u, v) in grid])
    Dr = [D[p, :] for p in 1:np]
    Ug = reduce(hcat, [collect(_uvw(u, v)) for (u, v) in grid])
    fine = range(-h, h; length=2refine + 1)
    kernels = Dict{Int,Any}()
    klock = ReentrantLock()
    kernel_for(L) = lock(() -> get!(() -> _phat_kernel_table(L, fs, band), kernels, L), klock)

    n = length(peaks)
    uv = fill(NaN, n, K, 2); score = fill(NaN, n, K); rscore = fill(NaN, n, K); pmin = fill(NaN, n, K)
    tdoa = fill(NaN, n, K, m); tarr = fill(NaN, n, K); nfound = zeros(Int, n)
    closure = fill(NaN, n); ncc_min = fill(NaN, n); gcc_peak_min = fill(NaN, n); gcc_npeaks = zeros(Int, n)
    function one(i)
        win = _clip_window(window .+ peaks[i], size(data, 1))
        snip = Float64.(data[win, :])
        curves = [gcc_curve(snip[:, j], snip[:, k]; fs=fs, maxlag=τmax[j, k], weighting=:phat, band=band, interp=interp) for (j, k) in pairs]
        q = _gcc_quality(curves, snip, pairs, m)
        closure[i], ncc_min[i], gcc_peak_min[i], gcc_npeaks[i] = q.closure, q.ncc_min, q.gcc_peak_min, q.gcc_npeaks
        kτ, kv = kernel_for(length(win))
        kern(t) = (y = clamp((t - first(kτ)) / step(kτ) + 1, 1, length(kv)); k = min(floor(Int, y), length(kv) - 1); kv[k] + (y - k) * (kv[k+1] - kv[k]))
        res = [(cv[1], copy(cv[2])) for cv in curves]
        sc(cvs, d) = sum(p -> curve_at(cvs[p], d[p]), 1:np)
        S = zeros(G); s1 = NaN; acc = Vector{Vector{Float64}}()
        for kk in 1:K
            fill!(S, 0.0)
            for p in 1:np; _add_curve!(S, res[p], Dr[p]); end
            for a in acc
                @inbounds for g in 1:G
                    acosd(clamp(a[1] * Ug[1, g] + a[2] * Ug[2, g] + a[3] * Ug[3, g], -1, 1)) < min_sep_deg && (S[g] = -Inf)
                end
            end
            g = argmax(S)
            isfinite(S[g]) || break
            best, (bu, bv) = S[g], grid[g]
            u0, v0 = bu, bv
            for du in fine, dv in fine
                u, v = u0 + du, v0 + dv
                u^2 + v^2 <= 1 || continue
                kk > 1 && any(a -> acosd(clamp(sum(a .* _uvw(u, v)), -1, 1)) < min_sep_deg, acc) && continue
                s = sc(res, pairdelays(u, v))
                s > best && ((best, bu, bv) = (s, u, v))
            end
            d = pairdelays(bu, bv)
            amp = [max(curve_at(res[p], d[p]), 0.0) for p in 1:np]
            rs = best / np
            if kk == 1
                s1 = rs
            elseif rs < min_rel * s1 || minimum(amp) < min_pair
                break
            end
            t = arrival(bu, bv)
            uv[i, kk, :] = [bu, bv]; rscore[i, kk] = rs; score[i, kk] = sc(curves, d) / np; pmin[i, kk] = minimum(amp)
            tdoa[i, kk, :] = t .- t[ref_channel]
            beam = sum(ch -> _frac_delay(snip[:, ch], -(t[ch] - t[ref_channel])), 1:m)
            tarr[i, kk] = win[argmax(abs.(hilbert(beam)))] - peaks[i]
            for p in 1:np
                lags, r = res[p]
                for (q_, l) in enumerate(lags); r[q_] -= amp[p] * kern(l - d[p]); end
            end
            push!(acc, collect(_uvw(bu, bv)))
            nfound[i] = kk
        end
    end
    if threaded
        Threads.@threads for i in 1:n; one(i); end
    else
        for i in 1:n; one(i); end
    end
    return (; uv, score, rscore, pmin, tdoa, tarr, nfound, closure, ncc_min, gcc_peak_min, gcc_npeaks)
end

"""
    track_lags(track, tarr, confirmed) -> Dict(id => (lag, n))

For each confirmed track: the median of (its candidate's beam arrival time - the earliest arrival among the confirmed-
track candidates of the same click), over the `n` clicks where it shares the click with another confirmed track.
A direct path arrives first (lag near 0); an image path lags. `confirmed` maps track id => Bool.
"""
function track_lags(track, tarr, confirmed)
    acc = Dict{Int,Vector{Float64}}()
    for i in axes(track, 1)
        ks = [k for k in axes(track, 2) if track[i, k] > 0 && confirmed[track[i, k]] && !isnan(tarr[i, k])]
        length(ks) < 2 && continue
        t0 = minimum(tarr[i, ks])
        for k in ks; push!(get!(acc, track[i, k], Float64[]), tarr[i, k] - t0); end
    end
    return Dict(id => (median(v), length(v)) for (id, v) in acc)
end

"""
    select_candidates(track, tarr, score, confirmed, hits; lags=track_lags(...), lag_echo=3.0, min_shared=5)

One candidate per click. Among the candidates on confirmed tracks, take the one whose track has the smallest
`track_lags` lag (tracks sharing fewer than `min_shared` clicks count as lag 0: no evidence), then more hits, then
a higher score. Returns (; sel, status): `sel[i]` the candidate index (rank 1 when no candidate is on a confirmed
track), `status[i]` "direct" (chosen lag <= `lag_echo` samples), "echo" (only lagging tracks) or "untracked".
"""
function select_candidates(track, tarr, score, confirmed, hits; lags=track_lags(track, tarr, confirmed), lag_echo=3.0, min_shared=5)
    n, K = size(track)
    lagof(id) = (l = get(lags, id, (0.0, 0)); l[2] >= min_shared ? l[1] : 0.0)
    sel = ones(Int, n); status = fill("untracked", n)
    for i in 1:n
        ks = [k for k in 1:K if track[i, k] > 0 && confirmed[track[i, k]]]
        isempty(ks) && continue
        best = argmin(k -> (lagof(track[i, k]), -hits[track[i, k]], -score[i, k]), ks)
        sel[i] = best
        status[i] = lagof(track[i, best]) <= lag_echo ? "direct" : "echo"
    end
    return (; sel, status)
end

"""
    run_multiarrival2(res_folder; out_dir, viewer_dir=res_folder, K=3, block_s=10.0, pad_s=0.05, echo_window_ms=1.6,
        candidate_kwargs=(;), tracker_kwargs=(;), select_kwargs=(;), tag="clean")

Round-2 version of `run_multiarrival`: CLEAN candidates (get_srp_candidates_clean), parallel over audio blocks,
tracks labelled by their arrival lag (track_lags) and one selected candidate per click (select_candidates).
Writes `out_dir/data/click_features_<tag>_<folder>.csv`, `out_dir/data/tracks_<tag>_<folder>.csv`,
`out_dir/multiarrival_summary_<tag>_<folder>.md/.csv`, and in `viewer_dir` the webui files Impulse7c_srpCleanTracks,
Impulse8c_srpCleanCandidates and Impulse9_selected (selected direction, click colours; echo / untracked as triangles).
"""
function run_multiarrival2(res_folder; out_dir, viewer_dir=res_folder, K=3, block_s=10.0, pad_s=0.05, echo_window_ms=1.6,
    candidate_kwargs=(;), tracker_kwargs=(;), select_kwargs=(;), tag="clean")

    folder = basename(rstrip(res_folder, '/'))
    mkpath(joinpath(out_dir, "data"))
    files = joinpath.(res_folder, sort(filter(f -> occursin("__cps", f) && endswith(f, ".jld2") && !endswith(f, "_angles.jld2"), readdir(res_folder))))
    w6 = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__Impulse6_srpPhat.csv"), DataFrame)
    w5 = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__Impulse5_gccPhatGuided.csv"), DataFrame)
    rx = rx_vect[:, get_relevant_channels(rx_vect)]
    parts = DataFrame[]; trk_parts = DataFrame[]; cands = []; g0 = 0; next_id = 1
    for (fi, f) in enumerate(files)
        j = jldopen(f); tr = j["res_impulsetrain"]; fs = j["fs"]; ri = j["res_impulse"]; close(j)
        aufname, len_data = ri.aufname, ri.len_data
        P = tr.pind_good; n = length(P)
        @info "multiarrival2: $(basename(aufname)), $n clicks"
        st = [tr.train_start_ind; n + 1]; blocks = UnitRange{Int}[]
        for t in 1:length(tr.train_start_ind)
            a = st[t]
            while a < st[t+1]
                b = a
                while b + 1 < st[t+1] && P[b+1] - P[a] <= block_s * fs; b += 1; end
                push!(blocks, a:b); a = b + 1
            end
        end
        C = (; uv=fill(NaN, n, K, 2), score=fill(NaN, n, K), rscore=fill(NaN, n, K), pmin=fill(NaN, n, K), tarr=fill(NaN, n, K),
            nfound=zeros(Int, n), closure=fill(NaN, n), ncc_min=fill(NaN, n), gcc_peak_min=fill(NaN, n), gcc_npeaks=zeros(Int, n), env=zeros(Int, n))
        pad = round(Int, pad_s * fs)
        Threads.@threads :dynamic for bl in blocks[sortperm(length.(blocks); rev=true)]   # longest first
            lo = max(P[bl[1]] - pad, 1)
            d = read_filtered_block(aufname, lo:min(P[bl[end]] + pad, len_data); fs=fs)
            pk = P[bl] .- lo .+ 1
            c = get_srp_candidates_clean(d, pk; window=window_impulsive, ref_channel=ref_channel, fs=fs, rx_vect=rx, K=K, threaded=false, candidate_kwargs...)
            C.uv[bl, :, :] = c.uv
            for nm in (:score, :rscore, :pmin, :tarr); getfield(C, nm)[bl, :] = getfield(c, nm); end
            for nm in (:nfound, :closure, :ncc_min, :gcc_peak_min, :gcc_npeaks); getfield(C, nm)[bl] = getfield(c, nm); end
            ea = zeros(Int, length(bl))
            for (q, p) in enumerate(pk)   # envelope_arrivals without its inner threads
                win = _clip_window(window_impulsive .+ p, size(d, 1))
                ea[q] = maximum(ch -> (pkk = findmaxima(abs.(hilbert(d[win, ch])), 8); isempty(pkk.heights) ? 0 : count(pkk.heights .>= 0.3maximum(pkk.heights))), axes(d, 2))
            end
            C.env[bl] = ea
        end
        ff = follower_flags(P, fs; echo_window_ms=echo_window_ms)
        tsec = (P .- 1) ./ fs
        T = track_candidates(tsec, C.uv, C.score; first_id=next_id, tracker_kwargs...)
        next_id = isempty(T.tracks.id) ? next_id : maximum(T.tracks.id) + 1
        conf = Dict(zip(T.tracks.id, T.tracks.confirmed)); hits = Dict(zip(T.tracks.id, T.tracks.hits))
        lg = track_lags(T.track, C.tarr, conf)
        S = select_candidates(T.track, C.tarr, C.score, conf, hits; lags=lg, select_kwargs...)
        T.tracks.file .= fi
        T.tracks.lag = [get(lg, id, (NaN, 0))[1] for id in T.tracks.id]
        T.tracks.n_shared = [get(lg, id, (NaN, 0))[2] for id in T.tracks.id]
        push!(trk_parts, T.tracks)
        rows = DataFrame(file=basename(aufname), click=1:n, sample=P, t_s=tsec, frame=w6.frame[g0 .+ (1:n)],
            train=[searchsortedlast(tr.train_start_ind, i) for i in 1:n], ppeak=tr.ppeak,
            gap_prev_ms=ff.gap_prev, gap_next_ms=ff.gap_next, follower=ff.follower,
            closure_phat=C.closure, ncc_min=C.ncc_min, gcc_peak_min=C.gcc_peak_min, gcc_npeaks=C.gcc_npeaks,
            env_peaks=C.env, impulse5_flagged=w5.shape[g0 .+ (1:n)] .== "triangle", nfound=C.nfound,
            selected=S.sel, select_status=S.status)
        for k in 1:K
            ae = [isnan(C.uv[i, k, 1]) ? (NaN, NaN) : uv2azel(C.uv[i, k, 1], C.uv[i, k, 2]) for i in 1:n]
            rows[!, "az_$k"] = rad2deg.(first.(ae)); rows[!, "el_$k"] = rad2deg.(last.(ae))
            rows[!, "score_$k"] = C.score[:, k]; rows[!, "rscore_$k"] = C.rscore[:, k]; rows[!, "pmin_$k"] = C.pmin[:, k]
            rows[!, "tarr_$k"] = C.tarr[:, k]; rows[!, "track_$k"] = T.track[:, k]
        end
        ae = [uv2azel(C.uv[i, S.sel[i], 1], C.uv[i, S.sel[i], 2]) for i in 1:n]
        rows.az_sel = rad2deg.(first.(ae)); rows.el_sel = rad2deg.(last.(ae))
        push!(parts, rows); push!(cands, (; g0, n, filt=T.filt, uv=C.uv, sel=S.sel))
        g0 += n
    end
    feat = reduce(vcat, parts); tracks = reduce(vcat, trk_parts)
    N = nrow(feat)

    tracks.colour = track_colours(tracks)
    st = Dict(id => (az=Float64[], el=Float64[], s=Float64[]) for id in tracks.id)
    for k in 1:K, i in 1:N
        id = feat[i, "track_$k"]; id == 0 && continue
        push!(st[id].az, feat[i, "az_$k"]); push!(st[id].el, feat[i, "el_$k"]); push!(st[id].s, feat[i, "score_$k"])
    end
    tv(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
    angd2(a, b) = acosd(clamp(sum(a .* b), -1, 1))
    tracks.az_mean = [median(st[id].az) for id in tracks.id]; tracks.el_mean = [median(st[id].el) for id in tracks.id]
    tracks.jitter_deg = [length(st[id].az) < 2 ? NaN : median([angd2(tv(st[id].az[q], st[id].el[q]), tv(st[id].az[q+1], st[id].el[q+1])) for q in 1:length(st[id].az)-1]) for id in tracks.id]
    tracks.score_med = [median(st[id].s) for id in tracks.id]
    tracks.duration_s = tracks.t_end .- tracks.t_start
    conf = Dict(zip(tracks.id, tracks.confirmed)); colour = Dict(zip(tracks.id, tracks.colour))
    feat.track_confirmed = [any(k -> (id = feat[i, "track_$k"]; id > 0 && conf[id]), 1:K) for i in 1:N]
    CSV.write(joinpath(out_dir, "data", "click_features_$(tag)_$folder.csv"), feat)
    CSV.write(joinpath(out_dir, "data", "tracks_$(tag)_$folder.csv"), tracks)

    tonal = w6[N+1:end, :]; imp6 = w6[1:N, :]
    px_of(u, v) = angle2px(reshape(collect(uv2azel(u, v)), 1, 2), fov_angle)
    mk() = DataFrame(frame=Float64[], px=Float64[], py=Float64[], radius=Float64[], r=Float64[], g=Float64[], b=Float64[], a=Float64[], shape=String[])
    r7, r8, r9 = mk(), mk(), mk()
    for cd in cands, i in 1:cd.n
        g = cd.g0 + i
        for k in 1:K
            isnan(cd.uv[i, k, 1]) && continue
            p = px_of(cd.uv[i, k, 1], cd.uv[i, k, 2])
            push!(r8, (imp6.frame[g], p[1], p[2], imp6.radius[g], imp6.r[g], imp6.g[g], imp6.b[g],
                imp6.a[g] * clamp(feat[g, "rscore_$k"] / feat[g, "rscore_1"], 0, 1), k == 1 ? "circle" : "triangle"))
            id = feat[g, "track_$k"]
            if id > 0 && conf[id]
                q = px_of(cd.filt[i, k, 1], cd.filt[i, k, 2]); col = TRACK_PALETTE[colour[id]]
                push!(r7, (imp6.frame[g], q[1], q[2], imp6.radius[g], 255col[1], 255col[2], 255col[3], imp6.a[g], "circle"))
            end
        end
        k = cd.sel[i]; id = feat[g, "track_$k"]
        q = id > 0 && conf[id] ? px_of(cd.filt[i, k, 1], cd.filt[i, k, 2]) : px_of(cd.uv[i, k, 1], cd.uv[i, k, 2])
        push!(r9, (imp6.frame[g], q[1], q[2], imp6.radius[g], imp6.r[g], imp6.g[g], imp6.b[g], imp6.a[g],
            feat.select_status[g] == "direct" ? "circle" : "triangle"))
    end
    for (key, df) in [("Impulse7c_srpCleanTracks", r7), ("Impulse8c_srpCleanCandidates", r8), ("Impulse9_selected", r9)]
        df.type = df.shape
        CSV.write(joinpath(viewer_dir, "detection_pixels__color255_webui__$key.csv"), vcat(df, tonal[:, names(df)]))
    end

    # summary, including the consecutive-click agreement of Impulse6 (candidate 1) vs the selection
    inside(df) = mean((0 .<= df.px .<= imsize[2]) .& (0 .<= df.py .<= imsize[1]))
    ct = tracks[tracks.confirmed, :]
    function agreement(az, el, lo, hi)
        a = Float64[]
        for i in 2:N
            feat.file[i] == feat.file[i-1] && lo <= feat.gap_prev_ms[i] < hi && push!(a, angd2(tv(az[i], el[i]), tv(az[i-1], el[i-1])))
        end
        return median(a), 100mean(a .< 2)
    end
    summ = DataFrame(metric=String[], value=Float64[])
    add(nm, v) = push!(summ, (nm, v))
    add("clicks", N)
    match1 = [angd2(tv(feat.az_1[g], feat.el_1[g]), tv(px2angle(imp6.px[g], imp6.py[g])...)) for g in 1:N]
    add("candidate1_vs_Impulse6_within_0.5deg_pct", 100mean(match1 .<= 0.5))
    for k in 1:K; add("clicks_with_ge$(k)_candidates_pct", 100mean(feat.nfound .>= k)); end
    add("rscore2_over_score1_median", median(filter(!isnan, feat.rscore_2 ./ feat.score_1)))
    add("tracks_total", nrow(tracks)); add("tracks_confirmed", nrow(ct))
    add("confirmed_duration_s_median", median(ct.duration_s)); add("confirmed_hits_median", median(ct.hits))
    add("confirmed_jitter_deg_median", median(filter(!isnan, ct.jitter_deg)))
    add("confirmed_tracks_with_lag_evidence", count(ct.n_shared .>= 5)); add("confirmed_tracks_lagging_gt3", count((ct.n_shared .>= 5) .& (ct.lag .> 3)))
    add("clicks_on_confirmed_track_pct", 100mean(feat.track_confirmed))
    for s in ("direct", "echo", "untracked"); add("selected_$(s)_pct", 100mean(feat.select_status .== s)); end
    add("selected_not_rank1_pct", 100mean(feat.selected .!= 1))
    for (lo, hi) in [(0.4, 1.6), (1.6, 2.2), (2.2, 10.0), (10.0, 30.0)]
        m6, s6 = agreement(feat.az_1, feat.el_1, lo, hi); m9, s9 = agreement(feat.az_sel, feat.el_sel, lo, hi)
        add("consecutive_$(lo)-$(hi)ms_median_deg_Impulse6", m6); add("consecutive_$(lo)-$(hi)ms_median_deg_Impulse9", m9)
        add("consecutive_$(lo)-$(hi)ms_lt2deg_pct_Impulse6", s6); add("consecutive_$(lo)-$(hi)ms_lt2deg_pct_Impulse9", s9)
    end
    add("inside_frame_pct_Impulse6", 100inside(imp6)); add("inside_frame_pct_Impulse9", 100inside(r9))
    add("inside_frame_pct_Impulse9_direct", 100inside(r9[r9.shape .== "circle", :]))
    add("rows_Impulse7c", nrow(r7)); add("rows_Impulse8c", nrow(r8)); add("rows_Impulse9", nrow(r9))
    CSV.write(joinpath(out_dir, "multiarrival_summary_$(tag)_$folder.csv"), summ)
    open(joinpath(out_dir, "multiarrival_summary_$(tag)_$folder.md"), "w") do io
        println(io, "# Multi-arrival summary ($tag): $folder\n\nSource: `$res_folder`\n\n| metric | value |\n|---|---|")
        for r in eachrow(summ); println(io, "| ", r.metric, " | ", @sprintf("%.4g", r.value), " |"); end
    end
    show(stdout, MIME("text/plain"), summ); println()
    return (; feat, tracks, summary=summ)
end

# support of given directions on the plain GCC-PHAT curves: per click and candidate, the smallest pair-curve value at
# the direction's lags (a true arrival is supported by every pair, a sidelobe ridge by one); NaN where uv is NaN
function srp_support(data, peaks, uv; window=default_window, fs, rx_vect=rx_vect, c=default_c, band=nothing, interp=8)
    m = size(data, 2); rx = _rx_for(data, rx_vect); τmax = tau_max_samples(rx, fs; c=c); band = default_gcc_band(band, fs)
    pairs = [(j, k) for j in 1:m for k in j+1:m]
    arrival(u, v) = -(rx' * [u, v, sqrt(max(1 - u^2 - v^2, 0.0))]) ./ c .* fs
    out = fill(NaN, size(uv, 1), size(uv, 2))
    Threads.@threads for i in eachindex(peaks)
        win = _clip_window(window .+ peaks[i], size(data, 1)); snip = Float64.(data[win, :])
        curves = [gcc_curve(snip[:, j], snip[:, k]; fs=fs, maxlag=τmax[j, k], weighting=:phat, band=band, interp=interp) for (j, k) in pairs]
        for k in axes(uv, 2)
            isnan(uv[i, k, 1]) && continue
            t = arrival(uv[i, k, 1], uv[i, k, 2])
            out[i, k] = minimum(p -> curve_at(curves[p], t[pairs[p][1]] - t[pairs[p][2]]), eachindex(pairs))
        end
    end
    return out
end

"""
    get_srp_candidates_sic(data, peaks; window, ref_channel, fs, rx_vect, c, band=nothing, interp=8, grid_deg=1.0,
        refine=10, K=3, min_sep_deg=8.0, gate=10, min_rel=0.0, min_pair=-Inf, threaded=true)

SRP-PHAT candidates by successive interference cancellation in the time domain. Candidate 1 is the SRP maximum
(= `get_tdoa_srp`). For each accepted candidate the channels are aligned to it and averaged (beam), the beam is gated
to ±`gate` samples around its envelope peak (Hann taper) to get that arrival's waveform, and the waveform, delayed
by the candidate's TDOAs and scaled by least squares, is subtracted from every channel. GCC-PHAT curves are then
recomputed on the residual channels and the next candidate is the residual SRP maximum, at least `min_sep_deg`
away. Unlike subtracting kernels from the PHAT curves, this stays valid with several arrivals in the window
(PHAT is not additive). Same outputs as `get_srp_candidates_clean`, plus `efrac[i, k]`: share of the window energy
(before that step) removed with candidate k. `rscore` = residual SRP when picked, `score` and `pmin` are evaluated
on the original curves, `tarr` comes from the residual beam (so a weak later arrival is timed without the strong one).
`min_pair` is off by default: a weaker arrival can have negative support on the original PHAT curves, where the
strong arrival dominates.
"""
function get_srp_candidates_sic(data, peaks; window=default_window, ref_channel=ref_channel, fs, rx_vect=rx_vect, c=default_c,
    band=nothing, interp=8, grid_deg=1.0, refine=10, K=3, min_sep_deg=8.0, gate=10, min_rel=0.0, min_pair=-Inf, threaded=true)

    m = size(data, 2)
    rx = _rx_for(data, rx_vect)
    τmax = tau_max_samples(rx, fs; c=c)
    band = default_gcc_band(band, fs)
    pairs = [(j, k) for j in 1:m for k in j+1:m]
    np = length(pairs)
    arrival(u, v) = -(rx' * [u, v, sqrt(max(1 - u^2 - v^2, 0.0))]) ./ c .* fs
    pairdelays(u, v) = (t = arrival(u, v); [t[j] - t[k] for (j, k) in pairs])
    h, grid, _, _ = _srp_lattice(grid_deg)
    G = length(grid)
    D = reduce(hcat, [pairdelays(u, v) for (u, v) in grid])
    Dr = [D[p, :] for p in 1:np]
    Ug = reduce(hcat, [collect(_uvw(u, v)) for (u, v) in grid])
    fine = range(-h, h; length=2refine + 1)
    gccs(x) = [gcc_curve(x[:, j], x[:, k]; fs=fs, maxlag=τmax[j, k], weighting=:phat, band=band, interp=interp) for (j, k) in pairs]
    sc(cvs, d) = sum(p -> curve_at(cvs[p], d[p]), 1:np)
    far(acc, w) = all(a -> acosd(clamp(sum(a .* w), -1, 1)) >= min_sep_deg, acc)

    n = length(peaks)
    uv = fill(NaN, n, K, 2); score = fill(NaN, n, K); rscore = fill(NaN, n, K); pmin = fill(NaN, n, K); efrac = fill(NaN, n, K)
    tdoa = fill(NaN, n, K, m); tarr = fill(NaN, n, K); nfound = zeros(Int, n)
    closure = fill(NaN, n); ncc_min = fill(NaN, n); gcc_peak_min = fill(NaN, n); gcc_npeaks = zeros(Int, n)
    function one(i)
        win = _clip_window(window .+ peaks[i], size(data, 1)); L = length(win)
        snip = Float64.(data[win, :]); resid = copy(snip)
        curves0 = gccs(snip); curves = curves0
        q = _gcc_quality(curves0, snip, pairs, m)
        closure[i], ncc_min[i], gcc_peak_min[i], gcc_npeaks[i] = q.closure, q.ncc_min, q.gcc_peak_min, q.gcc_npeaks
        S = zeros(G); s1 = NaN; acc = Vector{NTuple{3,Float64}}()
        for kk in 1:K
            fill!(S, 0.0)
            for p in 1:np; _add_curve!(S, curves[p], Dr[p]); end
            if !isempty(acc)
                @inbounds for g in 1:G
                    far(acc, (Ug[1, g], Ug[2, g], Ug[3, g])) || (S[g] = -Inf)
                end
            end
            g = argmax(S)
            isfinite(S[g]) || break
            best, (bu, bv) = S[g], grid[g]
            u0, v0 = bu, bv
            for du in fine, dv in fine
                u, v = u0 + du, v0 + dv
                u^2 + v^2 <= 1 || continue
                isempty(acc) || far(acc, _uvw(u, v)) || continue
                s = sc(curves, pairdelays(u, v))
                s > best && ((best, bu, bv) = (s, u, v))
            end
            d = pairdelays(bu, bv)
            rs = best / np
            supp = minimum(p -> curve_at(curves0[p], d[p]), 1:np)
            if kk == 1
                s1 = rs
            elseif rs < min_rel * s1 || supp < min_pair
                break
            end
            t = arrival(bu, bv); τ = t .- t[ref_channel]
            beam = sum(ch -> _frac_delay(resid[:, ch], -τ[ch]), 1:m) ./ m
            pk = argmax(abs.(hilbert(beam)))
            tmpl = zeros(L)
            for s in max(1, pk - gate):min(L, pk + gate); tmpl[s] = beam[s] * 0.5 * (1 + cos(π * (s - pk) / (gate + 1))); end
            E = sum(abs2, resid)
            for ch in 1:m
                model = _frac_delay(tmpl, τ[ch])
                resid[:, ch] .-= (dot(resid[:, ch], model) / max(dot(model, model), floatmin(Float64))) .* model
            end
            uv[i, kk, :] = [bu, bv]; rscore[i, kk] = rs; score[i, kk] = sc(curves0, d) / np; pmin[i, kk] = supp
            efrac[i, kk] = 1 - sum(abs2, resid) / max(E, floatmin(Float64))
            tdoa[i, kk, :] = τ; tarr[i, kk] = win[pk] - peaks[i]
            push!(acc, _uvw(bu, bv)); nfound[i] = kk
            kk < K && (curves = gccs(resid))
        end
    end
    if threaded
        Threads.@threads for i in 1:n; one(i); end
    else
        for i in 1:n; one(i); end
    end
    return (; uv, score, rscore, pmin, efrac, tdoa, tarr, nfound, closure, ncc_min, gcc_peak_min, gcc_npeaks)
end
