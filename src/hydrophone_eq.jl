# Equalizing the hydrophones: estimate each hydrophone's response relative to a reference hydrophone from clicks, and
# remove it before TDOA / SRP. Round 4 found that a fixed short response per hydrophone (mount and housing scattering)
# reproduces the real data's incoherence between channels; being fixed relative to the array, it can be estimated.
# plan: reviews/plans/multiarrival-round5-hydrophone-eq_2026-10-24.md
# needs src/localization.jl (get_tdoa_srp, _rx_for, default_gcc_band); works for any number of hydrophones
using FFTW, DSP, Statistics, LinearAlgebra

"""
    estimate_hydrophone_eq(data, peaks; fs, rx_vect=rx_vect, ref=ref_channel, c=default_c, half=150, window=-50:50,
        N=256, band=nothing, iters=3, min_coh=0.3, reg=0.1, magnitude=true, smooth_bins=1, taper=0.25, max_clicks=5000,
        remove_delay=true)

Relative transfer function G_m(f) = H_m(f) / H_ref(f) of each hydrophone (columns of `data`, geometry `rx_vect`) from
clicks at `peaks` (strong single-arrival clicks work best). Each iteration: the direction of every click is estimated by
SRP-PHAT on the snippet equalized with the current estimate, the raw channels are aligned by the plane-wave delays of that
direction, and G_m(f) = Σ A_m conj(A_ref) / Σ |A_ref|² (each click normalized by its energy; ±`smooth_bins` bins
averaged) is accumulated with its coherence. The linear-phase part of G_m - a pure delay, which data alone cannot tell
from a direction bias - is removed and reported in `delay_removed` (samples). The equalizer is E_m = conj(G_m) /
(|G_m|² + reg), normalized; bins outside `band` or with coherence below `min_coh` are left uncorrected (E = 1);
`magnitude=false` corrects the phase only. Returns (; N, fs, f, G, coh, E, ref, band, n_clicks, delay_removed).
"""
function estimate_hydrophone_eq(data, peaks; fs, rx_vect=rx_vect, ref=ref_channel, c=default_c, half=150, window=-50:50,
    N=256, band=nothing, iters=3, min_coh=0.3, reg=0.1, magnitude=true, smooth_bins=1, taper=0.25, max_clicks=5000,
    remove_delay=true)

    m = size(data, 2); rx = _rx_for(data, rx_vect)
    band = default_gcc_band(band, fs)
    f = rfftfreq(N, fs); nb = length(f); inb = (f .>= band[1]) .& (f .<= band[2])
    pk = [p for p in peaks if p - half >= 1 && p + half <= size(data, 1)]
    length(pk) > max_clicks && (pk = pk[round.(Int, range(1, length(pk); length=max_clicks))])
    L = length(window); w = tukey(L, taper)
    arrival(u, v) = -(rx' * [u, v, sqrt(max(1 - u^2 - v^2, 0.0))]) ./ c .* fs
    eq = (; N, fs, f, G=ones(ComplexF64, nb, m), coh=ones(nb, m), E=ones(ComplexF64, nb, m), ref, band, n_clicks=length(pk),
        delay_removed=zeros(m))
    nchunk = max(1, min(4Threads.nthreads(), length(pk)))
    chunks = collect(Iterators.partition(eachindex(pk), cld(length(pk), nchunk)))
    for it in 1:iters
        firs = it == 1 ? nothing : hydrophone_eq_firs(eq)
        acc = [(num=zeros(ComplexF64, nb, m), den=zeros(nb, m), pw=zeros(nb, m)) for _ in chunks]
        Threads.@threads for ci in eachindex(chunks)
            a = acc[ci]
            for i in chunks[ci]
                seg = Float64.(data[pk[i]-half:pk[i]+half, :])
                segq = isnothing(firs) ? seg : apply_hydrophone_eq(seg, eq; firs=firs)
                _, q = get_tdoa_srp(segq, [half + 1]; window=window, ref_channel=ref, fs=fs, rx_vect=rx, c=c, return_quality=true)
                τ = arrival(q.uv[1, 1], q.uv[1, 2]); τ .-= τ[ref]
                A = [rfft([seg[half+1 .+ window, k] .* w; zeros(N - L)]) .* cis.(2π .* f .* τ[k] ./ fs) for k in 1:m]
                e = sum(abs2, A[ref]) + floatmin(Float64)
                for k in 1:m
                    @views a.num[:, k] .+= A[k] .* conj.(A[ref]) ./ e
                    @views a.den[:, k] .+= abs2.(A[ref]) ./ e
                    @views a.pw[:, k] .+= abs2.(A[k]) ./ e
                end
            end
        end
        num = sum(x -> x.num, acc); den = sum(x -> x.den, acc); pw = sum(x -> x.pw, acc)
        if smooth_bins > 0   # average neighbouring bins
            sm(M) = (h = smooth_bins; [mean(M[max(1, i - h):min(nb, i + h), k]) for i in 1:nb, k in 1:m])
            num, den, pw = sm(num), sm(den), sm(pw)
        end
        G = num ./ max.(den, floatmin(Float64))
        coh = abs.(num) ./ sqrt.(max.(pw .* den, floatmin(Float64)))
        E = ones(ComplexF64, nb, m); dl = zeros(m)
        for k in 1:m
            k == ref && (G[:, k] .= 1; coh[:, k] .= 1; continue)
            ok = inb .& (coh[:, k] .>= min_coh)
            count(ok) < 4 && continue
            # remove the linear phase (pure delay), fitted with coherence weights on the unwrapped phase in the band
            ib = findall(inb)
            φ = unwrap(angle.(G[ib, k])); wt = coh[ib, k] .^ 2 .* ok[ib]
            X = [ones(length(ib)) f[ib]]
            ab = (X' * (wt .* X)) \ (X' * (wt .* φ))
            dl[k] = -ab[2] * fs / 2π
            remove_delay && (G[:, k] .*= cis.(-ab[2] .* f))
            Gk = magnitude ? G[:, k] ./ median(abs.(G[ok, k])) : G[:, k] ./ max.(abs.(G[:, k]), floatmin(Float64))
            E[ok, k] .= conj.(Gk[ok]) ./ (abs2.(Gk[ok]) .+ reg) .* (1 + reg)
        end
        eq = (; N, fs, f, G, coh, E, ref, band, n_clicks=length(pk), delay_removed=dl)
    end
    return eq
end

# zero-phase FIR per hydrophone (centred at tap N/2 + 1) from the equalizer on the N-point grid
hydrophone_eq_firs(eq) = [fftshift(irfft(eq.E[:, k], eq.N)) for k in axes(eq.E, 2)]

"""
    apply_hydrophone_eq(data, eq; firs=hydrophone_eq_firs(eq))

Filters each hydrophone (column) with its equalizer FIR (zero-phase, same length out); the reference column is copied.
"""
function apply_hydrophone_eq(data, eq; firs=hydrophone_eq_firs(eq))
    out = Matrix{Float64}(undef, size(data))
    c0 = eq.N ÷ 2 + 1; n = size(data, 1)
    for k in axes(data, 2)
        if k == eq.ref
            out[:, k] .= data[:, k]
        else
            out[:, k] .= conv(Float64.(data[:, k]), firs[k])[c0:c0+n-1]
        end
    end
    return out
end

"""
    estimate_hydrophone_eq_blind(data, peaks; fs, rx_vect=rx_vect, ref=ref_channel, c=default_c, window=-50:50, N=256,
        band=nothing, min_coh=0.3, reg=0.1, magnitude=true, taper=0.25, max_clicks=20000, n_phase=72)

Direction-free estimate of the relative responses G_m = H_m / H_ref, for strong distortion where click directions from
the raw data are unreliable (estimate_hydrophone_eq then aligns clicks wrongly and fails):
- |G_m| from the ratio of the average power spectra (each click normalized by its reference-channel energy);
- the dispersive phase from the averaged increments between adjacent frequency bins of each click's PHAT-normalized
  cross-spectrum Y = X_m conj(X_ref): Y(f + df) conj(Y(f)) ~ G(f + df) conj(G(f)) e^{-i 2π df Δτ}, so the click's own
  delay Δτ only adds a near-linear term; the integrated phase is detrended (removing that term and any pure delay);
- a constant phase per hydrophone (e.g. inverted polarity) by maximizing the summed GCC peaks over `n_phase` angles.
A pure delay per hydrophone cannot be estimated from data alone (it equals a shift of all directions): it is zero
here; set it from calibration data with `with_delays`. `coh` is the coherence of the phase increments.
Returns the same NamedTuple as estimate_hydrophone_eq.
"""
function estimate_hydrophone_eq_blind(data, peaks; fs, rx_vect=rx_vect, ref=ref_channel, c=default_c, window=-50:50, N=256,
    band=nothing, min_coh=0.3, reg=0.1, magnitude=true, taper=0.25, max_clicks=20000, n_phase=72)

    m = size(data, 2); rx = _rx_for(data, rx_vect)
    band = default_gcc_band(band, fs)
    f = rfftfreq(N, fs); nb = length(f); inb = (f .>= band[1]) .& (f .<= band[2]); ib = findall(inb)
    pk = [p for p in peaks if p + first(window) >= 1 && p + last(window) <= size(data, 1)]
    length(pk) > max_clicks && (pk = pk[round.(Int, range(1, length(pk); length=max_clicks))])
    L = length(window); w = tukey(L, taper)
    τmax = tau_max_samples(rx, fs; c=c)
    X = [hcat([rfft([Float64.(data[p .+ window, k]) .* w; zeros(N - L)]) for k in 1:m]...) for p in pk]   # per click: nb x m
    P = zeros(nb, m); D = zeros(ComplexF64, nb - 1, m); cnt = zeros(nb - 1, m)
    for Xc in X
        e = sum(abs2, Xc[:, ref]) + floatmin(Float64)
        P .+= abs2.(Xc) ./ e
        for k in 1:m
            Y = Xc[:, k] .* conj.(Xc[:, ref]); Y ./= abs.(Y) .+ floatmin(Float64)
            D[:, k] .+= Y[2:end] .* conj.(Y[1:end-1]); cnt[:, k] .+= 1
        end
    end
    G = ones(ComplexF64, nb, m); coh = ones(nb, m); E = ones(ComplexF64, nb, m)
    for k in 1:m
        k == ref && continue
        cohk = abs.(D[:, k]) ./ max.(cnt[:, k], 1)
        coh[2:end, k] .= cohk; coh[1, k] = cohk[1]
        # integrate the phase increments across the band, then detrend (click delays and any pure delay)
        φ = zeros(nb)
        for i in ib[2:end]; φ[i] = φ[i-1] + angle(D[i-1, k]); end
        wt = coh[ib, k] .^ 2
        Xd = [ones(length(ib)) f[ib]]; ab = (Xd' * (wt .* Xd)) \ (Xd' * (wt .* φ[ib]))
        φ[ib] .-= Xd * ab
        mag = sqrt.(P[:, k] ./ max.(P[:, ref], floatmin(Float64)))
        G[ib, k] .= mag[ib] .* cis.(φ[ib])
        # constant phase: maximize the summed GCC peaks (physical lags) of the corrected pair over all clicks
        lags = -ceil(Int, τmax[k, ref]):0.25:ceil(Int, τmax[k, ref])
        Ek = zeros(ComplexF64, nb); Ek[ib] .= conj.(G[ib, k]) ./ abs.(G[ib, k])
        best_θ, best = 0.0, -Inf
        for θ in range(0, 2π; length=n_phase + 1)[1:end-1]
            s = 0.0
            for Xc in X[1:min(end, 500)]
                Y = Xc[ib, k] .* conj.(Xc[ib, ref]) .* Ek[ib] .* cis(-θ); Y ./= abs.(Y) .+ floatmin(Float64)
                s += maximum(τ -> real(sum(Y .* cis.(2π .* f[ib] .* τ ./ fs))), lags)
            end
            s > best && ((best_θ, best) = (θ, s))
        end
        G[ib, k] .*= cis(best_θ)
        ok = inb .& (coh[:, k] .>= min_coh)
        count(ok) < 4 && (G[:, k] .= 1; continue)
        Gk = magnitude ? G[:, k] ./ median(abs.(G[ok, k])) : G[:, k] ./ max.(abs.(G[:, k]), floatmin(Float64))
        E[ok, k] .= conj.(Gk[ok]) ./ (abs2.(Gk[ok]) .+ reg) .* (1 + reg)
    end
    return (; N, fs, f, G, coh, E, ref, band, n_clicks=length(pk), delay_removed=fill(NaN, m))
end

"""
    with_delays(eq, delays; reg=0.1)

The equalizer `eq` with a pure delay `delays[m]` (samples, positive = hydrophone m's response is later than the
reference's) put back into each relative response, e.g. per-hydrophone delays from a calibration recording with known
source directions. Bins left uncorrected (E = 1) stay uncorrected.
"""
function with_delays(eq, delays; reg=0.1)
    G = copy(eq.G); E = copy(eq.E)
    for k in axes(G, 2)
        (k == eq.ref || delays[k] == 0) && continue
        G[:, k] .*= cis.(-2π .* eq.f .* delays[k] ./ eq.fs)
        ok = eq.E[:, k] .!= 1
        Gk = G[:, k] ./ median(abs.(G[ok, k]))
        E[ok, k] .= conj.(Gk[ok]) ./ (abs2.(Gk[ok]) .+ reg) .* (1 + reg)
    end
    return merge(eq, (; G, E))
end
