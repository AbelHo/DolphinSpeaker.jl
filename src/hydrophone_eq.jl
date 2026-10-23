# Equalizing the hydrophones: estimate each hydrophone's response relative to a reference hydrophone from clicks, and
# remove it before TDOA / SRP. Round 4 found that a fixed short response per hydrophone (mount and housing scattering)
# reproduces the real data's incoherence between channels; being fixed relative to the array, it can be estimated.
# plan: reviews/plans/multiarrival-round5-hydrophone-eq_2026-10-24.md
# needs src/localization.jl (get_tdoa_srp, _rx_for, default_gcc_band); works for any number of hydrophones
using FFTW, DSP, Statistics, LinearAlgebra

"""
    estimate_hydrophone_eq(data, peaks; fs, rx_vect=rx_vect, ref=ref_channel, c=default_c, half=150, window=-50:50,
        N=256, band=nothing, iters=3, min_coh=0.3, reg=0.1, magnitude=true, smooth_bins=1, taper=0.25, max_clicks=5000)

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
    N=256, band=nothing, iters=3, min_coh=0.3, reg=0.1, magnitude=true, smooth_bins=1, taper=0.25, max_clicks=5000)

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
            G[:, k] .*= cis.(-ab[2] .* f)
            dl[k] = -ab[2] * fs / 2π
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
