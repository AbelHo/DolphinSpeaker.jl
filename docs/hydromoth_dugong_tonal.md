# HydroMoth dugong chirp / trill tonal detector preset

`set_device__hydromoth_dugong` configures candidate tonal-frame detection
for **chirps and trills** in HydroMoth audio. Whistles are **not a target**
of this revision. Frequency filtering alone cannot exclude all whistles
because their frequency ranges overlap; this is not a call-type classifier.

This is a follow-up to the broad 800–16,000 Hz / 20 ms / +5 dB preset.
The revised default focuses on the annotated chirp/trill frequency region,
uses 10 ms frames for short chirps, and reinstates a correctly scaled
frequency-width filter to reject broad-band background frames.

## Usage

```julia
using DolphinSpeaker

# A single pass for chirps and trills (balanced false-alarm control).
set_device__hydromoth_dugong()

# To favor short chirps:
set_device__hydromoth_dugong(call_type=:chirp)

# To favor the lower-frequency, longer trills:
set_device__hydromoth_dugong(call_type=:trill)

# If recall is more important than false alarms:
set_device__hydromoth_dugong(call_type=:trill, sensitivity=:high_recall)

# Then call the existing tonal detection pipeline:
# DolphinSpeaker.detect_tonal("/path/to/recording.WAV", "/path/to/results")
```

Configuration globals are applied to subsequent detector calls. An explicit
keyword parameter in `detect_tonal` overrides the corresponding global.
The original `set_device__hydromoth()` remains unchanged.

## Proposed parameter sets

| Parameter | Chirp + trill (default) | Chirp | Trill |
|---|---:|---:|---:|
| `tonal_band_pass` (Hz) | 1,200–10,500 | 1,200–10,500 | 1,200–6,500 |
| `nfft_inS` | 0.010 s | 0.010 s | 0.020 s |
| FFT frequency resolution at 48 kHz | ~100 Hz | ~100 Hz | ~50 Hz |
| `tonal_threshold_median` balanced | +7 dB | +7 dB | +7 dB |
| `freq_width_db` | 6 dB | 6 dB | 6 dB |
| `freq_maxbandwidth` balanced | 3,000 Hz | 2,500 Hz | 2,000 Hz |
| `percent_quiet` | 0.05 | 0.05 | 0.05 |
| `threshold_tonal` | `nothing` | `nothing` | `nothing` |

Set `sensitivity=:high_recall` in any profile to use +4 dB above the
median frame score and `freq_maxbandwidth=Inf` (bandwidth veto disabled).
This will usually produce more candidates, including more false positives.
The default is `:balanced`.

## Evidence from manually drawn Raven selections

The three Raven tables contain 957 annotated chirps and 311 trills,
for 1,268 primary calls of the two current target types.

| Annotation statistic | Chirp (n=957) | Trill (n=311) |
|---|---:|---:|
| Duration median | 0.118 s | 0.638 s |
| Duration p05–p95 | 0.069–0.196 s | 0.164–2.770 s |
| Raven peak frequency median | 3.938 kHz | 3.375 kHz |
| Peak frequency p05–p95 | 2.344–8.531 kHz | 2.156–5.438 kHz |
| Drawn lower bound p05 | 1.269 kHz | 1.322 kHz |
| Drawn upper bound p95 | 9.402 kHz | 6.087 kHz |
| Drawn frequency width median | 2.108 kHz | 1.601 kHz |

Values are **Raven box metadata**, not measured instantaneous spectral
bandwidth, acoustic SNR, intra-trill pulse rates or validation results.
Some chirps are faint and trills can be shorter than 300 ms. Hard event-
duration gates would miss some genuine annotations. Long-duration trills
may require combining successive tonal frames downstream.

## Detector fixes in this branch

The old per-frame bandwidth test multiplied the number of bins by
`freqss[2]` (the absolute frequency of the second bin), which is not
the FFT-bin spacing. It also mixed a vector score with a transposed
bandwidth mask, potentially creating unintended broadcasting, and
unconditionally dropped the final detected frame. These are corrected.
A short-recording edge case in the quiet-frame baseline now uses at least
one frame instead of the entire recording. The frequency-bin indexing is
one-based and clipped at Nyquist. Returned sample indices are also
one-based for downstream array indexing.

Even with these mechanical fixes, this is an **unvalidated candidate
detector preset**. The provided Raven labels contain no matching raw WAV
samples, so no true precision, recall, false-positive rate, or optimal
threshold can be inferred. Verify on held-out labeled HydroMoth WAVs,
especially faint chirps, overlapping chirp/trill morphologies and
background sounds. Compare per-event recall and false positives/hour.

For the default 10,500 Hz upper edge, use at least a 24 kHz sampling
rate (48 kHz is a useful target). A lower sample rate causes the detector
to clip the band at Nyquist; a lower band edge above Nyquist raises an
error. `window_tonal` is set for consistency but `nfft_inS` actually
controls the spectrogram window in the current implementation.

## Tests

```bash
julia --project=. test/test_hydromoth_dugong_preset.jl
```

The test covers profile switching, invalid options, one-channel synthetic
chirp and trill signals, score/mask shapes and zero-detection behavior.
Synthetic tests are only software regression checks, not field validation.
