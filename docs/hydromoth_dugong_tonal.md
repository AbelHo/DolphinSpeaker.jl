# HydroMoth dugong tonal detector preset

The new `set_device__hydromoth_dugong(; sensitivity=:high_recall)` preset
is defined in `src/config.jl` and exported from `DolphinSpeaker`.
It calls the existing `set_device__hydromoth()` first to retain the
HydroMoth filename timestamp parser, reference channel and non-tonal
configurations. It does **not** alter `set_device__hydromoth()`.

## Usage

```julia
using DolphinSpeaker

set_device__hydromoth_dugong()                       # default high recall
# set_device__hydromoth_dugong(sensitivity=:balanced) # more selective

# Start the existing detection pipeline after applying this preset:
# DolphinSpeaker.detect_tonal("/path/to/20260408_100000.WAV", "/path/to/output")
```

**Required sampling rate:** more than 32 kHz to cover the upper search
boundary of 16 kHz. Prefer 48 kHz or higher. Confirm the actual WAV
header before running. Configuration globals affect calls made *after*
the setter; any explicit keyword arguments in detector calls take
precedence.

## Preset values

| Setting | High recall (default) | Balanced | Why |
|---|---:|---:|---|
| `tonal_band_pass` | [800, 16000] Hz | same | Includes most chirp, trill and whistle frequency regions with margins |
| `nfft_inS` | 0.02 s | same | ~50 Hz Fourier-bin spacing; frames short enough for most labeled chirps |
| `threshold_tonal` | `nothing` | same | Use existing automatic relative threshold |
| `tonal_threshold_median` | +5 dB | +8 dB | Lower threshold improves candidate recall at the cost of more false alarms |
| `percent_quiet` | 0.01 | same | Existing 1% low-energy baseline calculation |
| `freq_width_db` | 6 (veto off) | 3 | Existing per-frame narrowband comparison |
| `freq_maxbandwidth` | `Inf` | 15,000 | Recall mode bypasses the current narrowband veto; balanced mode enables it |

`window_tonal` is set to 0.02 s to match the FFT window but the current
`detect_tonal` implementation does not actually use this keyword to
compute the spectrogram; `nfft_inS` is the operative parameter.
`tonal_segment_threshold` and `tonal_segment_thresh_continue` remain
unchanged, because annotation duration distributions cannot calibrate
their Gaussian-convolution units.

**Important implementation detail:** the existing
`detect_tonal` code's narrowband test uses the second frequency-bin
**value** rather than the **frequency-bin spacing** when scaling
the number of bins in its width estimate. Consequently,
`freq_maxbandwidth` should not be interpreted as a properly calibrated
physical bandwidth without correcting that implementation. The default
high-recall mode bypasses this filter rather than applying a possibly
misleading box-width cutoff. The balanced mode is a conservative
heuristic, not a calibrated bandwidth filter.

## Basis in the Raven annotation data

Three Sibu Tengah Raven selection tables contain 1,743 selections,
including **1,329 primary chirp/trill/whistle boxes**:

| Type | Boxes | Duration p05–p95 | Median peak frequency |
|---|---:|---:|---:|
| Chirp | 957 | ~0.069–0.196 s | ~3.94 kHz |
| Trill | 311 | ~0.164–2.770 s | ~3.38 kHz |
| Whistle | 61 | ~0.072–1.108 s | ~3.28 kHz |

These ranges explain the lower search boundary, the need to avoid
duration-based hard cutoffs and the interest in high recall on faint calls.
Some drawn frequency rectangles extend beyond 12 kHz, justifying a
broader candidate-search region than the typical peak values.
The previously supplied JSON acoustic profile contains the complete
quantiles and original notes.

This is **not** a validated dugong-only acoustic classifier:
the Raven shapes are hand-drawn selection rectangles, not instantaneous
frequency contours, and the original WAV audio was not available in
this task. Some whistle annotations are uncertain, harmonic rectangles
may overlap the fundamentals, and acoustic SNR/source identity are
unknown. The output of the existing detector is tonal frames/segments;
chirp/trill/whistle classification and false-positive rejection need
additional validated audio examples. The detector thresholds should be
benchmarked on held-out WAVs (event precision/recall and false positives
per verified recording hour).

## Test

```bash
julia --project=. test/test_hydromoth_dugong_preset.jl
```

This tests both sensitivity modes, invalid-mode rejection, the retained
HydroMoth timestamp function and repeatability. It is an isolated preset
test, not an evaluation on raw audio.
