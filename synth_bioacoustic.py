"""
Synthetic Bioacoustics Dataset Scaffold
File: synth_bioacoustics.py

Purpose:
  Minimal, single-file Python package to generate synthetic bioacoustic signals (chirps/whistles,
  gaussian pulses/clicks, click trains, harmonic bird-like syllables), mix them into scenes with
  background beds, apply optional room/underwater propagation (pyroomacoustics / bellhop hooks),
  and export labeled clips + a manifest for ML training.

Design goals:
  - Ready-to-run examples and dataset generation CLI function.
  - Optional integrations: scaper, pyroomacoustics, ARLpy/Bellhop (used if installed; otherwise
    there are lightweight fallbacks so you can still generate useful synthetic data).
  - Explicit metadata/label generation for supervised tasks (start/end, class, SNR, params).

Notes:
  - Default sample rate is 192000 Hz to support high-frequency bat/echolocation signals. Change
    via `sr` parameter if you don't need that bandwidth.
  - Dependencies: numpy, scipy, soundfile (pysoundfile). Optional: scaper, pyroomacoustics.

Usage (quick):
  python -c "from synth_bioacoustics import example_generate; example_generate('out', n_clips=50)"

"""

from __future__ import annotations
import os
import json
import csv
import math
import random
from dataclasses import dataclass, asdict
from typing import Optional, List, Dict, Any, Tuple

import numpy as np
from scipy.signal import chirp, gaussian, fftconvolve
import soundfile as sf

# Optional imports (handled gracefully)
try:
    import scaper
    HAS_SCAPER = True
except Exception:
    HAS_SCAPER = False

try:
    import pyroomacoustics as pra
    HAS_PRA = True
except Exception:
    HAS_PRA = False

# ----------------------------- Utilities ---------------------------------

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def db_to_amp(db: float) -> float:
    return 10 ** (db / 20.0)

def normalize_sig(x: np.ndarray, peak: float = 0.99) -> np.ndarray:
    maxv = np.max(np.abs(x))
    if maxv == 0:
        return x
    return (x / maxv) * peak

# ----------------------------- Signal generators -------------------------

def generate_chirp_signal(duration: float,
                          sr: int,
                          f0: float,
                          f1: float,
                          method: str = 'linear',
                          amplitude: float = 1.0,
                          window_ms: float = 1.0) -> np.ndarray:
    """Generate an FM chirp/whistle using scipy.signal.chirp

    duration: seconds
    f0: start frequency (Hz)
    f1: end frequency (Hz)
    method: 'linear'|'quadratic'|'log'
    window_ms: apply a short Hann window ramp to avoid clicks
    """
    t = np.linspace(0, duration, int(round(sr * duration)), endpoint=False)
    y = chirp(t, f0=f0, f1=f1, t1=duration, method=method)

    # optional amplitude envelope (short ramps)
    ramp = int(round(sr * window_ms / 1000.0))
    if ramp > 0 and len(y) > 2 * ramp:
        env = np.ones_like(y)
        env[:ramp] = np.linspace(0, 1, ramp)
        env[-ramp:] = np.linspace(1, 0, ramp)
        y *= env

    y = y * amplitude
    return y.astype(np.float32)


def generate_gausspulse_click(duration: float,
                               sr: int,
                               fc: float,
                               bw: float = 0.5,
                               amplitude: float = 1.0) -> np.ndarray:
    """Generate a gaussian-modulated pulse centered at fc (Hz).

    Note: scipy.signal.gausspulse exists but for explicit control we build a modulated
    Gaussian envelope here.
    """
    n = int(round(sr * duration))
    t = np.linspace(-duration / 2, duration / 2, n)
    sigma = duration * (1 - bw) * 0.25 + (duration * 0.0005)
    env = np.exp(-0.5 * (t / sigma) ** 2)
    carrier = np.cos(2 * np.pi * fc * t)
    y = env * carrier
    y = y * amplitude
    return y.astype(np.float32)


def generate_click_train(sr: int,
                         num_clicks: int = 10,
                         click_duration: float = 0.0005,
                         mean_ici: float = 0.05,
                         ici_jitter: float = 0.01,
                         click_fc: float = 40000.0,
                         amplitude: float = 1.0) -> np.ndarray:
    """Generate a click-train (e.g., odontocete/bat-like) with randomized ICIs.

    Returns a 1D numpy array.
    """
    icis = np.maximum(0.0005, np.random.normal(loc=mean_ici, scale=ici_jitter, size=num_clicks))
    total_dur = click_duration * num_clicks + icis.sum()
    sig = np.zeros(int(math.ceil(total_dur * sr)), dtype=np.float32)

    pos = 0
    for i in range(num_clicks):
        click = generate_gausspulse_click(click_duration, sr, click_fc, amplitude=amplitude)
        idx = int(round(pos * sr))
        end = idx + len(click)
        if end > len(sig):
            # extend
            more = np.zeros(end - len(sig), dtype=np.float32)
            sig = np.concatenate([sig, more])
        sig[idx:idx + len(click)] += click
        pos += click_duration + icis[i]

    return sig


def generate_harmonic_tone(sr: int,
                           duration: float,
                           f0: float,
                           n_harmonics: int = 3,
                           harmonic_decay: float = 0.7,
                           amplitude: float = 1.0) -> np.ndarray:
    """Generate a short harmonic syllable (bird-like)."""
    t = np.linspace(0, duration, int(round(sr * duration)), endpoint=False)
    y = np.zeros_like(t)
    for h in range(1, n_harmonics + 1):
        y += (harmonic_decay ** (h - 1)) * np.sin(2 * np.pi * f0 * h * t)
    # Apply a gentle envelope
    env = np.sin(np.pi * np.linspace(0, 1, len(t)))
    y = y * env * amplitude
    return y.astype(np.float32)

# ----------------------------- Mixing / Scenes ---------------------------

@dataclass
class EventMeta:
    class_label: str
    start_time: float
    end_time: float
    params: Dict[str, Any]


def mix_events_into_background(background: Optional[np.ndarray],
    background_sr: Optional[int],
    events: List[Tuple[np.ndarray, float, str, Dict[str, Any]]],
    clip_duration: float,
    sr: int,
    target_snr_db: float = -5.0) -> Tuple[np.ndarray, List[EventMeta]]:
    """Mix a set of events (signal, start_time, class_label, params) into a background audio.

    If background is None, create silent bed. target_snr_db is event-level SNR relative
    to background RMS; events are mixed one by one using that SNR.

    Returns (mixed_clip, metadata_list)
    """
    n_samples = int(round(clip_duration * sr))
    if background is None:
        bed = np.zeros(n_samples, dtype=np.float32)
        bg_rms = 1e-9
    else:
        # if background sampling rate differs, resample? For simplicity we assume same sr.
        if background_sr != sr:
            raise ValueError('Background SR must equal target SR or pass background as None')
        bed = background[:n_samples].copy().astype(np.float32)
        if len(bed) < n_samples:
            pad = np.zeros(n_samples - len(bed), dtype=np.float32)
            bed = np.concatenate([bed, pad])
        bg_rms = np.sqrt(np.mean(bed ** 2) + 1e-12)

    out = bed.copy()
    meta_list: List[EventMeta] = []

    for sig, start_time, label, params in events:
        start_idx = int(round(start_time * sr))
        end_idx = start_idx + len(sig)
        if end_idx > len(out):
            # truncate the event to fit
            sig = sig[:len(out) - start_idx]
            end_idx = len(out)
        # scale sig to achieve target SNR
        sig_rms = np.sqrt(np.mean(sig ** 2) + 1e-12)
        if bg_rms < 1e-8:
            # if silent bed, scale to absolute amplitude based on target SNR as dBFS
            scale = db_to_amp(target_snr_db)
        else:
            desired_sig_rms = bg_rms * db_to_amp(target_snr_db)
            scale = desired_sig_rms / (sig_rms + 1e-12)
        out[start_idx:end_idx] += sig * scale
        meta_list.append(EventMeta(class_label=label,
                                   start_time=start_time,
                                   end_time=(start_idx + len(sig)) / sr,
                                   params=params))

    out = normalize_sig(out)
    return out, meta_list

# ----------------------------- Propagation / RIR -------------------------

def apply_rir_convolution(sig: np.ndarray, rir: np.ndarray) -> np.ndarray:
    """Convolve signal with an impulse response using FFT convolution.
    RIR should be 1D array.
    """
    convolved = fftconvolve(sig, rir, mode='full')
    # trim or keep full length
    return convolved.astype(np.float32)


def make_simple_underwater_attenuation(sig: np.ndarray, sr: int, range_m: float, alpha_db_per_km: float = 10.0) -> np.ndarray:
    """Apply a simple frequency-independent attenuation model as a fallback for underwater.

    alpha_db_per_km: absorption in dB per km (approx). This is a simplistic model.
    """
    # Convert to linear scale for amplitude loss
    loss_db = -alpha_db_per_km * (range_m / 1000.0)
    factor = db_to_amp(loss_db)
    return (sig * factor).astype(np.float32)

# If pyroomacoustics is available, provide a helper to generate an RIR and convolve
if HAS_PRA:
    def apply_room_rir(sig: np.ndarray, sr: int, room_dim: Tuple[float, float, float],
                       src_pos: Tuple[float, float, float], mic_pos: Tuple[float, float, float],
                       rt60: float = 0.3) -> np.ndarray:
        e_absorption, max_order = pra.inverse_sabine(rt60, room_dim)
        room = pra.ShoeBox(room_dim, fs=sr, materials=pra.Material(e_absorption), max_order=12)
        room.add_source(src_pos, signal=sig)
        room.add_microphone_array(pra.MicrophoneArray(np.array([mic_pos]).T, room.fs))
        room.compute_rir()
        room.simulate()
        out = room.mic_array.signals[0]
        return out.astype(np.float32)

# ----------------------------- Dataset generation ------------------------

def _random_event_for_species(species: str, sr: int, clip_duration: float) -> Tuple[np.ndarray, float, str, Dict[str, Any]]:
    """Return (signal, start_time, class_label, params) for a random event of `species`.

    species: a string like 'dolphin', 'bat', 'bird_chirp', 'click_train'
    """
    if species == 'dolphin_whistle':
        dur = random.uniform(0.2, 2.0)
        f0 = random.uniform(4000, 8000)
        f1 = random.uniform(1000, 20000)
        sig = generate_chirp_signal(dur, sr, f0, f1, method=random.choice(['linear', 'quadratic']), amplitude=1.0)
    elif species == 'bat_sweep':
        dur = random.uniform(0.002, 0.01)
        f0 = random.uniform(80000, 150000)
        f1 = random.uniform(20000, 50000)
        sig = generate_chirp_signal(dur, sr, f0, f1, method='linear', amplitude=1.0)
    elif species == 'click_train':
        num_clicks = random.randint(5, 50)
        mean_ici = random.uniform(0.005, 0.1)
        click_fc = random.uniform(20000, min(200000, sr // 2 - 1000))
        sig = generate_click_train(sr, num_clicks=num_clicks, mean_ici=mean_ici, click_fc=click_fc, amplitude=1.0)
        dur = len(sig) / sr
    elif species == 'bird_syllable':
        dur = random.uniform(0.05, 0.5)
        f0 = random.uniform(800, 4000)
        sig = generate_harmonic_tone(sr, dur, f0, n_harmonics=random.randint(1, 5))
    else:
        # fallback: a simple tone
        dur = random.uniform(0.1, 1.0)
        f0 = random.uniform(300, min(30000, sr // 2 - 1000))
        t = np.linspace(0, dur, int(round(sr * dur)), endpoint=False)
        sig = (np.sin(2 * np.pi * f0 * t) * np.sin(np.pi * t / dur)).astype(np.float32)

    # choose random start_time
    max_start = max(0.0, clip_duration - (len(sig) / sr) - 1e-6)
    start_time = random.uniform(0.0, max_start)
    params = {'species': species, 'duration': len(sig) / sr}
    return sig, start_time, species, params


def generate_dataset(output_dir: str,
                     n_clips: int = 100,
                     clip_duration: float = 10.0,
                     sr: int = 192000,
                     species_list: Optional[List[str]] = None,
                     backgrounds: Optional[List[str]] = None,
                     target_snr_db_range: Tuple[float, float] = (-12.0, -3.0),
                     events_per_clip: Tuple[int, int] = (1, 4),
                     seed: Optional[int] = None) -> None:
    """Generate a dataset of synthetic clips with labels.

    Outputs:
      output_dir/wavs/*.wav
      output_dir/manifest.jsonl  (one json per clip: {wav, sr, events: [...]})
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    ensure_dir(output_dir)
    wav_dir = os.path.join(output_dir, 'wavs')
    ensure_dir(wav_dir)

    manifest_path = os.path.join(output_dir, 'manifest.jsonl')

    if species_list is None:
        species_list = ['dolphin_whistle', 'click_train', 'bat_sweep', 'bird_syllable']

    # Preload background audio files if provided
    bg_signals = []
    if backgrounds is not None:
        for bg in backgrounds:
            data, bsr = sf.read(bg)
            if data.ndim > 1:
                data = data.mean(axis=1)
            bg_signals.append((data.astype(np.float32), bsr))

    with open(manifest_path, 'w', encoding='utf8') as mf:
        for i in range(n_clips):
            # choose background randomly
            if len(bg_signals) > 0:
                bg_idx = random.randrange(len(bg_signals))
                bg_data, bg_sr = bg_signals[bg_idx]
                # if background sample rate differs from target SR, skip for simplicity
                if bg_sr != sr:
                    # naive resample by simple interpolation (small overhead)
                    src_times = np.linspace(0, len(bg_data) / bg_sr, num=len(bg_data), endpoint=False)
                    tgt_len = int(round(clip_duration * sr))
                    tgt_times = np.linspace(0, clip_duration, num=tgt_len, endpoint=False)
                    bg = np.interp(tgt_times, src_times, bg_data).astype(np.float32)
                else:
                    # slice or tile background to clip_duration
                    need = int(round(clip_duration * sr))
                    if len(bg_data) >= need:
                        st = random.randint(0, len(bg_data) - need)
                        bg = bg_data[st:st + need].astype(np.float32)
                    else:
                        # tile
                        reps = int(math.ceil(need / len(bg_data)))
                        bg = np.tile(bg_data, reps)[:need].astype(np.float32)
            else:
                bg = None
                bg_sr = None

            n_events = random.randint(events_per_clip[0], events_per_clip[1])
            events = []
            for _ in range(n_events):
                species = random.choice(species_list)
                sig, start_time, label, params = _random_event_for_species(species, sr, clip_duration)
                events.append((sig, start_time, label, params))

            target_snr_db = random.uniform(target_snr_db_range[0], target_snr_db_range[1])
            mixed, metas = mix_events_into_background(bg, bg_sr, events, clip_duration, sr, target_snr_db=target_snr_db)

            # optional: small random RIR or underwater attenuation -- hooks
            # (Here we leave it as mixed; user can apply apply_rir_convolution or apply_room_rir manually)

            wav_name = f'clip_{i:06d}.wav'
            wav_path = os.path.join(wav_dir, wav_name)
            sf.write(wav_path, mixed, sr, subtype='FLOAT')

            manifest_item = {
                'wav': wav_path,
                'sr': sr,
                'duration': clip_duration,
                'events': [asdict(m) for m in metas]
            }
            mf.write(json.dumps(manifest_item) + '\n')

            if (i + 1) % 10 == 0:
                print(f'Generated {i+1}/{n_clips} clips')

    print('Dataset generation complete. Manifest written to', manifest_path)

# ----------------------------- Example / CLI -----------------------------

def example_generate(out_dir: str = 'synth_dataset', n_clips: int = 100):
    """Quick example that generates n_clips to out_dir using defaults."""
    species = ['dolphin_whistle', 'click_train', 'bat_sweep', 'bird_syllable']
    generate_dataset(out_dir, n_clips=n_clips, clip_duration=10.0, sr=192000, species_list=species, backgrounds=None)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='Synthetic bioacoustics dataset generator (single-file scaffold)')
    p.add_argument('--out', type=str, default='synth_dataset')
    p.add_argument('--n', type=int, default=100)
    p.add_argument('--sr', type=int, default=192000)
    p.add_argument('--dur', type=float, default=10.0)
    args = p.parse_args()
    generate_dataset(args.out, n_clips=args.n, clip_duration=args.dur, sr=args.sr)
