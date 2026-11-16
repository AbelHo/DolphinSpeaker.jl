#!/usr/bin/env python3
"""Statistics for GoPro CORI/IORI orientation CSVs written by tools/gpmf_lens_info.py --orientation-dir.

Conventions (documented; the angle is the same for either composition order):
- quaternions are (w, x, y, z), unit norm, body/camera frame.
- reference orientation of a series = component-wise median after sign alignment to the first sample,
  then normalized (a practical "median orientation").
- angle(q) = 2*acos(|w|), in degrees. The angle of a relative rotation is the rotation angle regardless
  of whether it is written CORI*conj(IORI) or conj(IORI)*CORI (they are conjugate, same angle).
- rel = CORI * conj(IORI): the rotation that takes the image orientation to the camera orientation.
  Its angle is the applied image rotation relative to the body (when CORI is the body orientation and
  IORI the stabilized image orientation); the sign/order convention is not verified against GoPro docs.
- axis split: rotation vector r = angle * (x,y,z)/|(x,y,z)| (body frame, w >= 0), per axis in degrees.
  We label the axes x, y, z (roll/pitch/yaw is a labelling we do not verify).

Usage (repo .venv):
  .venv/bin/python -I tools/orientation_stats.py --csv-dir calibration_runs/intrinsics/orientation \
      [--events A=path/to/events.csv --events B=path/to/events.csv] [--json-out OUT.json]
"""
import argparse
import bisect
import csv
import json
import math
import os
import sys


# ---------------------------------------------------------------- quaternion helpers

def qnorm(q):
    n = math.sqrt(sum(x * x for x in q))
    return [x / n for x in q] if n > 0 else [1.0, 0.0, 0.0, 0.0]


def qmul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return [aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw]


def qconj(q):
    return [q[0], -q[1], -q[2], -q[3]]


def qangle_deg(q):
    """Rotation angle of a unit quaternion in degrees, in [0, 180]."""
    w = min(1.0, abs(q[0]))
    return math.degrees(2.0 * math.acos(w))


def qrotvec_deg(q):
    """Rotation vector (x, y, z) in degrees: angle times unit axis, with w >= 0."""
    if q[0] < 0:
        q = [-x for x in q]
    s = math.sqrt(q[1] ** 2 + q[2] ** 2 + q[3] ** 2)
    ang = math.degrees(2.0 * math.atan2(s, q[0]))
    if s < 1e-12:
        return [0.0, 0.0, 0.0]
    return [ang * q[1] / s, ang * q[2] / s, ang * q[3] / s]


def qdot(a, b):
    return sum(x * y for x, y in zip(a, b))


def qmedian_ref(qs, max_candidates=300):
    """Orientation 'median' = medoid: the sample whose median angular distance to all samples is smallest.

    Component-wise medians are not used: when a series spans tens of degrees, the averaged quaternion is
    not a good centre. Candidates are an even subsample (at most max_candidates) of the samples.
    """
    if not qs:
        return [1.0, 0.0, 0.0, 0.0]
    stride = max(1, len(qs) // max_candidates)
    best, best_val = qs[0], float('inf')
    for k in range(0, len(qs), stride):
        cand = qs[k]
        d = [2.0 * math.acos(min(1.0, abs(qdot(cand, q)))) for q in qs]
        m = median(d)
        if m < best_val:
            best, best_val = cand, m
    return qnorm(best)


def median(vals):
    v = sorted(vals)
    n = len(v)
    if n == 0:
        return float('nan')
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def percentile(vals, p):
    """Linear-interpolated percentile, p in [0, 100]."""
    v = sorted(vals)
    if not v:
        return float('nan')
    k = (len(v) - 1) * p / 100.0
    lo = math.floor(k)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


# ---------------------------------------------------------------- CSV in

def load_orientation_csv(path):
    """Return list of (t_s, cori[4], iori[4]) from a *_orientation.csv file."""
    rows = []
    with open(path, newline='') as fh:
        rd = csv.reader(fh)
        next(rd)
        for r in rd:
            t = float(r[0])
            c = [float(x) for x in r[1:5]]
            i = [float(x) for x in r[5:9]]
            rows.append((t, c, i))
    return rows


# ---------------------------------------------------------------- per-clip statistics

def clip_stats(rows):
    cori = [qnorm(c) for _, c, _ in rows]
    iori = [qnorm(i) for _, _, i in rows]
    ref_c = qmedian_ref(cori)
    ref_i = qmedian_ref(iori)
    rel = [qnorm(qmul(c, qconj(i))) for c, i in zip(cori, iori)]
    ref_rel = qmedian_ref(rel)
    ang_c = [qangle_deg(qmul(c, qconj(ref_c))) for c in cori]
    ang_i = [qangle_deg(qmul(i, qconj(ref_i))) for i in iori]
    ang_i_ident = [qangle_deg(i) for i in iori]
    ang_rel_abs = [qangle_deg(r) for r in rel]
    ang_rel_vs_med = [qangle_deg(qmul(r, qconj(ref_rel))) for r in rel]
    rv = [qrotvec_deg(qmul(r, qconj(ref_rel))) for r in rel]
    axis = {}
    for k, name in enumerate(('x', 'y', 'z')):
        comp = [abs(v[k]) for v in rv]
        axis[name] = {'p95_abs_deg': percentile(comp, 95), 'max_abs_deg': max(comp) if comp else float('nan'),
                      'rms_deg': math.sqrt(sum(c * c for c in comp) / len(comp)) if comp else float('nan')}
    total = sum(a['rms_deg'] ** 2 for a in axis.values())
    for name in axis:
        axis[name]['share_of_rms2'] = (axis[name]['rms_deg'] ** 2 / total) if total > 0 else float('nan')
    # IORI is the image orientation relative to the body (the applied image rotation); axis split of IORI
    iv = [qrotvec_deg(qmul(i, qconj(ref_i))) for i in iori]
    iori_axis = {}
    for k, name in enumerate(('x', 'y', 'z')):
        comp = [abs(v[k]) for v in iv]
        iori_axis[name] = {'p95_abs_deg': percentile(comp, 95), 'max_abs_deg': max(comp) if comp else float('nan'),
                           'rms_deg': math.sqrt(sum(c * c for c in comp) / len(comp)) if comp else float('nan')}
    itot = sum(a['rms_deg'] ** 2 for a in iori_axis.values())
    for name in iori_axis:
        iori_axis[name]['share_of_rms2'] = (iori_axis[name]['rms_deg'] ** 2 / itot) if itot > 0 else float('nan')
    return {
        'iori_axes_vs_median_deg': iori_axis,
        'n': len(rows),
        'duration_s': rows[-1][0] - rows[0][0] if rows else None,
        'cori_vs_median': _summ(ang_c),
        'iori_vs_median': _summ(ang_i),
        'iori_vs_identity': _summ(ang_i_ident),
        'rel_abs_angle': _summ(ang_rel_abs),
        'rel_vs_median': _summ(ang_rel_vs_med),
        'rel_axes_vs_median_deg': axis,
        '_ref': {'cori': ref_c, 'iori': ref_i, 'rel': ref_rel},
    }


def _summ(vals):
    return {'p50': percentile(vals, 50), 'p95': percentile(vals, 95), 'max': max(vals) if vals else float('nan')}


# ---------------------------------------------------------------- events

def event_stats(rows, clip_ref, event_times):
    """Per-event values at the sample nearest each event time (no interpolation; 60 Hz samples)."""
    ts = [r[0] for r in rows]
    cori = [qnorm(c) for _, c, _ in rows]
    iori = [qnorm(i) for _, _, i in rows]
    rel = [qnorm(qmul(c, qconj(i))) for c, i in zip(cori, iori)]
    ev_c, ev_i, ev_i_ident, ev_rel_abs, ev_rel_q = [], [], [], [], []
    ev_idx = []
    iori_q_at = iori
    for t in event_times:
        k = bisect.bisect_left(ts, t)
        cands = [j for j in (k - 1, k) if 0 <= j < len(ts)]
        j = min(cands, key=lambda j: abs(ts[j] - t))
        ev_idx.append(j)
        ev_c.append(qangle_deg(qmul(cori[j], qconj(clip_ref['cori']))))
        ev_i.append(qangle_deg(qmul(iori[j], qconj(clip_ref['iori']))))
        ev_i_ident.append(qangle_deg(iori[j]))
        ev_rel_abs.append(qangle_deg(rel[j]))
        ev_rel_q.append(rel[j])
    if not ev_rel_q:
        return None
    ref_ev = qmedian_ref(ev_rel_q)
    spread_vs_event_median = [qangle_deg(qmul(r, qconj(ref_ev))) for r in ev_rel_q]
    ev_iori_q = [iori_q_at[j] for j in ev_idx]
    ref_iev = qmedian_ref(ev_iori_q)
    iori_spread = [qangle_deg(qmul(q, qconj(ref_iev))) for q in ev_iori_q]
    return {
        'n_events': len(event_times),
        'cori_vs_clip_median_deg': _summ(ev_c),
        'iori_vs_clip_median_deg': _summ(ev_i),
        'iori_vs_identity_deg': _summ(ev_i_ident),
        'iori_vs_event_median_deg': _summ(iori_spread),
        'rel_abs_angle_deg': _summ(ev_rel_abs),
        'rel_vs_event_median_deg': _summ(spread_vs_event_median),
    }


def load_events(path):
    """Return {basename_of_video_file: [video_time_s, ...]} from an events.csv."""
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh):
            base = os.path.basename(r['video_file'])
            out.setdefault(base, []).append(float(r['video_time_s']))
    return out


# ---------------------------------------------------------------- main

def fmt(d):
    if not isinstance(d, dict):
        return str(d)
    return ' '.join('%s=%.3f' % (k, v) for k, v in d.items() if isinstance(v, float))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--csv-dir', required=True)
    ap.add_argument('--clips', nargs='*', help='basenames without _orientation.csv (default: all in dir)')
    ap.add_argument('--events', action='append', default=[],
                    help='LABEL=path/to/events.csv; events are matched to clips by video file basename')
    ap.add_argument('--json-out')
    args = ap.parse_args(argv)

    stems = args.clips or sorted(f[:-len('_orientation.csv')] for f in os.listdir(args.csv_dir)
                                 if f.endswith('_orientation.csv'))
    ev_maps = {}
    for spec in args.events:
        label, path = spec.split('=', 1)
        ev_maps[label] = load_events(path)

    report = {}
    for stem in stems:
        rows = load_orientation_csv(os.path.join(args.csv_dir, stem + '_orientation.csv'))
        st = clip_stats(rows)
        video = stem + '.MP4'
        ev = None
        for label, m in ev_maps.items():
            if video in m:
                ev = event_stats(rows, st['_ref'], m[video])
                ev['label'] = label
        st_out = {k: v for k, v in st.items() if k != '_ref'}
        st_out['events'] = ev
        report[stem] = st_out
        print('== %s  n=%d  dur=%.1fs' % (stem, st['n'], st['duration_s'] or 0))
        print('  CORI vs median  deg: %s' % fmt(st['cori_vs_median']))
        print('  IORI vs median  deg: %s' % fmt(st['iori_vs_median']))
        print('  IORI vs identity deg: %s' % fmt(st['iori_vs_identity']))
        print('  IORI axes vs median (rms deg): %s' % ' '.join(
            '%s=%.3f' % (k, v['rms_deg']) for k, v in st['iori_axes_vs_median_deg'].items()))
        print('  rel=CORI*IORI^-1 angle deg: %s' % fmt(st['rel_abs_angle']))
        print('  rel vs its median deg: %s' % fmt(st['rel_vs_median']))
        print('  rel axes (deg, rms per axis): %s' % ' '.join(
            '%s=%.3f' % (k, v['rms_deg']) for k, v in st['rel_axes_vs_median_deg'].items()))
        print('  rel axes share of rms^2: %s' % ' '.join(
            '%s=%.2f' % (k, v['share_of_rms2']) for k, v in st['rel_axes_vs_median_deg'].items()))
        if ev:
            print('  events (%s, n=%d):' % (ev['label'], ev['n_events']))
            print('    CORI vs clip median deg: %s' % fmt(ev['cori_vs_clip_median_deg']))
            print('    IORI vs clip median deg: %s' % fmt(ev['iori_vs_clip_median_deg']))
            print('    IORI vs event median deg: %s' % fmt(ev['iori_vs_event_median_deg']))
            print('    rel abs angle deg: %s' % fmt(ev['rel_abs_angle_deg']))
            print('    rel vs event-median deg: %s' % fmt(ev['rel_vs_event_median_deg']))
    if args.json_out:
        with open(args.json_out, 'w') as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
