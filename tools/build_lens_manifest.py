#!/usr/bin/env python3
"""Build the GoPro lens/identity manifest and the checkerboard clip survey (read-only on calf).

Step 2  -> calibration_runs/intrinsics/lens_manifest.jsonl and lens_manifest.csv
           One record per selected file, from tools/gpmf_lens_info.py.
           Selection per target folder:
             'gh01' : GoPro chapter-01 clips (GHcc+counter, cc=01); a counter with no chapter 01 keeps
                      its lowest chapter. Non-GoPro names are skipped.
             'all'  : every .mp4/.mov/.lrv file (all chapters, all names).
Step 3  -> calibration_runs/intrinsics/gopro_clips_filter.csv   (every GoPro clip under calf matching
                                                                   counter 0001-0300 or creation_time >= 2026-01-15)
           calibration_runs/intrinsics/checkerboard_clips.csv   (A checkerboard clips + flagged B-era candidates)

A GoPro clip is a file whose ffprobe stream list has a 'gpmd' data stream (name-independent), or a GH/GX name.
Only ffprobe (headers) and the 5 s GPMF head + udta are read. No frames are decoded. Nothing is written
next to the data. Usage:  python3 -I tools/build_lens_manifest.py
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gpmf_lens_info as g  # noqa: E402

CALF = '/media/spin/anas/data/calf'
OUT_DIR = os.path.join(os.path.dirname(HERE), 'calibration_runs', 'intrinsics')

TARGETS = [
    # (relative folder, rule, device label from folder, group)
    ('Calibration/device3_calibration/checkerboard_2025-09-19', 'all', 'A', 'A_checker_0919'),
    ('Calibration/device3_calibration/checkerboard_2025-09-22', 'all', 'A', 'A_checker_0922'),
    ('Calibration/device3_calibration', 'gh01', 'A', 'A_device3'),
    ('upload/Calibration_05_11_2025', 'gh01', 'A', 'A_calib_0511'),
    ('upload/Calibration_03_10_2025', 'gh01', 'A', 'A_calib_0310'),
    ('Clicktest/clicker_20260209/A', 'gh01', 'A', 'A_clicktest'),
    ('Clicktest/clicker_20260209/B', 'gh01', 'B', 'B_clicktest'),
    ('upload/TT-shape/TDA_06022026/A', 'gh01', 'A', 'A_TT_0602'),
    ('upload/TT-shape/TDA_06022026/B', 'gh01', 'B', 'B_TT_0602'),
    ('upload/Single_dolphin_target/13_11_2026_Tutti/A', 'gh01', 'A', 'A_tutti_1311'),
    ('upload/Single_dolphin_target/13_11_2026_Tutti/B', 'gh01', 'B', 'B_tutti_1311'),
    ('temp/delete/from_mac/cam_calib', 'all', 'mar_may_2025_camera', 'cam_calib_0318'),
    ('Calibration/device_2_calibration/GoPro', 'gh01', 'mar_may_2025_camera', 'device_2'),
    ('Calibration/Housingv3_test_28_03_2025', 'all', 'mar_may_2025_camera', 'housingv3_0328'),
]

VIDEO_RE = re.compile(r'\.(mp4|mov|lrv)$', re.IGNORECASE)
GH_RE = re.compile(r'^G[HX](\d{2})(\d{4})\.MP4$', re.IGNORECASE)
NAME_FLAG_RE = re.compile(r'checker|chess|board|calib|cam', re.IGNORECASE)

# Chain by counter (from the device facts): A = 1265-2051, B = 0059-0133, Mar-May 2025 = 98xx-99xx.
A_FOLDERS = ('device3_calibration', 'Clicktest/clicker_20260209/A', 'TDA_06022026/A', '13_11_2026_Tutti/A',
             'Calibration_05_11_2025', 'Calibration_03_10_2025')
B_FOLDERS = ('Clicktest/clicker_20260209/B', 'TDA_06022026/B', '13_11_2026_Tutti/B')


def list_video_files(root_rel):
    root = os.path.join(CALF, root_rel)
    out = []
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if VIDEO_RE.search(fn):
                out.append(os.path.join(dp, fn))
    return sorted(out)


def select_gh01(paths):
    """Keep chapter 01 per counter; fall back to the lowest chapter when no 01 exists."""
    best = {}
    others = []
    for p in paths:
        m = GH_RE.match(os.path.basename(p))
        if not m:
            continue
        key = (os.path.dirname(p), m.group(2))
        ch = int(m.group(1))
        if key not in best or ch < best[key][0]:
            best[key] = (ch, p)
    return sorted(v[1] for v in best.values())


def probe_light(path):
    """ffprobe header info only: gpmd presence, creation_time, duration."""
    pr = g.probe(path)
    fmt = pr.get('format', {})
    has_gpmd = any(s.get('codec_tag_string') == 'gpmd' for s in pr.get('streams', []))
    ctime = (fmt.get('tags') or {}).get('creation_time')
    if ctime is None:
        for s in pr.get('streams', []):
            ctime = ctime or (s.get('tags') or {}).get('creation_time')
    dur = round(float(fmt['duration']), 3) if fmt.get('duration') else None
    return has_gpmd, ctime, dur


def infer_device(path, counter):
    rel = os.path.relpath(path, CALF)
    if any(f in rel for f in B_FOLDERS):
        return 'B', 'folder'
    if any(f in rel for f in A_FOLDERS):
        return 'A', 'folder'
    if counter is not None:
        if 1265 <= counter <= 2051:
            return 'A', 'counter_chain_gamma'
        if 59 <= counter <= 133:
            return 'B', 'counter_chain_beta'
        if 9000 <= counter <= 9999:
            return 'mar_may_2025_camera', 'counter_chain_alpha'
    return 'unknown', 'none'


def write_jsonl_csv(recs, base, extra_cols=None):
    extra_cols = extra_cols or []
    jl = base + '.jsonl'
    with open(jl, 'w') as fh:
        for r in recs:
            fh.write(json.dumps(r, sort_keys=True) + '\n')
    cols = g.CSV_FIELDS + extra_cols
    with open(base + '.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in recs:
            row = g._csv_row(r)
            for c in extra_cols:
                row[c] = r.get(c, '')
            w.writerow(row)
    return jl, base + '.csv'


def step2():
    os.makedirs(OUT_DIR, exist_ok=True)
    seen = set()
    recs = []
    for rel, rule, dev, group in TARGETS:
        paths = list_video_files(rel)
        if rule == 'gh01':
            paths = select_gh01(paths)
        for p in paths:
            if p in seen:
                continue
            seen.add(p)
            rec = g.info_for(p)
            rec['target'] = group
            rec['device_folder_label'] = dev
            recs.append(rec)
            print('step2', group, os.path.basename(p), rec['counter'], rec['VFOV'], rec['errors'], flush=True)
    jl, cs = write_jsonl_csv(recs, os.path.join(OUT_DIR, 'lens_manifest'),
                             extra_cols=['target', 'device_folder_label'])
    print('wrote', jl, cs, len(recs), 'records', flush=True)
    return recs


def step3():
    cand = []
    for dp, _, fns in os.walk(CALF):
        for fn in fns:
            if VIDEO_RE.search(fn):
                cand.append(os.path.join(dp, fn))
    cand.sort()
    rows = []
    probe_fail = []
    for p in cand:
        m = GH_RE.match(os.path.basename(p))
        counter = int(m.group(2)) if m else None
        chapter = int(m.group(1)) if m else None
        try:
            has_gpmd, ctime, dur = probe_light(p)
        except (RuntimeError, OSError, ValueError) as e:
            probe_fail.append((p, str(e)[:200]))
            continue
        if not (has_gpmd or m):
            continue  # not a GoPro clip
        in_counter = counter is not None and 1 <= counter <= 300
        b_era = ctime is not None and ctime[:10] >= '2026-01-15'
        if not (in_counter or b_era):
            continue
        rel = os.path.relpath(p, CALF)
        flag = bool(NAME_FLAG_RE.search(rel))
        rows.append({'path': p, 'counter': counter, 'chapter': chapter, 'creation_time': ctime,
                     'duration_s': dur, 'name_flag': 'checker/chess/board/calib/cam' if flag else '',
                     'device': None, 'device_source': None, 'has_gpmd': has_gpmd,
                     'in_counter_0001_0300': in_counter, 'b_era': b_era})
    print('step3 candidates', len(cand), 'GoPro-filtered', len(rows), 'probe failures', len(probe_fail),
          flush=True)
    for r in rows:
        r['device'], r['device_source'] = infer_device(r['path'], r['counter'])
        full = g.info_for(r['path'])
        for k in ('VFOV', 'ZFOV', 'CAME', 'LENS', 'FIRM', 'CASN', 'LINF', 'FMWR', 'VFOV_source'):
            r[k] = full.get(k)
        r['errors'] = full.get('errors')
        r['_full'] = full

    cols = ['device', 'device_source', 'path', 'counter', 'chapter', 'creation_time', 'duration_s',
            'VFOV', 'ZFOV', 'CAME', 'LENS', 'FIRM', 'name_flag', 'in_counter_0001_0300', 'b_era', 'has_gpmd',
            'CASN', 'LINF', 'FMWR', 'VFOV_source', 'errors']

    def cell(r, k):
        v = r.get(k)
        if k in ('CAME', 'LENS', 'FIRM') and isinstance(v, dict):
            v = v.get('ascii') or v.get('hex')
        if k == 'ZFOV' and isinstance(v, float):
            v = round(v, 3)
        if k == 'errors':
            v = '; '.join(v or [])
        return '' if v is None else v

    with open(os.path.join(OUT_DIR, 'gopro_clips_filter.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([cell(r, k) for k in cols])

    # checkerboard list: A checkerboard clips (from step 2 manifest selection, all chapters) + flagged candidates
    a_check_roots = ('Calibration/device3_calibration/checkerboard_2025-09-19',
                     'Calibration/device3_calibration/checkerboard_2025-09-22')
    sel = []
    for r in rows:
        rel = os.path.relpath(r['path'], CALF)
        if rel.startswith(a_check_roots) or (r['name_flag'] and r['b_era'] and r['device'] != 'A'):
            sel.append(r)
    # A checkerboard clips are always in the list even if the filter above missed them
    have = {r['path'] for r in sel}
    for root in a_check_roots:
        for p in list_video_files(root):
            if p in have:
                continue
            has_gpmd, ctime, dur = probe_light(p)
            m = GH_RE.match(os.path.basename(p))
            full = g.info_for(p)
            sel.append({'path': p, 'counter': int(m.group(2)) if m else None,
                        'chapter': int(m.group(1)) if m else None, 'creation_time': ctime,
                        'duration_s': dur, 'name_flag': 'checker/chess/board/calib/cam',
                        'device': 'A', 'device_source': 'folder', 'has_gpmd': has_gpmd,
                        'in_counter_0001_0300': False, 'b_era': False,
                        **{k: full.get(k) for k in ('VFOV', 'ZFOV', 'CAME', 'LENS', 'FIRM', 'CASN', 'LINF',
                                                   'FMWR', 'VFOV_source')},
                        'errors': full.get('errors'), '_full': full})
    sel.sort(key=lambda r: (r['device'], r['path']))
    with open(os.path.join(OUT_DIR, 'checkerboard_clips.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        cols2 = ['device', 'path', 'counter', 'chapter', 'creation_time', 'duration_s', 'VFOV', 'ZFOV',
                 'CAME', 'LENS', 'FIRM', 'name_flag', 'device_source', 'VFOV_source', 'CASN', 'LINF', 'FMWR']
        w.writerow(cols2)
        for r in sel:
            w.writerow([cell(r, k) if k != 'path' else r['path'] for k in cols2])
    for p, e in probe_fail:
        print('probe failed', p, e, flush=True)
    print('wrote gopro_clips_filter.csv', len(rows), 'rows;', 'checkerboard_clips.csv', len(sel), 'rows',
          flush=True)
    return rows, sel, probe_fail


if __name__ == '__main__':
    which = sys.argv[1:] or ['2', '3']
    if '2' in which:
        step2()
    if '3' in which:
        step3()
