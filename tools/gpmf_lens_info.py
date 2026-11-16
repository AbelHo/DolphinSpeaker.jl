#!/usr/bin/env python3
"""Read-only GoPro lens-mode and camera-identity extractor (Python 3 stdlib + ffprobe/ffmpeg).

For each input MP4 prints one JSON line:
  (a) ffprobe: width, height, fps, duration_s, creation_time, gpmd_index (data stream tagged 'gpmd');
  (b) the first 5 s of the gpmd stream is read into memory (ffmpeg -t 5 -c copy, pipe output);
  (c) the GPMF KLV tree is parsed; VFOV (lens-mode letter), ZFOV (float) and every string tag are collected;
  (d) the MP4 box tree is walked by header only (seek, no full read) to moov/udta, from which the
      GoPro boxes FIRM, LENS, CAME, SETT, MUID and HMMT are reported (hex; FIRM/LENS also as ascii).

Nothing is written next to the input files. Optional outputs go only where --jsonl/--csv point.

Usage:
  python3 tools/gpmf_lens_info.py [--jsonl OUT.jsonl] [--csv OUT.csv] FILE.MP4 [FILE2.MP4 ...]
"""
import argparse
import csv
import json
import os
import re
import struct
import subprocess
import sys

GPMF_SECONDS = 5

# ---------------------------------------------------------------- GPMF (KLV) parsing

_NUM_FMT = {
    'b': ('b', 1), 'B': ('B', 1), 's': ('>h', 2), 'S': ('>H', 2),
    'l': ('>i', 4), 'L': ('>I', 4), 'j': ('>q', 8), 'J': ('>Q', 8),
    'f': ('>f', 4), 'd': ('>d', 8),
}


def _decode_payload(type_char, size, repeat, payload):
    """Decode one GPMF payload. Returns a python value (str, number, list or hex string)."""
    if type_char == 'c':
        return payload.rstrip(b'\x00').decode('latin-1', errors='replace').strip()
    if type_char in _NUM_FMT:
        fmt, width = _NUM_FMT[type_char]
        vals = []
        for i in range(repeat):
            chunk = payload[i * size:(i + 1) * size]
            if len(chunk) < width:
                break
            vals.append(struct.unpack(fmt, chunk[:width])[0])
        return vals[0] if len(vals) == 1 else vals
    # U (date), F (fourcc), G (guid), ? (complex) and unknown types: keep the bytes
    return payload.hex()


def parse_klv(data, out, depth=0, max_depth=12):
    """Walk a GPMF byte string. Collects into out[key] = list of distinct decoded values (first 20)."""
    pos = 0
    n = len(data)
    while pos + 8 <= n:
        key = data[pos:pos + 4].decode('latin-1', errors='replace')
        type_char = chr(data[pos + 4])
        size = data[pos + 5]
        repeat = struct.unpack('>H', data[pos + 6:pos + 8])[0]
        pos += 8
        length = size * repeat
        padded = (length + 3) & ~3
        if pos + length > n:
            out.setdefault('_truncated', []).append(key)
            break
        payload = data[pos:pos + length]
        if type_char == '\x00':
            if depth < max_depth:
                parse_klv(payload, out, depth + 1, max_depth)
        else:
            value = _decode_payload(type_char, size, repeat, payload)
            vals = out.setdefault(key, [])
            if value not in vals and len(vals) < 20:
                vals.append(value)
        pos += padded
    return out


# ---------------------------------------------------------------- MP4 box walking (headers only)

def iter_boxes(f, start, end):
    """Yield (type, payload_start, payload_end, box_start) for boxes in [start, end). Seeks only."""
    pos = start
    while pos + 8 <= end:
        f.seek(pos)
        hdr = f.read(8)
        if len(hdr) < 8:
            return
        size, btype = struct.unpack('>I4s', hdr)
        btype = btype.decode('latin-1', errors='replace')
        header_len = 8
        if size == 1:
            ext = f.read(8)
            if len(ext) < 8:
                return
            size = struct.unpack('>Q', ext)[0]
            header_len = 16
        elif size == 0:
            size = end - pos
        if size < header_len:
            return
        yield btype, pos + header_len, min(pos + size, end), pos
        pos += size


def find_udta_boxes(f, file_size):
    """Return {child_type: payload_bytes} for moov/udta children. Reads only small payloads."""
    result = {}
    for btype, p0, p1, _ in iter_boxes(f, 0, file_size):
        if btype != 'moov':
            continue
        for mtype, m0, m1, _ in iter_boxes(f, p0, p1):
            if mtype != 'udta':
                continue
            for utype, u0, u1, _ in iter_boxes(f, m0, m1):
                if u1 - u0 > 4 * 1024 * 1024:
                    result[utype] = None
                    continue
                f.seek(u0)
                result[utype] = f.read(u1 - u0)
    return result


def _box_info(payload, btype):
    if payload is None:
        return None
    info = {'hex': payload.hex(), 'bytes': len(payload)}
    if btype in ('FIRM', 'LENS'):
        txt = payload.rstrip(b'\x00').decode('latin-1', errors='replace')
        if txt and all(32 <= ord(c) < 127 for c in txt):
            info['ascii'] = txt.strip()
    return info


def read_udta(path, file_size=None):
    """Return (box_info, child_types, udta_gpmf_tags).

    udta_gpmf_tags is the GPMF KLV summary stored in moov/udta/GPMF (small, in the file header).
    On HERO9 it holds VFOV/ZFOV (lens mode and zoom FOV), MINF, PRJT, CASN, LINF, CINF, FMWR, ...
    The lens-mode tags are not in the first seconds of the gpmd stream on this camera.
    """
    if file_size is None:
        file_size = os.path.getsize(path)
    with open(path, 'rb') as f:
        boxes = find_udta_boxes(f, file_size)
    out = {}
    for key in ('FIRM', 'LENS', 'CAME', 'SETT', 'MUID', 'HMMT'):
        if key in boxes:
            out[key] = _box_info(boxes[key], key)
    udta_gpmf = None
    if boxes.get('GPMF'):
        udta_gpmf = parse_klv(boxes['GPMF'], {})
        udta_gpmf.pop('_truncated', None)
    return out, sorted(boxes.keys()), udta_gpmf


# ---------------------------------------------------------------- ffprobe / ffmpeg

def probe(path):
    cmd = ['ffprobe', '-v', 'error', '-print_format', 'json', '-show_format', '-show_streams', path]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError('ffprobe failed: ' + res.stderr.strip()[:300])
    return json.loads(res.stdout)


def _fps(stream):
    rate = stream.get('avg_frame_rate') or stream.get('r_frame_rate') or '0/1'
    try:
        num, den = rate.split('/')
        return round(float(num) / float(den), 4) if float(den) else None
    except ValueError:
        return None


def read_gpmf_head(path, idx, seconds=GPMF_SECONDS):
    """Return the first `seconds` of the gpmd data stream as bytes (in memory only)."""
    base = ['ffmpeg', '-v', 'error', '-t', str(seconds), '-i', path, '-map', '0:%d' % idx, '-c', 'copy']
    for fmt in ('rawvideo', 'data'):
        res = subprocess.run(base + ['-f', fmt, 'pipe:1'], capture_output=True)
        if res.returncode == 0 and res.stdout:
            return res.stdout, fmt, None
        last_err = res.stderr.decode('utf-8', 'replace').strip()[:300]
    return b'', None, last_err


# ---------------------------------------------------------------- per-file record

_NAME_RE = re.compile(r'^G[HX](\d{2})(\d{4})\.MP4$', re.IGNORECASE)


def info_for(path):
    rec = {'path': path, 'filename': os.path.basename(path), 'chapter': None, 'counter': None,
           'width': None, 'height': None, 'fps': None, 'duration_s': None, 'creation_time': None,
           'gpmd_index': None, 'VFOV': None, 'VFOV_distinct': None, 'ZFOV': None, 'ZFOV_distinct': None,
           'DVNM': None, 'MINF': None, 'FIRM': None, 'LENS': None, 'CAME': None, 'SETT': None,
           'MUID': None, 'HMMT': None, 'CASN': None, 'LINF': None, 'FMWR': None, 'PRJT': None,
           'VFOV_source': None, 'ZFOV_source': None,
           'gpmf_string_tags': None, 'udta_gpmf': None, 'udta_children': None,
           'gpmf_format': None, 'errors': []}
    m = _NAME_RE.match(rec['filename'])
    if m:
        rec['chapter'] = int(m.group(1))
        rec['counter'] = int(m.group(2))
    try:
        pr = probe(path)
    except (RuntimeError, OSError, ValueError) as e:
        rec['errors'].append(str(e))
        return rec
    fmt = pr.get('format', {})
    rec['duration_s'] = round(float(fmt['duration']), 3) if fmt.get('duration') else None
    rec['creation_time'] = (fmt.get('tags') or {}).get('creation_time')
    for s in pr.get('streams', []):
        if s.get('codec_type') == 'video' and rec['width'] is None:
            rec['width'], rec['height'] = s.get('width'), s.get('height')
            rec['fps'] = _fps(s)
            if rec['creation_time'] is None:
                rec['creation_time'] = (s.get('tags') or {}).get('creation_time')
        if s.get('codec_tag_string') == 'gpmd' and rec['gpmd_index'] is None:
            rec['gpmd_index'] = s.get('index')
    if rec['gpmd_index'] is None:
        rec['errors'].append('no gpmd stream')
    else:
        data, used_fmt, err = read_gpmf_head(path, rec['gpmd_index'])
        rec['gpmf_format'] = used_fmt
        if used_fmt is None:
            rec['errors'].append('gpmf read failed: ' + str(err))
        else:
            tags = parse_klv(data, {})
            tags.pop('_truncated', None)
            _set_gpmf_fields(rec, tags, 'stream_head')
            rec['gpmf_string_tags'] = {k: v for k, v in tags.items()
                                       if all(isinstance(x, str) for x in v) and k not in ('DVNM', 'MINF')}
    try:
        udta, children, udta_gpmf = read_udta(path)
        rec['udta_children'] = children
        for k, v in udta.items():
            rec[k] = v
        if udta_gpmf is not None:
            rec['udta_gpmf'] = udta_gpmf
            # the udta summary is authoritative for lens mode; fill only the gaps from the stream head
            _set_gpmf_fields(rec, udta_gpmf, 'udta/GPMF', overwrite=True)
            rec['CASN'] = _first(udta_gpmf.get('CASN'))
            rec['LINF'] = _first(udta_gpmf.get('LINF'))
            rec['FMWR'] = _first(udta_gpmf.get('FMWR'))
            rec['PRJT'] = _first(udta_gpmf.get('PRJT'))
    except OSError as e:
        rec['errors'].append('udta read failed: ' + str(e))
    return rec


def _first(vals):
    return vals[0] if vals else None


def _set_gpmf_fields(rec, tags, source, overwrite=False):
    """Copy VFOV/ZFOV (with distinct values and source) and DVNM/MINF from a parsed GPMF tag dict."""
    for k in ('VFOV', 'ZFOV'):
        vals = tags.get(k)
        if vals and (overwrite or rec.get(k) is None):
            rec[k] = vals[0]
            rec[k + '_distinct'] = vals
            rec[k + '_source'] = source
    for k in ('DVNM', 'MINF'):
        if tags.get(k) and (overwrite or rec.get(k) is None):
            rec[k] = tags[k][0]


# ---------------------------------------------------------------- CLI

CSV_FIELDS = ['path', 'filename', 'chapter', 'counter', 'width', 'height', 'fps', 'duration_s',
              'creation_time', 'gpmd_index', 'VFOV', 'VFOV_source', 'ZFOV', 'ZFOV_source', 'DVNM', 'MINF',
              'FIRM', 'LENS', 'CAME', 'SETT', 'MUID', 'HMMT', 'CASN', 'LINF', 'FMWR', 'PRJT', 'errors']


def _csv_row(rec):
    row = {}
    for k in CSV_FIELDS:
        v = rec.get(k)
        if k in ('FIRM', 'LENS', 'CAME', 'SETT', 'MUID', 'HMMT'):
            # ascii when printable (FIRM, LENS), otherwise the full hex payload
            v = ((v or {}).get('ascii') or (v or {}).get('hex')) if v else None
        if k == 'ZFOV' and isinstance(v, float):
            v = round(v, 3)
        if k == 'errors':
            v = '; '.join(rec['errors'])
        row[k] = '' if v is None else v
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+')
    ap.add_argument('--jsonl', help='also write one JSON line per file to this path')
    ap.add_argument('--csv', help='also write a CSV summary to this path')
    args = ap.parse_args(argv)

    records = []
    for p in args.files:
        rec = info_for(os.path.abspath(p))
        records.append(rec)
        print(json.dumps(rec, sort_keys=True))
        sys.stdout.flush()

    if args.jsonl:
        with open(args.jsonl, 'w') as fh:
            for rec in records:
                fh.write(json.dumps(rec, sort_keys=True) + '\n')
    if args.csv:
        with open(args.csv, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
            w.writeheader()
            for rec in records:
                w.writerow(_csv_row(rec))
    return 0


if __name__ == '__main__':
    sys.exit(main())
