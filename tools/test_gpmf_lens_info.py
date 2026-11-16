#!/usr/bin/env python3
"""Synthetic tests for tools/gpmf_lens_info.py (no real video, no temp files).

Run: python3 -I tools/test_gpmf_lens_info.py   (or python3 tools/test_gpmf_lens_info.py)
"""
import io
import os
import struct
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gpmf_lens_info as g  # noqa: E402


def klv(key, type_char, size, repeat, payload):
    """Build one GPMF KLV entry; payload must already be size*repeat bytes."""
    assert len(payload) == size * repeat
    body = payload + b'\x00' * ((-len(payload)) % 4)
    return key.encode('ascii') + type_char.encode('ascii') + bytes([size]) + struct.pack('>H', repeat) + body


def container(key, inner):
    return klv(key, '\x00', 1, len(inner), inner)


def box(btype, payload):
    return struct.pack('>I4s', 8 + len(payload), btype.encode('ascii')) + payload


class TestKLV(unittest.TestCase):
    def test_flat_entries_and_types(self):
        data = (klv('VFOV', 'c', 1, 1, b'S')
                + klv('ZFOV', 'f', 4, 1, struct.pack('>f', 118.5))
                + klv('MINF', 'c', 1, 9, b'HERO9 Blk')
                + klv('ISO ', 'S', 2, 2, struct.pack('>HH', 100, 200)))
        out = g.parse_klv(data, {})
        self.assertEqual(out['VFOV'], ['S'])
        self.assertAlmostEqual(out['ZFOV'][0], 118.5, places=3)
        self.assertEqual(out['MINF'], ['HERO9 Blk'])
        self.assertEqual(out['ISO '], [[100, 200]])

    def test_nested_container_and_padding(self):
        # 'c' payload of length 3 is padded to 4 bytes; the following key must still parse.
        inner = klv('VFOV', 'c', 1, 1, b'W') + klv('DVNM', 'c', 1, 3, b'abc')
        data = container('DEVC', inner) + klv('ZFOV', 'f', 4, 1, struct.pack('>f', 2.5))
        out = g.parse_klv(data, {})
        self.assertEqual(out['VFOV'], ['W'])
        self.assertEqual(out['DVNM'], ['abc'])
        self.assertEqual(out['ZFOV'], [2.5])

    def test_distinct_values_and_truncation(self):
        data = klv('VFOV', 'c', 1, 1, b'W') + klv('VFOV', 'c', 1, 1, b'W') + klv('VFOV', 'c', 1, 1, b'L')
        out = g.parse_klv(data, {})
        self.assertEqual(out['VFOV'], ['W', 'L'])
        # cut the stream inside a payload: the truncated key is reported, nothing crashes
        cut = klv('ZFOV', 'f', 4, 2, struct.pack('>ff', 1.0, 2.0))[:-3]
        out2 = g.parse_klv(cut, {})
        self.assertIn('_truncated', out2)
        self.assertNotIn('ZFOV', out2)

    def test_unknown_type_kept_as_hex(self):
        data = klv('CAME', '?', 4, 1, b'\x01\x02\x03\x04')
        out = g.parse_klv(data, {})
        self.assertEqual(out['CAME'], ['01020304'])


class TestBoxes(unittest.TestCase):
    def _mp4(self):
        udta_payload = (box('FIRM', b'HD9.01.01.50.00\x00')
                        + box('LENS', b'\x00\x01')
                        + box('CAME', bytes(range(16)))
                        + box('MUID', b'\xaa\xbb\xcc\xdd'))
        moov = box('moov', box('mvhd', b'\x00' * 100) + box('udta', udta_payload))
        ftyp = box('ftyp', b'mp41\x00\x00\x00\x00mp41')
        mdat = box('mdat', b'\x11' * 32)
        return ftyp + mdat + moov

    def test_udta_children_found_by_header_walk(self):
        raw = self._mp4()
        f = io.BytesIO(raw)
        udta = g.find_udta_boxes(f, len(raw))
        self.assertEqual(sorted(udta.keys()), ['CAME', 'FIRM', 'LENS', 'MUID'])
        self.assertEqual(udta['MUID'], b'\xaa\xbb\xcc\xdd')

    def test_box_info_hex_and_ascii(self):
        raw = self._mp4()
        udta = g.find_udta_boxes(io.BytesIO(raw), len(raw))
        firm = g._box_info(udta['FIRM'], 'FIRM')
        self.assertEqual(firm['ascii'], 'HD9.01.01.50.00')
        came = g._box_info(udta['CAME'], 'CAME')
        self.assertEqual(came['hex'], bytes(range(16)).hex())
        self.assertNotIn('ascii', came)

    def test_largesize_and_zero_size_boxes(self):
        payload = b'\x00' * 8
        large = struct.pack('>I4sQ', 1, b'free', 16 + len(payload)) + payload
        raw = large + box('moov', box('udta', box('HMMT', b'\x00\x00\x00\x01')))
        udta = g.find_udta_boxes(io.BytesIO(raw), len(raw))
        self.assertEqual(udta['HMMT'], b'\x00\x00\x00\x01')
        # size 0 means "to end of file"
        tail = struct.pack('>I4s', 0, b'moov') + box('udta', box('SETT', b'\x01'))
        udta2 = g.find_udta_boxes(io.BytesIO(tail), len(tail))
        self.assertEqual(udta2['SETT'], b'\x01')

    def test_no_moov(self):
        raw = box('ftyp', b'isom')
        self.assertEqual(g.find_udta_boxes(io.BytesIO(raw), len(raw)), {})


class TestUdtaGpmf(unittest.TestCase):
    """moov/udta/GPMF holds the lens-mode summary; the stream head value is only a fallback."""

    def test_udta_gpmf_lens_tags_override_stream_head(self):
        gpmf = (klv('VFOV', 'c', 1, 1, b'S')
                + klv('ZFOV', 'f', 4, 1, struct.pack('>f', 122.5))
                + klv('MINF', 'c', 1, 9, b'HERO9 Blk'))
        raw = box('moov', box('udta', box('GPMF', gpmf) + box('CAME', bytes(16))))
        udta_tags = g.parse_klv(g.find_udta_boxes(io.BytesIO(raw), len(raw))['GPMF'], {})
        rec = {'VFOV': None, 'ZFOV': None, 'DVNM': None, 'MINF': None}
        g._set_gpmf_fields(rec, {'VFOV': ['W'], 'DVNM': ['Highlights']}, 'stream_head')
        self.assertEqual(rec['VFOV'], 'W')
        self.assertEqual(rec['VFOV_source'], 'stream_head')
        g._set_gpmf_fields(rec, udta_tags, 'udta/GPMF', overwrite=True)
        self.assertEqual(rec['VFOV'], 'S')
        self.assertEqual(rec['VFOV_source'], 'udta/GPMF')
        self.assertAlmostEqual(rec['ZFOV'], 122.5, places=3)
        self.assertEqual(rec['MINF'], 'HERO9 Blk')
        self.assertEqual(rec['DVNM'], 'Highlights')

    def test_stream_head_fills_gap_when_udta_absent(self):
        rec = {'VFOV': None, 'ZFOV': None, 'DVNM': None, 'MINF': None}
        g._set_gpmf_fields(rec, {'ZFOV': [2.0]}, 'stream_head')
        self.assertEqual(rec['ZFOV'], 2.0)
        self.assertIsNone(rec['VFOV'])


class TestNames(unittest.TestCase):
    def test_gopro_name_pattern(self):
        m = g._NAME_RE.match('GH011265.MP4')
        self.assertEqual((m.group(1), m.group(2)), ('01', '1265'))
        m = g._NAME_RE.match('gh021284.mp4'.upper())
        self.assertEqual((m.group(1), m.group(2)), ('02', '1284'))
        self.assertIsNone(g._NAME_RE.match('GH011326_short.mp4'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
