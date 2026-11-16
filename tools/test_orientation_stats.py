#!/usr/bin/env python3
"""Synthetic tests for tools/orientation_stats.py (no files except in-memory data).

Run: .venv/bin/python -I tools/test_orientation_stats.py
"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orientation_stats as o  # noqa: E402


def axis_angle_q(axis, deg):
    a = math.radians(deg) / 2.0
    s = math.sin(a)
    return [math.cos(a), axis[0] * s, axis[1] * s, axis[2] * s]


class TestQuaternions(unittest.TestCase):
    def test_angle_of_known_rotations(self):
        self.assertAlmostEqual(o.qangle_deg(axis_angle_q((0, 0, 1), 90)), 90.0, places=9)
        self.assertAlmostEqual(o.qangle_deg(axis_angle_q((1, 0, 0), 30)), 30.0, places=9)
        # the sign of q does not change the angle
        q = axis_angle_q((0, 1, 0), 40)
        self.assertAlmostEqual(o.qangle_deg([-x for x in q]), 40.0, places=9)

    def test_rotvec_axis_and_angle(self):
        rv = o.qrotvec_deg(axis_angle_q((0, 0, 1), 25))
        self.assertAlmostEqual(rv[0], 0.0, places=9)
        self.assertAlmostEqual(rv[1], 0.0, places=9)
        self.assertAlmostEqual(rv[2], 25.0, places=9)
        # a negative rotation about z has a negative z component
        self.assertLess(o.qrotvec_deg(axis_angle_q((0, 0, 1), -10))[2], 0)

    def test_composition_angle_is_order_independent(self):
        c = axis_angle_q((1, 2, 0.5), 70)
        i = axis_angle_q((0, 1, 0), 13)
        a1 = o.qangle_deg(o.qmul(c, o.qconj(i)))
        a2 = o.qangle_deg(o.qmul(o.qconj(i), c))
        self.assertAlmostEqual(a1, a2, places=9)


class TestReferenceAndStats(unittest.TestCase):
    def test_medoid_ignores_outlier_and_is_a_sample(self):
        base = axis_angle_q((0, 0, 1), 10)
        qs = [base] * 50 + [axis_angle_q((0, 0, 1), 12)] * 10 + [axis_angle_q((0, 0, 1), 170)] * 2
        ref = o.qmedian_ref(qs)
        self.assertAlmostEqual(o.qangle_deg(o.qmul(ref, o.qconj(base))), 0.0, places=9)

    def test_percentile_and_median(self):
        self.assertEqual(o.median([3, 1, 2]), 2)
        self.assertAlmostEqual(o.percentile([0, 10], 50), 5.0)
        self.assertAlmostEqual(o.percentile(list(range(101)), 95), 95.0)

    def test_clip_stats_static_camera_has_zero_spread(self):
        rows = [(i / 60.0, axis_angle_q((0, 0, 1), 0)[:], axis_angle_q((0, 1, 0), 2.0)) for i in range(120)]
        st = o.clip_stats(rows)
        self.assertAlmostEqual(st['cori_vs_median']['max'], 0.0, places=9)
        self.assertAlmostEqual(st['iori_vs_median']['max'], 0.0, places=9)
        self.assertAlmostEqual(st['iori_vs_identity']['p50'], 2.0, places=9)
        # IORI is constant, so its spread about its own median is zero on every axis (share undefined: NaN)
        self.assertAlmostEqual(st['iori_axes_vs_median_deg']['y']['rms_deg'], 0.0, places=9)
        self.assertTrue(math.isnan(st['iori_axes_vs_median_deg']['y']['share_of_rms2']))

    def test_event_lookup_nearest_sample(self):
        rows = [(i / 60.0, axis_angle_q((0, 0, 1), i), axis_angle_q((0, 0, 1), 0)) for i in range(120)]
        st = o.clip_stats(rows)
        ev = o.event_stats(rows, st['_ref'], [0.51, 1.0])
        self.assertEqual(ev['n_events'], 2)
        self.assertGreater(ev['cori_vs_clip_median_deg']['max'], 0.0)

    def test_load_events_basename(self):
        # scratch goes under the git-ignored calibration_runs/, never /tmp
        here = os.path.dirname(os.path.abspath(__file__))
        scratch = os.path.join(os.path.dirname(here), 'calibration_runs', 'intrinsics', 'orientation')
        os.makedirs(scratch, exist_ok=True)
        name = os.path.join(scratch, '_test_events.csv')
        with open(name, 'w') as fh:
            fh.write('event_id,video_file,video_time_s\n')
            fh.write('e1,/media/x/GH011574.MP4,20.5\n')
            fh.write('e2,/media/x/GH011574.MP4,21.0\n')
        try:
            m = o.load_events(name)
        finally:
            os.unlink(name)
        self.assertEqual(m, {'GH011574.MP4': [20.5, 21.0]})


if __name__ == '__main__':
    unittest.main(verbosity=2)
