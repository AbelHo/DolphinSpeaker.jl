import importlib.util
import json
from html.parser import HTMLParser
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[2] / 'src/clicker_calibration/interactive_review.py'
spec = importlib.util.spec_from_file_location('interactive_review', SOURCE)
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


def row(index=459, **values):
    result = dict(review_index=str(index), frame_id=f'frame-{index}',
                  event_id=f'event-{index}', image=f'frames/{index}.png',
                  sync_status='mapped', frame_pts_s='2.5', uncertainty_px='35',
                  terra_x='100', terra_y='200', terra_status='uncertain',
                  calibrated_x='110', calibrated_y='220', calibrated_status='ok',
                  old_x='-500', old_y='900', old_status='outside_image')
    return dict(result, **values)


class InteractiveReviewTests(unittest.TestCase):
    def test_slider_zero_and_one_step_per_detection(self):
        payload = review.review_payload([row(), row(460)], 'A')
        figure = review.figure_spec(payload)
        class Inputs(HTMLParser):
            def handle_starttag(self, tag, attrs):
                attributes = dict(attrs)
                if tag == 'input' and attributes.get('id') == 'detection-slider':
                    self.slider = attributes
        parser = Inputs()
        parser.feed(review.document(payload, figure))
        self.assertEqual({k:parser.slider[k] for k in ('type','min','max','step','value')},
                         dict(type='range',min='0',max='2',step='1',value='0'))
        self.assertNotIn('sliders', figure['layout'])
        self.assertEqual([e['review_index'] for e in payload['events']], [459, 460])
        for trace in figure['data'][5:]:
            self.assertEqual(trace['x'], [])
            self.assertEqual(trace['y'], [])

    def test_missing_and_outside_coordinates_are_preserved(self):
        payload = review.review_payload([row(terra_x='', terra_y='', calibrated_x='nan')], 'B')
        points = payload['events'][0]['points']
        self.assertEqual([(p['x'], p['y']) for p in points], [(None, None), (None, None), (-500., 900.)])
        self.assertEqual(points[2]['location'], 'outside image')
        json.dumps(payload, allow_nan=False)

    def test_duplicate_identity_and_unsafe_image_path_rejected(self):
        with self.assertRaises(ValueError):
            review.review_payload([row(), row()], 'A')
        with self.assertRaises(ValueError):
            review.review_payload([row(image='../elsewhere.png')], 'A')

    def test_hover_click_identity_is_one_based_slider_identity(self):
        payload = review.review_payload([row(), row(460)], 'A')
        figure = review.figure_spec(payload)
        for trace in figure['data'][1:4]:
            self.assertEqual([p[0] for p in trace['customdata']], [1, 2])
            self.assertEqual(trace['customdata'][0][1:3], ['00459', 'frame-459'])
        self.assertEqual(figure['data'][3]['x'], [-500., -500.])
        self.assertEqual(figure['layout']['yaxis']['range'], [1520, 0])

    def test_embedded_payload_is_safe_and_offline(self):
        payload = review.review_payload([row(terra_status='</script><script>bad()</script>')], 'A')
        text = review.document(payload, review.figure_spec(payload))
        embedded = text.split('<script id="review-data" type="application/json">', 1)[1].split('</script>', 1)[0]
        self.assertEqual(json.loads(embedded), payload)
        self.assertNotIn('</script>', embedded)
        self.assertIn("byId('detection-slider').addEventListener('input'", text)
        self.assertNotIn('<script src=', text)


if __name__ == '__main__':
    unittest.main()
