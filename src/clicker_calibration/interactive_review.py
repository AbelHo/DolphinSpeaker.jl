#!/usr/bin/env python3
"""Create an offline Plotly slider review from a three-way positions.csv.

python src/clicker_calibration/interactive_review.py DEVICE_COMPARISON_DIRECTORY
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path

import plotly
from plotly.offline import get_plotlyjs

METHODS = (('terra', 'Terra', '#ffbf32', 'circle'),
           ('calibrated', 'Calibrated', '#ff4545', 'x'),
           ('old', 'Old D3', '#20e5f3', 'diamond'))
WIDTH, HEIGHT = 2704, 1520


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def review_payload(rows, device):
    events = []
    for i, row in enumerate(rows, 1):
        image = Path(row['image'])
        if image.is_absolute() or '..' in image.parts:
            raise ValueError('Image paths must stay within the review directory')
        points = []
        for key, label, color, symbol in METHODS:
            x, y = finite(row[key + '_x']), finite(row[key + '_y'])
            if x is None or y is None:
                x = y = None
            points.append(dict(method=key, label=label, color=color, symbol=symbol,
                               x=x, y=y, status=row[key + '_status'],
                               location='missing' if x is None else
                               'in image' if 0 <= x < WIDTH and 0 <= y < HEIGHT else 'outside image'))
        events.append(dict(index=i, review_index=int(row['review_index']),
                           frame_id=row['frame_id'], event_id=row['event_id'],
                           image=image.as_posix(), sync_status=row['sync_status'],
                           frame_pts_s=row['frame_pts_s'],
                           uncertainty_px=finite(row['uncertainty_px']), points=points))
    if not events or len({e['event_id'] for e in events}) != len(events):
        raise ValueError('Expected nonempty, unique event IDs')
    if len({e['frame_id'] for e in events}) != len(events):
        raise ValueError('Expected unique frame IDs')
    return dict(device=device, width=WIDTH, height=HEIGHT, events=events)


def figure_spec(payload):
    events = payload['events']
    traces = []
    # Fixed trace order: background line, three populations, image border,
    # selected line, three selected markers, selected Terra uncertainty circle.
    lines_x, lines_y = [], []
    for event in events:
        points = [p for p in event['points'] if p['x'] is not None]
        lines_x.extend([p['x'] for p in points] + [None])
        lines_y.extend([p['y'] for p in points] + [None])
    traces.append(dict(type='scatter', mode='lines', x=lines_x, y=lines_y,
                       line=dict(color='#64748b', width=.7), opacity=.14,
                       showlegend=False, hoverinfo='skip'))
    for j, (_, label, color, symbol) in enumerate(METHODS):
        points = [e['points'][j] for e in events]
        traces.append(dict(type='scatter', mode='markers', name=label,
            x=[p['x'] for p in points], y=[p['y'] for p in points],
            marker=dict(color=color, symbol=symbol, size=8), opacity=.72,
            customdata=[[e['index'], str(e['review_index']).zfill(5), e['frame_id'],
                         p['status'], p['location']] for e, p in zip(events, points)],
            hovertemplate=f'<b>{label}</b> · detection %{{customdata[0]}} · image %{{customdata[1]}}'
                          '<br>%{customdata[2]}<br>x=%{x:.2f}, y=%{y:.2f}'
                          '<br>%{customdata[3]} · %{customdata[4]}<extra></extra>'))
    traces.append(dict(type='scatter', mode='lines', x=[0, WIDTH, WIDTH, 0, 0],
        y=[0, 0, HEIGHT, HEIGHT, 0], line=dict(color='#94a3b8', width=1.5, dash='dot'),
        showlegend=False, hoverinfo='skip'))
    traces.append(dict(type='scatter', mode='lines', x=[], y=[],
        line=dict(color='#172839', width=4), showlegend=False, hoverinfo='skip'))
    for j, (_, label, color, symbol) in enumerate(METHODS):
        traces.append(dict(type='scatter', mode='markers+text', name=label + ' selected',
            x=[], y=[], text=[label], textposition=('top left', 'bottom center', 'top right')[j],
            textfont=dict(size=15, color='#172839'), cliponaxis=True,
            marker=dict(color=color, symbol=symbol, size=20,
                        line=dict(color='#172839', width=2)),
            showlegend=False, hovertemplate=label + '<br>x=%{x:.2f}, y=%{y:.2f}<extra></extra>'))
    traces.append(dict(type='scatter', mode='lines', x=[], y=[],
        line=dict(color=METHODS[0][2], width=2), showlegend=False, hoverinfo='skip'))
    n = len(events)
    tick_every = max(1, math.ceil(n / 12))
    steps = [dict(label=str(i) if i % tick_every == 0 or i == n else '',
                  value=str(i), method='skip', args=[i]) for i in range(n + 1)]
    layout = dict(autosize=True, height=680,
        margin=dict(l=65, r=25, t=20, b=110), paper_bgcolor='white', plot_bgcolor='white',
        xaxis=dict(title='Image x (pixels)', range=[0, WIDTH], zeroline=False),
        yaxis=dict(title='Image y (pixels)', range=[HEIGHT, 0], scaleanchor='x',
                   scaleratio=1, zeroline=False),
        legend=dict(orientation='h', x=0, y=1.09),
        hovermode='closest', dragmode='pan', uirevision='clicker-review',
        sliders=[dict(active=0, pad=dict(t=45), currentvalue=dict(visible=False),
                      steps=steps, ticklen=3, minorticklen=0)],
        font=dict(family='system-ui, sans-serif', color='#172839'))
    return dict(data=traces, layout=layout)


SCRIPT = r'''
const payload = JSON.parse(document.getElementById('review-data').textContent);
const original = JSON.parse(document.getElementById('figure-data').textContent);
const gd = document.getElementById('plot');
let selected = 0, busy = false, pending = null;
const byId = id => document.getElementById(id);
function bounds(points) {
  const valid = points.filter(p => Number.isFinite(p.x) && Number.isFinite(p.y));
  if (!valid.length) return [[0,payload.width],[payload.height,0]];
  const xs=valid.map(p=>p.x), ys=valid.map(p=>p.y);
  const left=Math.min(...xs), right=Math.max(...xs), top=Math.min(...ys), bottom=Math.max(...ys);
  const px=Math.max(180,(right-left)*.16), py=Math.max(110,(bottom-top)*.18);
  return [[left-px,right+px],[bottom+py,top-py]];
}
function rangeFor(event) {
  const mode=byId('view').value;
  if (mode==='all') return bounds(payload.events.flatMap(e=>e.points));
  if (mode==='follow' && event) return bounds(event.points);
  return [[0,payload.width],[payload.height,0]];
}
async function render(index) {
  const event=index ? payload.events[index-1] : null;
  const data=structuredClone(original.data), layout=structuredClone(original.layout);
  const backdrop=byId('background').checked;
  data[0].opacity=byId('connections').checked ? (index ? .035 : .14) : 0;
  for (let i=1;i<=3;i++) data[i].opacity=index ? (backdrop ? .11 : 0) : .72;
  if (event) {
    const points=event.points.filter(p=>p.x!==null);
    data[5].x=points.map(p=>p.x); data[5].y=points.map(p=>p.y);
    event.points.forEach((p,i)=>{
      data[6+i].x=p.x===null ? [] : [p.x]; data[6+i].y=p.y===null ? [] : [p.y];
      data[6+i].text=[p.label];
    });
    const terra=event.points[0], radius=event.uncertainty_px;
    if (terra.x!==null && radius>0) for (let k=0;k<=60;k++) {
      const a=2*Math.PI*k/60; data[9].x.push(terra.x+radius*Math.cos(a));
      data[9].y.push(terra.y+radius*Math.sin(a));
    }
  }
  const ranges=rangeFor(event); layout.xaxis.range=ranges[0]; layout.yaxis.range=ranges[1];
  layout.sliders[0].active=index;
  layout.uirevision=`${index}-${byId('view').value}`;
  await Plotly.react(gd,data,layout,{responsive:true,scrollZoom:true,displaylogo:false,
    toImageButtonOptions:{format:'png',filename:`${payload.device}_detection_${index}`}});
  selected=index; byId('detection').value=index;
  byId('previous').disabled=index===0; byId('next').disabled=index===payload.events.length;
  byId('selection-label').textContent=event ? `Detection ${index} / ${payload.events.length} · image ${String(event.review_index).padStart(5,'0')}` : `0 · All ${payload.events.length} detections — no highlight`;
  byId('frame-id').textContent=event ? event.frame_id : 'Slide right to highlight a frame’s three points.';
  byId('event-id').textContent=event ? event.event_id : '';
  byId('sync').textContent=event ? `Sync: ${event.sync_status} · source frame PTS ${event.frame_pts_s} s` : '';
  byId('sync').className=event && event.sync_status!=='mapped' ? 'warning' : '';
  byId('preview').hidden=!event; byId('point-table').hidden=!event;
  if (event) {
    byId('frame-link').href=event.image; byId('selected-image').src=event.image;
    byId('selected-image').alt=`Image ${event.review_index}: ${event.frame_id}`;
    event.points.forEach((p,i)=>{
      byId(`point-${i}-x`).textContent=p.x===null ? '—' : p.x.toFixed(2);
      byId(`point-${i}-y`).textContent=p.y===null ? '—' : p.y.toFixed(2);
      byId(`point-${i}-status`).textContent=`${p.location}; ${p.status}`;
    });
  } else byId('selected-image').removeAttribute('src');
  gd.dataset.selection=String(index);
}
function select(index) {
  pending=Math.max(0,Math.min(payload.events.length,Math.round(Number(index)||0)));
  if (busy) return;
  busy=true;
  (async()=>{
    try { while (pending!==null) { const next=pending; pending=null; await render(next); } }
    finally { busy=false; }
  })().catch(error=>{byId('selection-label').textContent=`Plot error: ${error.message}`; console.error(error);});
}
Plotly.newPlot(gd,structuredClone(original.data),structuredClone(original.layout),{responsive:true,scrollZoom:true,displaylogo:false}).then(()=>{
  gd.on('plotly_sliderchange',event=>select(Number(event.step.value)));
  gd.on('plotly_click',event=>{
    const index=event.points?.[0]?.customdata?.[0]; if (Number.isInteger(index)) select(index);
  });
  byId('detection').addEventListener('input',event=>select(event.target.value));
  byId('previous').onclick=()=>select(selected-1); byId('next').onclick=()=>select(selected+1);
  byId('reset').onclick=()=>select(0);
  for (const id of ['view','background','connections']) byId(id).onchange=()=>select(selected);
  window.clickerReview={select,get selected(){return selected;},get busy(){return busy;},payload};
  select(0); document.body.dataset.ready='true';
});
'''


def document(payload, figure):
    def encoded(value):
        return json.dumps(value, allow_nan=False, separators=(',', ':')).replace('</', '<\\/')
    n, device = len(payload['events']), html.escape(payload['device'])
    table = ''.join(f'<tr><th style="border-left:4px solid {color}">{label}</th><td id="point-{i}-x"></td>'
                    f'<td id="point-{i}-y"></td><td id="point-{i}-status"></td></tr>'
                    for i, (_, label, color, _) in enumerate(METHODS))
    return f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{device} · interactive three-way clicker review</title><style>
body{{font:16px system-ui,sans-serif;color:#172839;background:#f4f7fa;margin:0;padding:22px}}
main{{max-width:1800px;margin:auto}}h1{{font-size:25px;margin:0 0 8px}}p{{line-height:1.5}}a{{color:#075a9c}}
.controls{{display:flex;flex-wrap:wrap;align-items:center;gap:14px;background:white;padding:14px;border-radius:10px}}
button,select,input{{font:inherit;padding:7px}}input[type=number]{{width:80px}}input[type=checkbox]{{accent-color:#075a9c}}
.grid{{display:grid;grid-template-columns:minmax(0,1fr) 370px;gap:18px;margin-top:16px}}
#plot,.details{{background:white;border-radius:10px}}.details{{padding:18px;overflow-wrap:anywhere}}
#selection-label{{font-weight:700;font-size:18px}}#frame-id{{font-weight:600}}#selected-image{{width:100%;height:auto;margin-top:14px}}
table{{width:100%;border-collapse:collapse;font-size:13px}}td,th{{padding:8px 4px;border-bottom:1px solid #e2e8f0;text-align:left}}
.warning{{color:#a34c00;font-weight:600}}.hint{{color:#536373;font-size:14px}}[hidden]{{display:none!important}}
@media(max-width:1050px){{.grid{{grid-template-columns:1fr}}.details{{max-width:900px}}}}
</style><script>{get_plotlyjs()}</script></head><body><main>
<h1>Device {device} · Terra / calibrated / original D3</h1>
<p>Slider <b>0</b> shows every detection without a highlight. Slide right to emphasize one frame’s points and connecting line.
Click a background point to select its detection. Scroll to zoom or drag to pan.</p>
<div class="controls"><button id="previous" disabled>← Previous</button><label>Detection <input id="detection" type="number" min="0" max="{n}" value="0" aria-label="Detection index"></label>
<button id="next">Next →</button><button id="reset">0 · No highlight</button>
<label>View <select id="view"><option value="follow">Follow highlighted set</option><option value="image">Image bounds</option><option value="all">All estimates</option></select></label>
<label><input id="background" type="checkbox" checked>Other detections</label>
<label><input id="connections" type="checkbox">Background lines</label></div>
<div class="grid"><div id="plot"></div><aside class="details"><p id="selection-label">Loading plot…</p><p id="frame-id"></p><p id="event-id" class="hint"></p><p id="sync" class="hint"></p>
<table id="point-table" hidden><thead><tr><th>Method</th><th>x</th><th>y</th><th>Status</th></tr></thead><tbody>{table}</tbody></table>
<div id="preview" hidden><a id="frame-link" target="_blank" rel="noopener"><img id="selected-image" alt="Selected overlay frame"></a><p class="hint">Click the image for the original-resolution overlay.</p></div>
<p class="hint">Amber circle: original Terra annotation. Red X: calibrated acoustic estimate. Cyan diamond: original D3 estimate. Missing points remain missing; off-image estimates retain their actual coordinates. Follow view fits the highlighted points.</p>
<p class="hint">Descriptive comparison using deployment fits, including fitting data. Image coordinates use x right / y down.</p>
<p><a href="index.html">Device results</a> · <a href="positions.csv">All positions (CSV)</a> · <a href="../index.html">A/B comparison</a></p></aside></div>
</main><script id="review-data" type="application/json">{encoded(payload)}</script>
<script id="figure-data" type="application/json">{encoded(figure)}</script><script>{SCRIPT}</script></body></html>'''


def export_interactive(directory):
    directory = Path(directory).resolve()
    csv_path = directory / 'positions.csv'
    with csv_path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    manifest_path = directory / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    csv_hash = digest(csv_path)
    if manifest['outputs']['positions.csv'] != csv_hash:
        raise ValueError('Positions CSV does not match the completed review manifest')
    payload = review_payload(rows, manifest['device'])
    if len(payload['events']) != manifest['frame_count']:
        raise ValueError('Frame count mismatch')
    for event in payload['events']:
        if not (directory / event['image']).is_file():
            raise FileNotFoundError(event['image'])
    target = directory / 'positions_interactive.html'
    temporary = target.with_suffix('.html.tmp')
    temporary.write_text(document(payload, figure_spec(payload)))
    temporary.replace(target)
    metadata = dict(schema_version='clicker-interactive-review/1', device=payload['device'],
                    detections=len(payload['events']), slider_min=0,
                    slider_max=len(payload['events']), slider_zero='no highlight',
                    positions_sha256=csv_hash, html_sha256=digest(target),
                    generator_sha256=digest(__file__), plotly_version=plotly.__version__,
                    offline=True, image_links='relative to this directory')
    (directory / 'interactive_review.json').write_text(json.dumps(metadata, indent=2) + '\n')
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('device_directory', type=Path)
    print(json.dumps(export_interactive(parser.parse_args().device_directory), indent=2))
