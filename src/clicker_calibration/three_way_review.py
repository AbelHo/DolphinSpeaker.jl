#!/usr/bin/env python3
"""Render selected clicker frames using sealed annotations and two event predictions."""
import argparse
import collections
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clicker-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'calibration_runs/clicker_20260209'
W, H = 2704, 1520
COLORS = {'terra': '#ffbf32', 'calibrated': '#ff4545', 'old': '#20e5f3'}
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'

def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def unique(items, key):
    result = {}
    for row in items:
        k = row[key]
        if k in result: raise ValueError(f'duplicate {key}: {k}')
        result[k] = row
    return result

def number(value):
    try:
        v = float(value)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError): return None

def location(x, y):
    x, y = number(x), number(y)
    if x is None or y is None: return 'missing'
    return 'in image' if 0 <= x < W and 0 <= y < H else 'outside image'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()

def render_frame(row, dst):
    with Image.open(row['raw_image']) as src:
        if src.size != (W, H): raise ValueError(f"Unexpected image size: {row['raw_image']} {src.size}")
        image = src.convert('RGB')
    draw = ImageDraw.Draw(image)
    for method in ('terra', 'calibrated', 'old'):
        if row[method + '_location'] != 'in image': continue
        x, y = float(row[method + '_x']), float(row[method + '_y'])
        color = COLORS[method]
        # Black under-stroke preserves visibility on light and dark footage.
        if method == 'terra':
            radius = number(row['uncertainty_px'])
            if radius is not None and radius > 0:
                draw.ellipse((x-radius,y-radius,x+radius,y+radius), outline='black', width=8)
                draw.ellipse((x-radius,y-radius,x+radius,y+radius), outline=color, width=4)
            for width, c in ((9,'black'),(4,color)):
                draw.ellipse((x-12,y-12,x+12,y+12),outline=c,width=width)
        elif method == 'calibrated':
            for width, c in ((10,'black'),(5,color)):
                draw.line((x-18,y-18,x+18,y+18), fill=c,width=width)
                draw.line((x-18,y+18,x+18,y-18), fill=c,width=width)
        else:
            for width, c in ((10,'black'),(5,color)):
                draw.polygon([(x,y-23),(x+23,y),(x,y+23),(x-23,y)],outline=c,width=width)
    # Keep caption away from all image markers: choose the least occupied corner.
    boxes = [(12,12,1300,260), (1404,12,2692,260), (12,1260,1300,1508), (1404,1260,2692,1508)]
    points = [(number(row[m+'_x']),number(row[m+'_y'])) for m in COLORS if row[m+'_location']=='in image']
    def score(box):
        return sum(box[0]-45 <= x <= box[2]+45 and box[1]-45 <= y <= box[3]+45 for x,y in points)
    box = min(boxes,key=score)
    overlay = Image.new('RGBA',image.size); d = ImageDraw.Draw(overlay)
    d.rounded_rectangle(box, radius=12, fill=(0,0,0,220))
    f = ImageFont.truetype(FONT, 29); bold = ImageFont.truetype(BOLD, 31)
    x0,y0=box[0]+20,box[1]+12
    d.text((x0,y0), f"{row['review_index']:05d}  {row['frame_id']}", font=bold, fill='white')
    for i,(method,label) in enumerate((('terra','Terra circle'),('calibrated','Calibrated X'),('old','Old D3 diamond'))):
        status = row[method+'_status']; loc = row[method+'_location']
        if method == 'old': status += '; physically feasible=' + row['old_physical_feasible']
        d.text((x0,y0+47+i*42), f'{label}: {loc}; {status}'[:86], font=f, fill=COLORS[method])
    d.text((x0,y0+177), f"Sync: {row['sync_status']} | {row['frame_pts_s']} s", font=f, fill='white')
    Image.alpha_composite(image.convert('RGBA'),overlay).convert('RGB').save(dst,compress_level=3)

def plot(data, destination, full=False):
    fig, ax = plt.subplots(figsize=(13,8))
    labels={'terra':'Terra','calibrated':'Calibrated','old':'Old D3'}
    markers={'terra':'o','calibrated':'x','old':'D'}
    counts={}
    for method in COLORS:
        counts[method]=dict(collections.Counter(r[method+'_location'] for r in data))
        selected=[r for r in data if r[method+'_location'] != 'missing' and (full or r[method+'_location']=='in image')]
        ax.scatter([float(r[method+'_x']) for r in selected], [float(r[method+'_y']) for r in selected],
                   c=COLORS[method], marker=markers[method], s=28, alpha=.65, label=f"{labels[method]}: in {counts[method].get('in image',0)}, outside {counts[method].get('outside image',0)}, missing {counts[method].get('missing',0)}")
    if not full: ax.set_xlim(0,W); ax.set_ylim(H,0)
    else: ax.invert_yaxis(); ax.axvline(0,color='gray'); ax.axvline(W,color='gray'); ax.axhline(0,color='gray'); ax.axhline(H,color='gray')
    ax.set_aspect('equal',adjustable='box'); ax.set_xlabel('Image x (pixels)'); ax.set_ylabel('Image y (pixels)')
    ax.set_title('Click-frame positions — descriptive in-sample deployment model' + (' (full extent)' if full else ''))
    ax.grid(alpha=.2); ax.legend(loc='upper center',bbox_to_anchor=(.5,-.12), ncol=1, fontsize=9)
    fig.tight_layout(); fig.savefig(destination.with_suffix('.png'),dpi=160); fig.savefig(destination.with_suffix('.svg')); plt.close(fig)
    return counts

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('device',choices=['A','B']); ap.add_argument('destination',type=Path); args=ap.parse_args()
    d=args.device; source=RUN/'20261030_v2'/d; model=(RUN/'20261030_v2_fit/A/refinement_v1' if d=='A' else source/'refinement_v1'); legacy=RUN/'20261031_three_way_v1'/d/'legacy_predictions.csv'
    selection=json.loads((source/'human_review_selection.json').read_text()); seal=json.loads((source/'annotations.complete.json').read_text())
    for key,path in [('frames_sha256',source/'frames.csv'),('annotations_sha256',source/'annotations.csv'),('human_review_selection_sha256',source/'human_review_selection.json')]:
        if seal[key] != digest(path): raise ValueError(f'Seal mismatch: {path}')
    frames=rows(source/'frames.csv'); annotations=unique(rows(source/'annotations.csv'),'frame_id'); events=unique(rows(model/'inference_events.csv'),'frame_id'); old=unique(rows(legacy),'event_id')
    ids=selection['accepted_frame_ids']; accepted=set(ids)
    if len(ids)!=len(accepted): raise ValueError('Duplicate selected frame')
    ordered=[(i,r) for i,r in enumerate(frames) if r['frame_id'] in accepted]
    if [r['frame_id'] for _,r in ordered] != ids or len(ordered)!=len(ids): raise ValueError('Selected gallery order mismatch')
    if [i for i,_ in ordered] != list(range(selection['first_index'],selection['last_index']+1)): raise ValueError('Review indices mismatch')
    prepared=json.loads((source/'prepare_manifest.json').read_text())
    raw_hashes={str(Path(r['path']).resolve()):r['sha256'] for r in prepared['outputs']}
    data=[]
    for index,frame in ordered:
        fid=frame['frame_id']; a=annotations[fid]; e=events[fid]; o=old[e['event_id']]
        if o['frame_id'] != fid: raise ValueError('Legacy event/frame identity mismatch')
        raw=Path(frame['image_path']).resolve()
        if raw_hashes.get(str(raw)) != digest(raw): raise ValueError(f'Raw frame checksum mismatch: {raw}')
        row=dict(review_index=index, frame_id=fid,event_id=e['event_id'],source_video=frame['video_file'],source_frame_index=frame['frame_index'],frame_pts_s=frame['frame_pts_s'],raw_image=frame['image_path'],sync_status=e['sync_status'],sync_error_s=e['sync_error_s'],eligibility_status=e['eligibility_status'],terra_x=a['x'],terra_y=a['y'],terra_status=a['visibility'],medium=a['medium'],uncertainty_px=a['uncertainty_px'],annotation_note=a['note'],calibrated_x=e['final_x'],calibrated_y=e['final_y'],calibrated_status=e['final_status'],calibrated_valid=e['final_valid'],calibrated_estimator_status=e['estimator_status'])
        row.update({k:o[k] for k in o if k.startswith('old_')})
        for m in COLORS: row[m+'_location']=location(row[m+'_x'],row[m+'_y'])
        row['image']=f"frames/{index:05d}_{fid}.png"; data.append(row)
    dest=args.destination
    dest.mkdir(parents=True,exist_ok=False); (dest/'frames').mkdir()
    with (dest/'positions.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0])); writer.writeheader(); writer.writerows(data)
    for row in data: render_frame(row,dest/row['image'])
    counts=plot(data,dest/'positions_2d'); plot(data,dest/'positions_2d_full_extent',True)
    with tempfile.TemporaryDirectory(prefix='clicker-video-') as t:
        for i,row in enumerate(data): (Path(t)/f'{i:05d}.png').symlink_to((dest/row['image']).resolve())
        qualification=['ffmpeg','-hide_banner','-loglevel','error','-f','lavfi','-i',f'color=size={W}x{H}:rate=4','-frames:v','2','-c:v','h264_nvenc','-pix_fmt','yuv420p','-f','null','-']
        qualified=subprocess.run(qualification,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE).returncode == 0
        encoder='h264_nvenc' if qualified else 'libx264'
        command=['ffmpeg','-hide_banner','-loglevel','error','-y','-framerate','4','-i',f'{t}/%05d.png','-c:v',encoder,'-pix_fmt','yuv420p','-movflags','+faststart',str(dest/'comparison_4fps.mp4')]
        subprocess.run(command,check=True)
    probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-count_frames','-select_streams','v:0','-show_entries','stream=width,height,r_frame_rate,nb_read_frames,duration,codec_name,pix_fmt','-of','json',str(dest/'comparison_4fps.mp4')]))['streams'][0]
    if probe['r_frame_rate']!='4/1' or int(probe['nb_read_frames'])!=len(data) or (probe['width'],probe['height'])!=(W,H) or abs(float(probe['duration'])-len(data)/4)>0.001: raise ValueError(f'Video verification failed: {probe}')
    manifest=dict(device=d,frame_count=len(data),order='original review gallery order',review_index_first=data[0]['review_index'],review_index_last=data[-1]['review_index'],counts=counts,video=probe,encoding_command=command,description='Descriptive comparison using the in-sample deployment model. Silent 4 fps montage, not real-time or audio-synchronized. Old solver may report finite estimates for physically impossible delays.',input_sha256={str(p):digest(p) for p in [source/'frames.csv',source/'annotations.csv',source/'annotations.complete.json',source/'human_review_selection.json',model/'inference_events.csv',legacy]},outputs={'positions.csv':digest(dest/'positions.csv'),'positions_2d.png':digest(dest/'positions_2d.png'),'comparison_4fps.mp4':digest(dest/'comparison_4fps.mp4')})
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (dest/'index.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>{d} three-way clicker review</title><style>body{{font:18px system-ui;max-width:1100px;margin:40px auto;line-height:1.5}} img{{max-width:100%}}</style><h1>{d} three-way clicker review</h1><p>{html.escape(manifest['description'])}</p><p>Coordinates outside the image are retained in the CSV and full extent plot. Frame captions show missing and outside positions and recorded status.</p><ul><li><a href="comparison_4fps.mp4">Silent 4 fps montage</a></li><li><a href="positions.csv">Positions and statuses</a></li><li><a href="positions_2d_full_extent.png">Full extent plot</a></li><li><a href="manifest.json">Manifest</a></li></ul><img src="positions_2d.png" alt="Position comparison"><h2>Frames</h2>'''+''.join(f'<p><a href="{html.escape(r["image"])}">{r["review_index"]:05d} {html.escape(r["frame_id"])}</a></p>' for r in data))
    print(json.dumps(manifest,indent=2))

if __name__=='__main__': main()
