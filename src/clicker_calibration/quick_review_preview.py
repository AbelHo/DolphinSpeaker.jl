"""Bounded, downsampled Terra-only event excerpts; scientific inputs unchanged."""
from __future__ import annotations
import argparse
import bisect
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
import io
from pathlib import Path
import struct
import time
import wave

import numpy as np
from PIL import Image


def _exporter():
    spec=importlib.util.spec_from_file_location("_quick_preview_exporter",Path(__file__).with_name("video_export.py"))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def complete_float_wav(path):
    """Return the data span only for a complete declared IEEE mono48kHz WAV."""
    p=Path(path)
    if not p.is_file() or p.stat().st_size<56:return None
    with p.open("rb") as f:
        header=f.read(12)
        if header[:4]!=b"RIFF" or header[8:12]!=b"WAVE" or struct.unpack_from("<I",header,4)[0]+8!=p.stat().st_size:return None
        fmt=None;data=None
        while f.tell()+8<=p.stat().st_size:
            tag,n=struct.unpack("<4sI",f.read(8));at=f.tell()
            if at+n>p.stat().st_size:return None
            if tag==b"fmt ":fmt=struct.unpack("<HHIIHH",f.read(16))
            if tag==b"data":data=(at,n)
            f.seek(at+n+(n%2))
        if fmt is None or fmt[0:3]!=(3,1,48000) or fmt[5]!=32 or data is None or data[1]%4:return None
        return data


def _bounded_audio(metadata,global_start,duration,timing,work,ve,ffmpeg,threads):
    fs=int(metadata["wav_files"][0]["fs"])
    # The origin is on the same recorder48kHz sample grid as the full exporter.
    # 4096 native samples on both sides absorb bounded resampler edge support.
    first=math.floor(((global_start-timing["offset_s"])/timing["scale"])*48000)*4-4096
    stop=math.ceil(((global_start+duration-timing["offset_s"])/timing["scale"])*fs)+4096
    stop=4*math.ceil(stop/4)
    mono=np.zeros(stop-first,dtype=np.float32)
    used=[]
    for info in metadata["wav_files"]:
        lo=max(first,int(info["start_sample"]));hi=min(stop,int(info["start_sample"])+int(info["samples"]))
        if hi<=lo:continue
        with wave.open(str(info["path"]),"rb") as f:
            if (f.getframerate(),f.getnchannels(),f.getsampwidth())!=(192000,5,3):raise ValueError("Expected original five-channel192kHzPCM24")
            f.setpos(lo-int(info["start_sample"]));raw=f.readframes(hi-lo)
        if len(raw)!=(hi-lo)*15:raise ValueError("Partial original PCM read")
        b=np.frombuffer(raw,dtype=np.uint8).reshape(-1,5,3)[:,4,:].astype(np.int32)
        value=b[:,0]|(b[:,1]<<8)|(b[:,2]<<16);value=(value^0x800000)-0x800000
        mono[lo-first:hi-first]=value.astype(np.float64)/8388608.0
        used.append(dict(path=info["path"],first_sample=lo-int(info["start_sample"]),samples=hi-lo))
    raw=work/"source_ch5_192k.f32";raw.write_bytes(mono.astype("<f4").tobytes())
    sampled=work/"source_ch5_48k.f32"
    ve._run([ffmpeg,"-v","error","-nostdin","-y","-threads",threads,"-filter_threads",threads,"-f","f32le","-ar",fs,"-ac","1","-i",raw,
             "-af","aresample=48000:async=0","-c:a","pcm_f32le","-f","f32le",sampled],work,"bounded_resample")
    audio=work/"aligned.wav"
    local_offset=timing["offset_s"]-global_start+timing["scale"]*(first/fs)
    ve.write_affine_audio(sampled,audio,duration,local_offset,timing["scale"])
    span=complete_float_wav(audio)
    if span is None:raise ValueError("Incomplete bounded aligned WAV")
    with audio.open("rb") as f:f.seek(span[0]);values=np.frombuffer(f.read(span[1]),dtype="<f4")
    return audio,values,dict(method="bounded physical5PCM24;4096sample padding; gapless sample-count timeline;FFmpeg bandlimit48k then exporter Float64 affine interpolation",
                             recorder_first_sample=first,source_spans=used,local_affine_offset_s=local_offset)


def _checkpoint_snapshot(source,output,frames,ve):
    """Capture one opened CSV version; atomic source replacement cannot mix rows."""
    destination=output/"annotations.checkpoint.csv"
    if destination.exists():raise FileExistsError("Checkpoint snapshot already exists; use a fresh output directory")
    data=None;captured=None
    for _ in range(3):
        with Path(source).open("rb") as f:
            a=os.fstat(f.fileno());value=f.read();b=os.fstat(f.fileno())
        # Ignore ctime/link-count changes when another agent atomically replaces
        # the pathname: this descriptor still owns the complete opened version.
        if (a.st_dev,a.st_ino,a.st_size,a.st_mtime_ns)==(b.st_dev,b.st_ino,b.st_size,b.st_mtime_ns) and len(value)==a.st_size:
            data=value;captured=a;break
    if data is None:raise ValueError("Checkpoint is being written in place; no consistent snapshot")
    temporary=destination.with_suffix(".partial")
    with temporary.open("xb") as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(temporary,destination)
    reader=csv.DictReader(io.StringIO(data.decode("utf-8-sig"),newline=""))
    required={"frame_id","x","y","visibility","medium","uncertainty_px","note","annotator"}
    if reader.fieldnames is None or not required<=set(reader.fieldnames):raise ValueError("Incomplete checkpoint schema")
    rows=list(reader);ids=[r.get("frame_id") for r in rows]
    if not rows or not all(ids) or len(set(ids))!=len(ids) or not set(ids)<=set(frames):raise ValueError("Checkpoint identities are not a unique subset of prepared frames")
    if any(None in r or any(v is None for v in r.values()) for r in rows):raise ValueError("Partial checkpoint CSV row")
    reviewed={r["frame_id"]:r for r in rows if r["annotator"].strip() and not r["annotator"].strip().lower().startswith("pending")}
    if not reviewed:raise ValueError("Checkpoint contains no individually reviewed rows")
    return reviewed,dict(status="provisional_partial_annotation_checkpoint",snapshot_path=str(destination),
        snapshot_sha256=hashlib.sha256(data).hexdigest(),source_path=str(Path(source).resolve()),
        source_opened_inode=captured.st_ino,source_opened_size=captured.st_size,
        captured_utc=datetime.now(timezone.utc).isoformat(),snapshot_row_count=len(rows),
        inspected_count=len(reviewed),total_prepared_frames=len(frames),
        semantics="Immutable opened-version snapshot; no canonical annotations, seal or user selection created; unreviewed rows never filled")


def create_preview(device_dir,output_dir,*,count=20,ffmpeg="/usr/bin/ffmpeg",aligned_wave=None,checkpoint_preview=False):
    ve=_exporter();device=Path(device_dir).resolve();output=Path(output_dir).resolve();output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    metadata=ve._json(device/"prepare_metadata.json");sync=ve._json(device/"synchronization.json");qc=ve._json(device/"sync_eligibility.json")
    timing=ve._timing(metadata,sync,qc)
    events=ve._rows(device/"events.csv","event_id");frames=ve._rows(device/"frames.csv","frame_id")
    checkpoint=None
    if checkpoint_preview:
        annotations,checkpoint=_checkpoint_snapshot(device/"annotations.pending.csv",output,frames,ve)
        accepted=set(annotations);annotations_sha=checkpoint["snapshot_sha256"];completion_sha=None;selection_sha=None
    else:
        annotations=ve._rows(device/"annotations.csv","frame_id")
        seal=ve._json(device/"annotations.complete.json");selection=ve._json(device/"human_review_selection.json")
        for key,name in [("annotations_sha256","annotations.csv"),("frames_sha256","frames.csv"),("human_review_selection_sha256","human_review_selection.json")]:
            if seal[key]!=ve._sha(device/name):raise ValueError("Stale original human annotation seal")
        accepted=set(selection["accepted_frame_ids"]);annotations_sha=ve._sha(device/"annotations.csv")
        completion_sha=ve._sha(device/"annotations.complete.json");selection_sha=ve._sha(device/"human_review_selection.json")
    if qc["events_sha256_after"]!=ve._sha(device/"events.csv"):raise ValueError("Stale timing eligibility")
    saved=qc["input_binding"]["saved_inputs"]
    for key,name in [("metadata","prepare_metadata.json"),("synchronization","synchronization.json")]:
        if saved[key]["sha256"]!=ve._sha(device/name):raise ValueError("Stale timing inputs")
    candidates=[e for e in events.values() if e.get("frame_id") in accepted and e.get("sync_status")=="mapped"
                and ve._number(annotations[e["frame_id"]].get("x")) is not None and ve._number(annotations[e["frame_id"]].get("y")) is not None]
    if checkpoint_preview:
        candidates=[e for e in candidates if annotations[e["frame_id"]]["medium"]=="underwater"
            and annotations[e["frame_id"]]["visibility"] in {"visible","uncertain"}
            and 0<=ve._number(annotations[e["frame_id"]]["x"])<int(frames[e["frame_id"]]["width"])
            and 0<=ve._number(annotations[e["frame_id"]]["y"])<int(frames[e["frame_id"]]["height"])]
    candidates.sort(key=lambda e:float(e["recording_time_s"]))
    n=min(count,len(candidates));chosen=[candidates[round(i*(len(candidates)-1)/(n-1))] for i in range(n)] if n>1 else candidates[:n]
    if not chosen:raise ValueError("No selected synchronized coordinate events")
    byvideo={v["path"]:(i,v) for i,v in enumerate(metadata["videos"])}
    qualification=[]
    for path in dict.fromkeys(frames[e["frame_id"]]["video_file"] for e in chosen):
        chapter,video=byvideo[path]
        event=next(e for e in chosen if frames[e["frame_id"]]["video_file"]==path)
        frame=frames[event["frame_id"]]
        probe=ve._source_probe_saved(video,str(Path(ffmpeg).with_name("ffprobe")),output,f"source_header_{chapter}")
        qualification.append(dict(info=video,probe=probe,events=[dict(source_pts_s=float(frame["frame_pts_s"]),prepared_image_path=frame["image_path"])]))
    resources=ve._resources(ffmpeg,output,dict(videos=qualification),0,"auto","auto",23,"ultrafast")
    threads=resources["cpu_threads"]
    plan=[];parts=[];audios=[];full_span=complete_float_wav(aligned_wave) if aligned_wave else None
    for number,event in enumerate(chosen):
        work=output/f"excerpt_{number:02d}";work.mkdir(exist_ok=False)
        frame=frames[event["frame_id"]];chapter,video=byvideo[frame["video_file"]];pts=list(map(float,video["pts"]));index=int(frame["frame_index"]);target=pts[index]
        if abs(target-float(frame["frame_pts_s"]))>1e-9:raise ValueError("Prepared target PTS differs")
        first=max(0,bisect.bisect_left(pts,target-0.4));last=min(len(pts),bisect.bisect_left(pts,target+0.6))
        clip_start=pts[first];clip_end=pts[last] if last<len(pts) else float(video["duration_s"])
        duration=clip_end-clip_start;global_start=timing["bases"][chapter]+clip_start-float(video.get("audio_start_s",0))
        marks=ve._event_marks(event,annotations[event["frame_id"]],None,int(video["width"]),int(video["height"]),accepted=True,annotation_only=True)
        marks.update(event_id=event["event_id"],frame_id=event["frame_id"],frame_index=index)
        marks["label"]+=f" | source-global {timing['bases'][chapter]+target:.6f}s"
        ass=work/"overlay.ass";fallback=ve._write_ass(ass,int(video["width"]),int(video["height"]),pts,[marks],annotation_only=True)
        if fallback:raise ValueError("Selected target lacks exact ASS interval")
        text=f"PREVIEW downsampled1352x760 | excerpt{number+1}/{n} | physical5 hydrophone audio | original Terra uncertainty"
        if checkpoint_preview:
            text=f"PROVISIONAL annotation-checkpoint excerpts | {checkpoint['inspected_count']}/{len(frames)} inspected | excerpt{number+1}/{n} | downsampled1352x760"
        with ass.open("a") as f:f.write(f"Dialogue: 1,0:00:00.00,9:59:59.99,Text,,0,0,0,,{{\\an7\\pos(8,110)}}{text}\n")
        rawpng=work/"target_raw.png";part=work/"video.mp4";seek=max(0,clip_start-0.8)
        pixel_format=next(v["probe"]["pix_fmt"] for v in qualification if v["info"]["path"]==video["path"])
        download="hwdownload,format=nv12,format="+pixel_format+"," if resources["decoder"]!="cpu" else ""
        graph=(f"[0:v]{download}settb=1/1000000,trim=start={clip_start:.6f}:end={clip_end:.6f},split=2[raw][base];"
               f"[raw]select='eq(pts,{round(target*1e6)})'[check];"
               f"[base]ass=filename='{ve._filter_path(ass)}',scale=1352:760,setpts=PTS-STARTPTS[v]")
        ve._run([ffmpeg,"-v","error","-nostdin","-y",*ve._decode_options(resources),"-filter_complex_threads",threads,
                 "-ss",f"{seek:.6f}","-t",f"{clip_end-seek+0.2:.6f}","-copyts","-i",video["path"],"-filter_complex",graph,
                 "-map","[v]","-an",*ve._encode_options(resources,23,"ultrafast"),"-pix_fmt","yuv420p",
                 "-fps_mode","passthrough","-enc_time_base","1:1000000","-video_track_timescale","1000000",part,
                 "-map","[check]","-c:v","png","-threads",threads,"-frames:v","1",rawpng],work,"bounded_video")
        with Image.open(rawpng) as im,Image.open(frame["image_path"]) as original:
            a=np.array(im.convert("RGB"));b=np.array(original.convert("RGB"))
        if a.shape!=b.shape or not np.array_equal(a,b):raise ValueError("Bounded decoded target does not equal prepared PNG pixels")
        if full_span:
            f=Path(aligned_wave).open("rb");start_sample=round(global_start*48000);samples=round(duration*48000)
            f.seek(full_span[0]+start_sample*4);values=np.frombuffer(f.read(samples*4),dtype="<f4").copy();f.close()
            if len(values)!=samples:raise ValueError("Aligned reference excerpt truncated")
            audioinfo=dict(method="complete declared aligned full-reference WAV",path=str(aligned_wave))
        else:
            _,values,audioinfo=_bounded_audio(metadata,global_start,duration,timing,work,ve,ffmpeg,threads)
        audios.append(values);parts.append(part)
        plan.append(dict(event_id=event["event_id"],frame_id=event["frame_id"],source_video=video["path"],frame_index=index,
                         source_pts_s=target,source_global_time_s=timing["bases"][chapter]+target,
                         clip_source_start_s=clip_start,clip_source_end_s=clip_end,clip_global_start_s=global_start,
                         target_preview_local_s=target-clip_start,source_frame_count=last-first,
                         exact_prepared_target_pixels=True,prepared_png_sha256=ve._sha(frame["image_path"]),audio=audioinfo))
        ve._write_json(output/"progress.json",dict(status="rendering",completed=number+1,total=n,elapsed_s=time.monotonic()-started))
    values=np.concatenate(audios);peak=float(np.max(np.abs(values)));gain=min(512.0,0.6/peak) if peak else 1.0;values=values.astype(np.float64)*gain
    if np.max(np.abs(values))>0.600001:raise ValueError("Audio preview gain safety bound exceeded")
    audio=output/"preview_audio.wav"
    with wave.open(str(audio),"wb") as f:f.setnchannels(1);f.setsampwidth(2);f.setframerate(48000);f.writeframes(np.rint(values*32767).astype("<i2").tobytes())
    listing=output/"parts.ffconcat";listing.write_text("ffconcat version 1.0\n"+"".join(f"file {ve._quote_file(p)}\n" for p in parts))
    destination=output/f"{device.name}_terra_quick_review.mp4"
    ve._run([ffmpeg,"-v","error","-nostdin","-y","-threads",threads,"-f","concat","-safe","0","-i",listing,"-i",audio,
             "-map","0:v:0","-map","1:a:0","-c:v","copy","-c:a","aac","-b:a","128k","-threads",threads,"-movflags","+faststart",destination],output,"mux_preview")
    receipt=dict(status="complete",video=str(destination),video_sha256=ve._sha(destination),excerpts=plan,event_count=len(plan),
                 duration_audio_s=len(values)/48000,resolution=[1352,760],threads=threads,resources=resources,source_unchanged=True,
                 annotation_only=True,preview="downsampled excerpts; no acoustic model markers; original uncertainty shown in amber",
                 annotations_sha256=annotations_sha,completion_sha256=completion_sha,selection_sha256=selection_sha,
                 checkpoint_preview=bool(checkpoint_preview),annotation_checkpoint=checkpoint,
                 canonical_annotation_files_written=False,sync_eligibility_sha256=ve._sha(device/"sync_eligibility.json"),
                 timing=timing,audio_gain=gain,source_audio_peak=peak,output_audio_peak=float(np.max(np.abs(values))),
                 elapsed_s=time.monotonic()-started,helper_sha256=ve._sha(__file__))
    ve._write_json(output/"preview_manifest.json",receipt);return receipt


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("device_dir",type=Path);parser.add_argument("output_dir",type=Path)
    parser.add_argument("--count",type=int,default=20);parser.add_argument("--aligned-wave",type=Path);parser.add_argument("--checkpoint-preview",action="store_true")
    args=parser.parse_args()
    try:print(json.dumps(create_preview(**vars(args))))
    except BaseException as error:
        _exporter()._log_failure(Path(args.output_dir),"quick_review_preview",error);raise
