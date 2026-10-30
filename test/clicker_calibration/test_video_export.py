"""Short synthetic media controls; never encode real A/B recordings."""
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import wave

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("video_export", ROOT / "src/clicker_calibration/video_export.py")
VE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VE)


def write_csv(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def command(args, folder, label):
    return VE._run(args, folder, label)


def rgb(path, indices=None):
    args = ["/usr/bin/ffmpeg", "-v", "error", "-threads", "0", "-i", str(path)]
    if indices is not None:
        expression = "+".join(f"eq(n,{i})" for i in indices)
        args += ["-vf", f"select='{expression}'", "-fps_mode", "passthrough"]
    data = subprocess.check_output(args + ["-threads", "0", "-f", "rawvideo", "-pix_fmt", "rgb24", "pipe:1"])
    return np.frombuffer(data, np.uint8).reshape(-1,120,160,3)


def colored(image, x, y, channel):
    region = image[y-2:y+3,x-2:x+3].astype(int)
    others = [i for i in range(3) if i != channel]
    return bool(((region[:,:,channel] > 150) & (region[:,:,others[0]] < 110) & (region[:,:,others[1]] < 110)).any())


def fixture(root):
    device=root/"device";model=root/"model";output=root/"export"
    device.mkdir();model.mkdir()
    fs=192000;offset=0.025;scale=1.003
    lengths=[int(.19*fs),int(.41*fs)];total=sum(lengths)
    pulse_times=[(.11-offset)/scale,(.263333-offset)/scale]
    samples=np.arange(total)
    waveform=sum(.8*np.exp(-((samples-t*fs)/6)**2) for t in pulse_times)
    wavs=[];position=0
    for i,n in enumerate(lengths):
        path=root/f"recorder{i}.wav"
        data=np.zeros((n,5),dtype="<i2")
        data[:,4]=np.round(waveform[position:position+n]*30000).astype("<i2")
        # A different channel has a distractor; physical5 must remain selected.
        if i==0:data[int(.14*fs),0]=30000
        with wave.open(str(path),"wb") as stream:
            stream.setnchannels(5);stream.setsampwidth(2);stream.setframerate(fs);stream.writeframes(data.tobytes())
        wavs.append(dict(path=str(path),fs=fs,channels=5,samples=n,start_sample=position));position+=n
    videos=[]
    for i,(count,duration) in enumerate(zip([6,6,9],[.11,.12,.1])):
        raw=root/f"source_video{i}.mp4";path=root/f"chapter{i}.mp4";sound=root/f"camera{i}.wav"
        command(["/usr/bin/ffmpeg","-v","error","-nostdin","-y","-f","lavfi","-i","color=c=gray:s=160x120:r=60000/1001",
                 "-frames:v",count,"-an","-c:v","libx264","-threads","0","-pix_fmt","yuv420p","-video_track_timescale","1000000",raw],root,f"make_video{i}")
        with wave.open(str(sound),"wb") as stream:
            stream.setnchannels(1);stream.setsampwidth(2);stream.setframerate(48000);stream.writeframes(bytes(round(duration*48000)*2))
        command(["/usr/bin/ffmpeg","-v","error","-nostdin","-y","-threads","0","-i",raw,"-i",sound,
                 "-c:v","copy","-c:a","aac","-threads","0","-video_track_timescale","1000000",path],root,f"mux_source{i}")
        probe=VE._probe(path,"/usr/bin/ffprobe",root,f"fixture_probe{i}")
        videos.append(dict(path=str(path),pts=probe["pts"],width=160,height=120,duration_s=probe["duration"],
                           audio_duration_s=duration,audio_start_s=0,video_start_s=0,
                           global_start_s=[0,.112,.234][i]))
    metadata=dict(wav_files=wavs,videos=videos)
    VE._write_json(device/"prepare_metadata.json",metadata)
    VE._write_json(device/"synchronization.json",dict(status="valid",offset_s=.02,scale=1.002))
    positions=[(0,2,35),(1,0,70),(2,2,35),(2,5,70),(2,8,105),(2,7,120)]
    original_image=root/"original.png"
    command(["/usr/bin/ffmpeg","-v","error","-y","-threads","0","-i",videos[0]["path"],
             "-frames:v","1","-threads","0",original_image],root,"original_png")
    frames=[];events=[];annotations=[];predictions=[]
    bases=[0,.11,.23]
    for i,(chapter,index,x) in enumerate(positions):
        frame=f"frame{i}";event=f"event{i}";pts=videos[chapter]["pts"][index]
        frames.append(dict(frame_id=frame,video_file=videos[chapter]["path"],frame_index=index,frame_pts_s=pts,
                           video_global_time_s=bases[chapter]+pts,image_path=str(original_image),width=160,height=120))
        events.append(dict(event_id=event,frame_id=frame,video_file=videos[chapter]["path"],frame_index=index,
                           recording_time_s=(bases[chapter]+pts-offset)/scale,sync_status="mapped"))
        annotations.append(dict(frame_id=frame,x=x,y=90,visibility="excluded" if i==4 else "uncertain" if i==3 else "visible",
                                medium="underwater",uncertainty_px=35 if i==3 else 2,note="synthetic",annotator="synthetic"))
        if i!=5:
            predictions.append(dict(event_id=event,predicted_x=300 if i==2 else 105 if i==0 else 130,
                                    predicted_y=90,predicted_valid="false" if i==3 else "true",
                                    predicted_status="impossible_tdoa" if i==3 else "outside_image" if i==2 else "ok"))
    write_csv(device/"frames.csv",frames);write_csv(device/"events.csv",events);write_csv(device/"annotations.csv",annotations)
    write_csv(model/"predictions.csv",predictions)
    VE._write_json(device/"human_review_selection.json",dict(schema_version="clicker-human-review-selection-v1",
        frames_sha256=VE._sha(device/"frames.csv"),accepted_frame_ids=[row["frame_id"] for row in annotations if row["visibility"]!="excluded"]))
    VE._write_json(device/"annotations.complete.json",dict(schema_version="clicker-annotations-v1",n_frames=len(frames),
        frames_sha256=VE._sha(device/"frames.csv"),annotations_sha256=VE._sha(device/"annotations.csv"),
        human_review_selection_sha256=VE._sha(device/"human_review_selection.json")))
    (device/"event_snippets.f64").write_bytes(b"synthetic cached fixture, never used for fitting")
    sources=[Path(x["path"]) for x in wavs+videos]
    prepared=["events.csv","frames.csv","synchronization.json","prepare_metadata.json","event_snippets.f64"]
    VE._write_json(device/"prepare_manifest.json",dict(status="complete",frame_count=len(frames),event_count=len(events),
        events_csv=str(device/"events.csv"),frames_csv=str(device/"frames.csv"),
        snippet_cache=dict(path=str(device/"event_snippets.f64"),sha256=VE._sha(device/"event_snippets.f64"),bytes=(device/"event_snippets.f64").stat().st_size),
        outputs=[dict(path=str(device/name),sha256=VE._sha(device/name),bytes=(device/name).stat().st_size) for name in prepared],
        request=dict(inputs=[dict(path=str(path),sha256=VE._sha(path)) for path in sources])))
    VE._write_json(device/"sync_eligibility.json",dict(events_sha256_after=VE._sha(device/"events.csv"),input_binding=dict(saved_inputs={
        "metadata":dict(sha256=VE._sha(device/"prepare_metadata.json")),"synchronization":dict(sha256=VE._sha(device/"synchronization.json"))}),review=dict(
        counterfactual_fit=dict(offset_s=offset,scale=scale),
        chapter_bases=[dict(video_file=v["path"],exact_declared_start_s=base) for v,base in zip(videos,bases)],
        events=[dict(event_id=e["event_id"],final_sync_status="mapped") for e in events])))
    bindings={name:VE._sha(device/name) for name in ["annotations.csv","annotations.complete.json","frames.csv","events.csv"]}
    VE._write_json(model/"calibration.json",dict(hashes=bindings))
    VE._write_json(model/"report.json",dict(outputs_sha256={name:VE._sha(model/name) for name in ["calibration.json","predictions.csv"]}))
    return device,model,output


class VideoExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent=ROOT/"calibration_runs/clicker_20260209/20261030_v1/test_scratch"
        parent.mkdir(parents=True,exist_ok=True)
        cls.root=Path(tempfile.mkdtemp(prefix="video_export_",dir=parent))
        cls.device,cls.model,cls.output=fixture(cls.root)
        with patch.dict(os.environ,{"CLICKER_ERROR_LOG":str(cls.root/"expected_errors.jsonl")}):
            cls.result=VE.export_video(cls.device,cls.model,cls.output,threads=0)

    def test_auto_resources_nvenc_mapping_and_reference_gain(self):
        resources=self.result["resources"]
        self.assertEqual(resources["cpu_threads"],VE._cpu_count())
        self.assertGreater(resources["cpu_threads"],2)
        self.assertIn(resources["encoder"],["h264_nvenc","libx264"])
        if resources["encoder"]=="h264_nvenc":
            flags=resources["encoder_options"]
            self.assertIn("-cq",flags);self.assertIn("p4",flags);self.assertNotIn("-crf",flags)
        audio=self.result["review_audio_gain"]
        self.assertTrue(audio["reference_WAV_unscaled"])
        self.assertAlmostEqual(audio["constant_gain"],min(512,0.6/audio["whole_reference_peak"]))
        self.assertLessEqual(audio["predicted_MP4_PCM_peak"],0.600001)
        self.assertTrue(audio["no_peak_clipping"])
        self.assertLess(audio["aac_decoded_peak"],1.0)
        for chapter in VE._json(self.output/self.result["generation"]/"marker_plan.json")["chapters"]:
            self.assertTrue(chapter["events"])

    def test_saved_pts_probe_uses_headers_not_full_source_decode(self):
        metadata=VE._json(self.device/"prepare_metadata.json")
        info=metadata["videos"][0]
        header=VE._source_probe_saved(info,"/usr/bin/ffprobe",self.root,"saved_header_control")
        self.assertEqual(header["pts"],info["pts"])
        self.assertIn("freshly verified source media SHA256",header["pts_origin"])
        text=(self.root/"saved_header_control.log").read_text()
        self.assertNotIn("-show_frames",text)
        self.assertNotIn("-threads",text)

    def test_frame_count_pts_target_and_previous(self):
        r=self.result
        self.assertEqual(r["source_frames"],21);self.assertEqual(r["output_frames"],21)
        self.assertEqual(r["chapter_bases_s"],[0,.11,.23]);self.assertLessEqual(r["maximum_pts_error_s"],2.1e-6)
        # Last chapter video extends past its declared audio: all nine frames survive.
        images=rgb(r["video_path"],[1,2,5,6,20])
        self.assertFalse(colored(images[0],35,90,1));self.assertTrue(colored(images[1],35,90,1))
        self.assertTrue(colored(images[1],105,90,0))
        self.assertFalse(colored(images[2],70,90,1));self.assertTrue(colored(images[3],70,90,1))
        self.assertEqual(len(images),5)
        frames=VE._json(self.output/r["generation"]/"frame_manifest.json")["events"]
        self.assertEqual(len(frames),6)
        first=next(row for row in frames if row["event_id"]=="event0")
        raw=rgb(self.output/r["generation"]/first["raw"])[0]
        overlay=rgb(self.output/r["generation"]/first["overlay"])[0]
        self.assertFalse(colored(raw,35,90,1));self.assertTrue(colored(overlay,35,90,1))
        statuses={row["event_id"]:row["acoustic_status"] for row in frames}
        self.assertEqual(statuses["event2"],"outside_image")
        self.assertEqual(statuses["event3"],"impossible_tdoa")
        self.assertEqual(statuses["event4"],"user_excluded")
        self.assertEqual(statuses["event5"],"missing_prediction")
        uncertain=next(row for row in frames if row["event_id"]=="event3")
        amber=rgb(self.output/r["generation"]/uncertain["overlay"])[0,88:93,68:73].astype(int)
        self.assertTrue(((amber[:,:,0]>180)&(amber[:,:,1]>120)&(amber[:,:,2]<100)).any())
        outside=next(row for row in frames if row["event_id"]=="event2")
        pixels=rgb(self.output/r["generation"]/outside["overlay"])[0,80:105].astype(int)
        self.assertFalse(((pixels[:,:,0]>150)&(pixels[:,:,1]<100)&(pixels[:,:,2]<100)).any())

    def test_affine_physical5_lossless_review_audio(self):
        path=Path(self.result["reference_audio_path"]);data=path.read_bytes()
        self.assertEqual(struct.unpack_from("<H",data,20)[0],3)  # IEEE Float32 PCM, not AAC or quantized rate.
        self.assertEqual(data[36:40],b"fact")
        audio=np.frombuffer(data, dtype="<f4",offset=56)
        self.assertEqual(len(audio),self.result["review_audio"]["samples"])
        self.assertEqual(len(audio),math.ceil(self.result["video_duration_s"]*48000-1e-7))
        for target in [.11,.263333]:
            index=round(target*48000);lo=index-12;hi=index+13
            observed=lo+int(np.argmax(abs(audio[lo:hi])))
            self.assertLessEqual(abs(observed-index),1)
        self.assertLess(float(np.max(abs(audio[:1100]))),1e-6)  # before recorder origin: exact zero.
        distractor=round((.025+1.003*.14)*48000)
        self.assertLess(float(np.max(abs(audio[distractor-10:distractor+11]))),1e-5)
        cache=self.result["source_audio_cache"]
        self.assertLessEqual(abs(cache["samples"]-cache["expected_rounded_samples"]),1)
        self.assertEqual(self.result["timing"]["offset_s"],.025)
        self.assertEqual(self.result["timing"]["scale"],1.003)

    def test_resume_and_partial_output_rejection(self):
        before=VE._sha(self.output/"export_manifest.json")
        r=VE.export_video(self.device,self.model,self.output,threads=0)
        self.assertEqual(r,self.result);self.assertEqual(before,VE._sha(self.output/"export_manifest.json"))
        Path(r["video_path"]).unlink()
        rerun=VE.export_video(self.device,self.model,self.output,threads=0)
        self.assertNotEqual(rerun["generation"],r["generation"])
        self.assertTrue(Path(rerun["video_path"]).is_file());self.assertEqual(rerun["output_frames"],21)

    def test_stale_annotation_model_and_coarse_sync_rejected(self):
        ann=self.device/"annotations.csv";seal=self.device/"annotations.complete.json"
        old,oldseal=ann.read_bytes(),seal.read_bytes()
        try:
            ann.write_bytes(old+b"\n")
            value=VE._json(seal);value["annotations_sha256"]=VE._sha(ann);VE._write_json(seal,value)
            with patch.dict(os.environ,{"CLICKER_ERROR_LOG":str(self.root/"expected_errors.jsonl")}):
                with self.assertRaisesRegex(ValueError,"current annotations"):
                    VE.export_video(self.device,self.model,self.output,threads=0)
        finally:
            ann.write_bytes(old);seal.write_bytes(oldseal)
        metadata=VE._json(self.device/"prepare_metadata.json");sync=VE._json(self.device/"synchronization.json")
        with self.assertRaisesRegex(ValueError,"coarse chapter bases"):
            VE._timing(metadata,sync)

    def test_ass_quantization_and_exact_index_fallback(self):
        self.assertEqual(VE.ass_interval(.016683,.033367),(2,16))
        self.assertIsNone(VE.ass_interval(0,.005))
        event=dict(event_id="fallback",frame_index=1,label="fallback",terra=dict(x=70,y=90,color="00FF00",radius=11),acoustic=None)
        ass=self.root/"fallback.ass";fallback=VE._write_ass(ass,160,120,[0,.005,.010],[event])
        self.assertEqual(len(fallback),1)
        out=self.root/"fallback.mp4"
        graph=f"ass=filename='{VE._filter_path(ass)}',"+",".join(VE._fallback_filters(fallback))
        command(["/usr/bin/ffmpeg","-v","error","-y","-f","lavfi","-i","color=c=gray:s=160x120:r=200",
                 "-frames:v","3","-vf",graph,"-c:v","libx264","-threads","0","-filter_threads","1",out],self.root,"fallback_encode")
        images=rgb(out)
        self.assertFalse(colored(images[0],70,90,1));self.assertTrue(colored(images[1],70,90,1))

    def test_chunked_affine_linear_samples_and_zero_padding(self):
        raw=self.root/"linear.f32";np.array([0,1,2,3],dtype="<f4").tofile(raw)
        output=self.root/"linear.wav"
        VE.write_affine_audio(raw,output,6/48000,.5/48000,1.0,chunk_samples=2)
        values=np.frombuffer(output.read_bytes(),dtype="<f4",offset=56)
        np.testing.assert_array_equal(values,[0,.5,1.5,2.5,0,0])

    def test_annotation_only_without_model_precise_sound_and_uncertainty_circle(self):
        root=self.root/"annotation_only";root.mkdir()
        device,model,output=fixture(root)
        # No fake model/predictions are supplied or manufactured for this mode.
        import shutil
        shutil.rmtree(model)
        result=VE.export_video(device,output_dir=output,annotation_only=True,threads=0)
        self.assertEqual(result["mode"],"annotation_only")
        self.assertEqual(Path(result["video_path"]).name,"terra_review_synchronized.mp4")
        self.assertEqual(result["output_frames"],21)
        self.assertLessEqual(result["maximum_pts_error_s"],2.1e-6)
        self.assertEqual(result["timing"]["offset_s"],.025)
        self.assertEqual(result["timing"]["scale"],1.003)
        self.assertFalse(any("calibration.json" in name or "predictions.csv" in name for name in result["inputs_sha256"]))
        self.assertIsNone(result["annotation_fit_policy_bridge"])
        images=rgb(result["video_path"],[1,2,5,6])
        self.assertFalse(colored(images[0],35,90,1));self.assertTrue(colored(images[1],35,90,1))
        self.assertFalse(colored(images[1],105,90,0))
        self.assertFalse(colored(images[2],70,90,1));self.assertTrue(colored(images[3],70,90,1))
        mapping=VE._json(output/result["generation"]/"frame_manifest.json")["events"]
        uncertain=next(row for row in mapping if row["event_id"]=="event3")
        image=rgb(output/result["generation"]/uncertain["overlay"])[0].astype(int)
        amber=(image[:,:,0]>180)&(image[:,:,1]>100)&(image[:,:,2]<100)
        # The annotation radius is35px, rather than the fixed11px center marker.
        self.assertTrue(amber[88:93,102:108].any())
        self.assertFalse(amber[88:93,111:116].any())
        self.assertTrue(amber[52:59,68:73].any())
        for chapter in VE._json(output/result["generation"]/"marker_plan.json")["chapters"]:
            for event in chapter["events"]:
                self.assertIsNone(event["acoustic"])
                self.assertNotIn("Acoustic",event["label"])
        audio=np.frombuffer(Path(result["reference_audio_path"]).read_bytes(),dtype="<f4",offset=56)
        index=round(.263333*48000);found=index-12+int(np.argmax(abs(audio[index-12:index+13])))
        self.assertLessEqual(abs(found-index),1)
        self.assertEqual(VE.export_video(device,output_dir=output,annotation_only=True,threads=0),result)

    def test_selected_only_events_leave_full_timeline_but_omit_other_markers_and_pngs(self):
        root=self.root/"selected_only";root.mkdir()
        device,model,output=fixture(root)
        selection=VE._json(device/"human_review_selection.json")
        selection.update(publish_selected_only=True,accepted_frame_ids=["frame0","frame3"])
        VE._write_json(device/"human_review_selection.json",selection)
        seal=VE._json(device/"annotations.complete.json")
        seal["human_review_selection_sha256"]=VE._sha(device/"human_review_selection.json")
        VE._write_json(device/"annotations.complete.json",seal)
        result=VE.export_video(device,output_dir=output,annotation_only=True)
        self.assertEqual(result["source_frames"],21);self.assertEqual(result["output_frames"],21)
        self.assertEqual(result["event_frame_pairs"],2)
        self.assertEqual({r["event_id"] for r in result["user_omitted_events"]},{"event1","event2","event4","event5"})
        pairs=VE._json(output/result["generation"]/"frame_manifest.json")["events"]
        self.assertEqual({r["frame_id"] for r in pairs},{"frame0","frame3"})
        self.assertEqual(len(list((output/result["generation"]/"frames/raw").glob("*.png"))),2)
        images=rgb(result["video_path"],[2,6])
        self.assertTrue(colored(images[0],35,90,1))
        self.assertFalse(colored(images[1],70,90,1))
        self.assertLessEqual(result["maximum_pts_error_s"],2.1e-6)
        selection["accepted_frame_ids"]=["frame0"]
        VE._write_json(device/"human_review_selection.json",selection)
        with self.assertRaisesRegex(ValueError,"selection"):
            VE.export_video(device,output_dir=output,annotation_only=True)

    def test_real_policy_validator_bridge_keeps_original_uncertain_amber(self):
        helper=ROOT/"src/clicker_calibration/fit_annotation_policy.py"
        self.assertTrue(helper.is_file(),"The explicit policy helper must be ready before exporter integration")
        specification=importlib.util.spec_from_file_location("test_export_policy_helper",helper)
        policy=importlib.util.module_from_spec(specification);specification.loader.exec_module(policy)
        root=self.root/"policy_bridge";root.mkdir()
        original,old_model,output=fixture(root)
        rows=list(VE._rows(original/"annotations.csv","frame_id").values())
        rows[0]["visibility"]="uncertain"
        write_csv(original/"annotations.csv",rows)
        seal=VE._json(original/"annotations.complete.json");seal["annotations_sha256"]=VE._sha(original/"annotations.csv")
        VE._write_json(original/"annotations.complete.json",seal)
        original_bytes=(original/"annotations.csv").read_bytes()
        fit=root/"fit_input"
        qualified=policy.create_fit_annotation_input(original,fit,policy_name="user_confirmed_uncertain_coordinates",
                                                    reason="Explicit synthetic control of approved uncertain-coordinate bridge")
        self.assertTrue(qualified["valid"])
        model=fit/"refinement_v1";model.mkdir()
        (model/"predictions.csv").write_bytes((old_model/"predictions.csv").read_bytes())
        VE._write_json(model/"calibration.json",dict(hashes={name:VE._sha(fit/name) for name in
            ["annotations.csv","annotations.complete.json","events.csv","frames.csv"]}))
        VE._write_json(model/"report.json",dict(outputs_sha256={name:VE._sha(model/name) for name in ["calibration.json","predictions.csv"]}))
        with patch.dict(os.environ,{"CLICKER_ERROR_LOG":str(root/"expected_errors.jsonl")}):
            result=VE.export_video(original,model,output,threads=0)
            self.assertTrue(result["annotation_fit_policy_bridge"]["valid"])
            self.assertEqual((original/"annotations.csv").read_bytes(),original_bytes)
            self.assertEqual(VE._rows(fit/"annotations.csv","frame_id")["frame0"]["visibility"],"visible")
            mapping=VE._json(output/result["generation"]/"frame_manifest.json")["events"]
            event=next(row for row in mapping if row["event_id"]=="event0")
            self.assertEqual(event["terra_visibility"],"uncertain")
            image=rgb(output/result["generation"]/event["overlay"])[0,88:93,33:38].astype(int)
            self.assertTrue(((image[:,:,0]>180)&(image[:,:,1]>120)&(image[:,:,2]<100)).any())
            # Even resume must revalidate the explicit bridge, not trust a prior
            # successful output as permission to bypass a changed policy.
            path=fit/"annotation_fit_policy.json";saved=path.read_bytes()
            path.write_bytes(saved+b"\n")
            with self.assertRaises(ValueError):VE.export_video(original,model,output,threads=0)


if __name__=="__main__":
    unittest.main()
