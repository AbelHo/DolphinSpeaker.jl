#!/usr/bin/env bash

if [ $# -eq 0 ]
then
  echo 'vid_addspectro.sh [videoFile_path] [acousticFile_path] [output_newVideoFile_path]'
  exit 0
fi


#ffplay -f lavfi "amovie=$1,asplit=3[out1][a][b]; [out1]; [b]showspectrum=s=1440x400[spectrum]; [waves][spectrum] vstack[out0]"

# ffplay -f lavfi "amovie=$1,asplit=3[out1][a][b]; [a]showwaves=s=640x240:mode=cline:draw=full[waves]; [b]showspectrum=s=640x240[spectrum]; [waves][spectrum] vstack[out0]"

# ffplay -f lavfi "amovie=$1,asplit=2[out1][a]; [a]showspectrum=color=channel:scale=cbrt:orientation=vertical:overlap=1:s=2048x1024[out0]"


## works
# ffplay -f lavfi "amovie=$1,asplit=3[out1][a][b]; [a]showwaves=s=640x240:mode=cline:draw=full[waves]; [b]showspectrum=s=640x240[spectrum]; [waves][spectrum] vstack[out0]" 


# ffmpeg -f s16le -ar 400000 -ac 4 -i "$2" -i "$1" -filter_complex "[0:0]asplit[a][b]; [a]showwaves=s=1440x400:mode=cline:colors=blue[waves]; [b]showspectrum=s=1440x400[spectrum]; [waves][spectrum] vstack[outAu]; [1:v][outAu] vstack [out]" -map "[out]" -map 0:0 -acodec pcm_s16le -f matroska pipe:1 | ffplay -


# ffmpeg -f s16le -ar 400000 -ac 4 -i "$2" -i "$1" -filter_complex "[0:0]asplit[a][b]; [a]showwaves=s=1440x400:mode=cline:colors=blue[waves]; [b]showspectrum=s=1440x400[spectrum]; [waves][spectrum] vstack[outAu]; [1:v][outAu] vstack [out]" -map "[out]" -map 0:0 -acodec pcm_s16le -f matroska "$3"

ffmpeg -f s16le -ar 400000 -ac 5 -i "$2" -i "$1" -filter_complex "[0:0]asplit[a][b]; [a]showwaves=s=5760x400:mode=cline:colors=blue[waves]; [b]showspectrum=s=5760x400[spectrum]; [waves][spectrum] vstack[outAu]; [1:v][outAu] vstack [out]" -map "[out]" -map 0:0 -acodec pcm_s16le -vcodec libx264 -b 5m -vf format=yuv420p -f matroska "$3"
