#!/usr/bin/env bash
if [ $# -eq 0 ]
then
  echo 'combineVidAu_wav.sh [videoFile_path] [acousticFile_path] [output_newVideoFile_path]'
  exit 0
fi


# ffmpeg -i "$1" -i "$2" -map 0:v -map 1:a -vcodec copy -acodec copy -f matroska "$3.mkv"
ffmpeg -i "$1" -i "$2" -map 0:v -map 1:a -vf format=yuv420p -ar 96000 "$3" 
# ffmpeg -i "$1" -i "$2" -map 0:v -map 1:a -vf format=yuv420p -ar 96000 -af loudnorm=I=-16:LRA=11:TP=-1.5 "$3" 
#  -c:v libx264 -crf 25 -ar 96000 -af loudnorm=I=-16:LRA=11:TP=-1.5 