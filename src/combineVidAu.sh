#!/usr/bin/env bash
if [ $# -eq 0 ]
then
  echo 'combineVidAu.sh [videoFile_path] [acousticFile_path] [output_newVideoFile_path]'
  exit 0
fi


ffmpeg -i "$1" -f s16le -ar 400000 -ac 5 -i "$2" -map 0:v -map 1:a -vcodec copy -acodec copy -f matroska "$3.mkv"