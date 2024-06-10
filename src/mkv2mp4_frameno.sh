#!/bin/bash

src_dir="$1"
dst_dir="$2"

find "$src_dir" -name "*.mkv" -type f | while read -r file
do
    # echo $file
    rel_path=${file#"$src_dir"}
    dst_file="$dst_dir/${rel_path%.mkv}.mp4"
    mkdir -p "$(dirname "$dst_file")"
    echo ffmpeg -i "$file" -filter_complex "drawtext=text='%{n}': x=10: y=35: fontsize=48: fontcolor=red" "$dst_file"
    ffmpeg -i "$file" -filter_complex "drawtext=text='%{n}': x=10: y=35: fontsize=48: fontcolor=red" "$dst_file" -loglevel error
done