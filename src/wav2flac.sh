#!/bin/bash

# Navigate to the directory containing the .wav files
cd /path/to/your/directory

# Create the log files
touch t_converted_log.txt
touch t_error_log.txt

# Find all .wav files in the directory and its subdirectories
find . -type f -name "*.wav" -print0 | while IFS= read -r -d '' file; do
    # Use ffmpeg to convert the .wav file to .flac
    echo ffmpeg -i ""$file"" ""${file%.*}.flac"" -hide_banner '&&' rm "$file" '&&' echo "${file%.*}.flac" '>>' t_converted_log.txt '||' echo "$file" '>>' t_error_log.txt >> s.sh
done