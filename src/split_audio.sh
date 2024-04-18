input_file="$1"  # Replace with your input file
output_file="$2"  # Replace with your output file prefix
duration=360  # Duration for each split (in seconds)

# Get the length of the input file in seconds
length=$(ffprobe -v error -show_entries format=duration -of default=noprint_wrappers=1:nokey=1 "$input_file")
length=${length%.*}  # Remove the decimal part

# Split the file
for start in $(seq 0 $duration $length); do
    ffmpeg -ss $start -t $duration -i "$input_file" -acodec copy "${output_file}_$start.wav"
done