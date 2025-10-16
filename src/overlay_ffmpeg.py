import pandas as pd
import os

csv_path = "detection_pixels.csv"
df = pd.read_csv(csv_path)
filters = []
radius = 100 #20  # diameter of the circle

for _, row in df.iterrows():
    frame = int(round(row['frame']))
    x = int(round(row['px']))
    y = int(round(row['py']))
    filters.append(
        f"drawbox=x={x-radius//2}:y={y-radius//2}:w={radius}:h={radius}:color=red@0.7:t=fill:enable='eq(n,{frame})'"
        # f"drawbox=x={x-radius//2}:y={y-radius//2}:w={radius}:h={radius}:color=red@0.7:t=2:enable='eq(n,{frame})'"
    )
filters = ",".join(filters)
print(filters)

os.system(f'ffmpeg -i /media/spin/anas2/data_res/dolphin/calf/temp/delete/1/combined__1.GoPro_Clicker.MP4.mp4 -vf "{filters}" -codec:a copy output_with_circles_fill.mkv')

# ffmpeg -i /media/spin/anas2/data_res/dolphin/calf/temp/delete/1/combined__1.GoPro_Clicker.MP4.mp4 -vf "$FILTERS" -codec:a copy output_with_circles.mp4