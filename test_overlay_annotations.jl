# Test script for overlay_annotations_on_video with new array and :in_annotations features

using DataFrames, CSV

# Example 1: Using arrays for multiple annotation sets with different properties
println("Example 1: Multiple annotation sets with array-based properties")
println("=" ^ 60)

# Create two annotation DataFrames
df1 = DataFrame(frame=[0, 1, 2, 3, 4], px=[100, 110, 120, 130, 140], py=[100, 100, 100, 100, 100])
df2 = DataFrame(frame=[0, 1, 2, 3, 4], px=[200, 210, 220, 230, 240], py=[200, 200, 200, 200, 200])

annotations = [df1, df2]

# Use arrays to specify different properties for each annotation set
# df1 will have red circles with radius 30
# df2 will have blue rectangles with radius 50
println("""
overlay_annotations_on_video(
    annotations,
    "input_video.mp4",
    "output_array_example.mp4";
    default_color = ["red@0.5", "blue@0.5"],
    radius = [30, 50],
    default_shape = [:circle, :rect],
    default_alpha = [0.5, 0.7]
)
""")

# Example 2: Using :in_annotations to read properties from DataFrame columns
println("\nExample 2: Reading properties from DataFrame columns")
println("=" ^ 60)

# Create a DataFrame with per-row color, alpha, radius, and shape
df_with_props = DataFrame(
    frame = [0, 1, 2, 3, 4],
    px = [150, 160, 170, 180, 190],
    py = [150, 150, 150, 150, 150],
    r = [1.0, 0.0, 0.0, 1.0, 0.0],      # red, green, blue, yellow, cyan
    g = [0.0, 1.0, 0.0, 1.0, 1.0],
    b = [0.0, 0.0, 1.0, 0.0, 1.0],
    a = [0.5, 0.6, 0.7, 0.8, 0.9],
    radius = [20, 30, 40, 50, 60],
    shape = [:circle, :circle, :rect, :rect, :circle]
)

println("""
overlay_annotations_on_video(
    df_with_props,
    "input_video.mp4",
    "output_in_annotations.mp4";
    default_color = :in_annotations,
    default_alpha = :in_annotations,
    radius = :in_annotations,
    default_shape = :in_annotations
)
""")

# Example 3: Hybrid approach - some properties from columns, some uniform
println("\nExample 3: Hybrid - some properties from DataFrame, some uniform")
println("=" ^ 60)

df_hybrid = DataFrame(
    frame = [0, 1, 2, 3, 4],
    px = [250, 260, 270, 280, 290],
    py = [250, 250, 250, 250, 250],
    r = [1.0, 0.5, 0.0, 0.5, 1.0],      # gradient from red to blue
    g = [0.0, 0.0, 0.0, 0.0, 0.0],
    b = [0.0, 0.5, 1.0, 0.5, 0.0],
    radius = [25, 35, 45, 35, 25]       # growing then shrinking
)

println("""
overlay_annotations_on_video(
    df_hybrid,
    "input_video.mp4",
    "output_hybrid.mp4";
    default_color = :in_annotations,      # read r,g,b from DataFrame
    radius = :in_annotations,             # read radius from DataFrame
    default_alpha = 0.6,                  # uniform alpha for all
    default_shape = :circle               # uniform shape for all
)
""")

# Example 4: Backward compatibility - old usage still works
println("\nExample 4: Backward compatibility - simple uniform properties")
println("=" ^ 60)

df_simple = DataFrame(frame=[0, 1, 2, 3], px=[100, 150, 200, 250], py=[100, 150, 200, 250])

println("""
overlay_annotations_on_video(
    "annotations.csv",
    "input_video.mp4",
    "output_simple.mp4";
    default_color = "yellow@0.5",
    radius = 40,
    default_shape = :circle
)
""")

println("\n" * "=" ^ 60)
println("All examples defined! To run, uncomment and provide actual video files.")
println("=" ^ 60)

# Save example DataFrames to CSV for reference
CSV.write("example_annotations_simple.csv", df1)
CSV.write("example_annotations_with_props.csv", df_with_props)
CSV.write("example_annotations_hybrid.csv", df_hybrid)
println("\nExample CSV files saved:")
println("  - example_annotations_simple.csv")
println("  - example_annotations_with_props.csv")
println("  - example_annotations_hybrid.csv")
