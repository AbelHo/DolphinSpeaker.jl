# Enhanced overlay_annotations_on_video Function

## New Features

The `overlay_annotations_on_video` function has been enhanced with flexible property specification for overlays. You can now:

1. **Use arrays** to specify different properties for each annotation set
2. **Use `:in_annotations`** to read properties from DataFrame/CSV columns (per-row customization)
3. **Mix and match** - use arrays for some properties and `:in_annotations` for others
4. **Backward compatibility** - all existing code continues to work

## Property Arguments

The following arguments now support multiple modes:

- `default_color`: Color specification
- `default_alpha`: Alpha/opacity value (0.0 to 1.0)
- `radius`: Size of the overlay marker
- `default_shape`: Shape of the overlay (`:circle`, `:rect`, etc.)

## Usage Modes

### Mode 1: Single Value (Backward Compatible)

All annotations use the same properties:

```julia
overlay_annotations_on_video(
    annotations,
    "input.mp4",
    "output.mp4";
    default_color = "red@0.5",
    radius = 25,
    default_alpha = 0.5,
    default_shape = :circle
)
```

### Mode 2: Array of Values

Each annotation set gets its own properties:

```julia
# Two annotation sets with different properties
df1 = DataFrame(frame=[0,1,2], px=[100,110,120], py=[100,100,100])
df2 = DataFrame(frame=[0,1,2], px=[200,210,220], py=[200,200,200])

overlay_annotations_on_video(
    [df1, df2],
    "input.mp4",
    "output.mp4";
    default_color = ["red@0.5", "blue@0.7"],      # df1=red, df2=blue
    radius = [30, 50],                             # df1=30px, df2=50px
    default_alpha = [0.5, 0.7],                    # df1=0.5, df2=0.7
    default_shape = [:circle, :rect]               # df1=circles, df2=rectangles
)
```

**Note**: Array length must match the number of annotation sets.

### Mode 3: Read from DataFrame Columns (`:in_annotations`)

Each row in the DataFrame can have its own properties:

```julia
# DataFrame with per-row properties
df = DataFrame(
    frame = [0, 1, 2, 3, 4],
    px = [100, 110, 120, 130, 140],
    py = [100, 100, 100, 100, 100],
    r = [1.0, 0.0, 0.0, 1.0, 0.5],      # Red, Green, Blue, Yellow, Orange
    g = [0.0, 1.0, 0.0, 1.0, 0.5],
    b = [0.0, 0.0, 1.0, 0.0, 0.0],
    a = [0.5, 0.6, 0.7, 0.8, 0.9],
    radius = [20, 25, 30, 35, 40],
    shape = [:circle, :circle, :rect, :rect, :circle]
)

overlay_annotations_on_video(
    df,
    "input.mp4",
    "output.mp4";
    default_color = :in_annotations,    # Read r,g,b columns
    default_alpha = :in_annotations,    # Read a column
    radius = :in_annotations,           # Read radius column
    default_shape = :in_annotations     # Read shape column
)
```

**Required columns when using `:in_annotations`:**
- `default_color = :in_annotations`: Requires `r`, `g`, `b` columns (values 0.0-1.0). Optional `a` column for alpha.
- `default_alpha = :in_annotations`: Requires `a` column (values 0.0-1.0)
- `radius = :in_annotations`: Requires `radius` column (integer pixel values)
- `default_shape = :in_annotations`: Requires `shape` column (symbols like `:circle`, `:rect`)

### Mode 4: Hybrid Approach

Mix uniform properties with per-row customization:

```julia
df = DataFrame(
    frame = [0, 1, 2, 3],
    px = [100, 110, 120, 130],
    py = [100, 100, 100, 100],
    r = [1.0, 0.5, 0.0, 0.5],     # Gradient from red to blue
    g = [0.0, 0.0, 0.0, 0.0],
    b = [0.0, 0.5, 1.0, 0.5],
    radius = [20, 30, 40, 30]      # Growing then shrinking
)

overlay_annotations_on_video(
    df,
    "input.mp4",
    "output.mp4";
    default_color = :in_annotations,    # Per-row colors from r,g,b
    radius = :in_annotations,           # Per-row radius
    default_alpha = 0.6,                # Uniform alpha
    default_shape = :circle             # Uniform shape
)
```

## Color Specifications

Colors can be specified in multiple formats:

1. **Named colors with alpha**: `"red@0.5"`, `"blue@0.8"`, `"yellow@0.3"`
2. **Hex colors with alpha**: `"#ff0000@0.5"`, `"#00f@0.8"`
3. **RGB tuples**: `(1.0, 0.0, 0.0)` for red, `(0.0, 1.0, 0.0)` for green
4. **RGBA tuples**: `(1.0, 0.0, 0.0, 0.5)` for red with 50% opacity
5. **CSV format**: `"1.0,0.0,0.0"` for red
6. **From DataFrame**: Use `r`, `g`, `b` columns with values 0.0-1.0

Available named colors:
- red, green, blue, yellow, white, black, cyan, magenta, orange

## Examples

### Example 1: Multiple Tracking Targets

Track two objects with different colors and sizes:

```julia
using DataFrames, CSV

# Object 1: small red circles
obj1 = DataFrame(frame=[0,1,2,3], px=[100,105,110,115], py=[100,100,100,100])

# Object 2: large blue rectangles  
obj2 = DataFrame(frame=[0,1,2,3], px=[300,305,310,315], py=[200,200,200,200])

overlay_annotations_on_video(
    [obj1, obj2],
    "tracking.mp4",
    "tracking_annotated.mp4";
    default_color = ["red@0.6", "blue@0.6"],
    radius = [15, 40],
    default_shape = [:circle, :rect]
)
```

### Example 2: Confidence-Based Visualization

Show detection confidence using alpha values:

```julia
df = DataFrame(
    frame = [0, 1, 2, 3, 4],
    px = [100, 110, 120, 130, 140],
    py = [100, 105, 110, 115, 120],
    r = fill(1.0, 5),
    g = fill(0.0, 5),
    b = fill(0.0, 5),
    a = [0.3, 0.5, 0.8, 0.9, 1.0]      # Increasing confidence
)

overlay_annotations_on_video(
    df,
    "detections.mp4",
    "detections_confidence.mp4";
    default_color = :in_annotations,
    default_alpha = :in_annotations
)
```

### Example 3: Growing/Shrinking Markers

Animate marker size over time:

```julia
df = DataFrame(
    frame = 0:20,
    px = fill(200, 21),
    py = fill(200, 21),
    r = fill(0.0, 21),
    g = fill(1.0, 21),
    b = fill(0.0, 21),
    radius = [20 + 2*i for i in 0:20]  # Growing from 20 to 60 pixels
)

overlay_annotations_on_video(
    df,
    "video.mp4",
    "video_growing_marker.mp4";
    default_color = :in_annotations,
    radius = :in_annotations
)
```

## Backward Compatibility

All existing code continues to work without modification:

```julia
# This still works exactly as before
overlay_annotations_on_video(
    "annotations.csv",
    "input.mp4", 
    "output.mp4";
    radius = 25,
    default_color = "red@0.5"
)
```

## Testing

Run the test script to see examples:

```julia
include("test_overlay_annotations.jl")
```

This will create example CSV files demonstrating the different usage modes.
