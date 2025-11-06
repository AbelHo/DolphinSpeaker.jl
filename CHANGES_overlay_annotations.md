# Summary of Changes to overlay_annotations_on_video

## Overview
Enhanced the `overlay_annotations_on_video` function in `src/video.jl` to support:
1. **Array-based properties** for multiple annotation sets
2. **Per-row properties** from DataFrame columns using `:in_annotations`
3. **Full backward compatibility** with existing code

## Changes Made

### 1. Function Signature Updated
- Changed parameter types from strictly typed (e.g., `radius::Int`) to flexible unions
- Now accepts scalars, arrays, or `:in_annotations` symbol
- Parameters affected: `radius`, `default_color`, `default_alpha`, `default_shape`

### 2. New Preprocessing Logic
- Detects whether `:in_annotations` is used for each property
- Converts scalar values to arrays for uniform handling
- Validates array lengths match number of annotation sets

### 3. Enhanced normalize_item Function
- Now takes `item_idx` parameter to access correct array element
- Checks for property columns in DataFrame (r, g, b, a, radius, shape)
- Returns tuple indicating whether properties are per-row or uniform
- Tuple format: `(frames, pts, per_row::Bool, color_data, alpha_data, radius_data, shape_data)`

### 4. Updated Drawing Loops
- Modified both VideoIO mode and stream/frames mode loops
- Handles both per-row and uniform property cases
- Extracts correct properties based on `per_row` flag

## Key Features

### Feature 1: Array-Based Properties
```julia
overlay_annotations_on_video(
    [df1, df2],
    video_path,
    output_path;
    default_color = ["red@0.5", "blue@0.7"],
    radius = [30, 50],
    default_shape = [:circle, :rect]
)
```

### Feature 2: DataFrame Column Properties
```julia
# DataFrame with r, g, b, a, radius, shape columns
overlay_annotations_on_video(
    df,
    video_path,
    output_path;
    default_color = :in_annotations,
    radius = :in_annotations,
    default_alpha = :in_annotations,
    default_shape = :in_annotations
)
```

### Feature 3: Backward Compatibility
```julia
# Old code still works
overlay_annotations_on_video(
    "annotations.csv",
    video_path,
    output_path;
    radius = 25,
    default_color = "red@0.5"
)
```

## Testing

Created two new files for testing and documentation:
1. `test_overlay_annotations.jl` - Test script with multiple examples
2. `docs/overlay_annotations_features.md` - Comprehensive documentation

## DataFrame Column Requirements

When using `:in_annotations`, the following columns are expected:

| Argument | Required Columns | Format |
|----------|------------------|--------|
| `default_color` | `r`, `g`, `b` (optional: `a`) | Float 0.0-1.0 |
| `default_alpha` | `a` | Float 0.0-1.0 |
| `radius` | `radius` | Integer pixels |
| `default_shape` | `shape` | Symbol (`:circle`, `:rect`) |

## Benefits

1. **Flexibility**: Support multiple annotation styles in one video
2. **Per-row control**: Each detection/annotation can have unique properties
3. **Confidence visualization**: Use alpha/size to show detection confidence
4. **Multi-object tracking**: Different colors/shapes for different objects
5. **Backward compatible**: No breaking changes to existing code

## Code Quality

- All changes maintain existing code structure
- Proper error handling with informative messages
- Array length validation
- Graceful fallbacks (e.g., missing columns use defaults)
- Clear documentation in docstring
