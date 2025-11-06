# Visual Guide to overlay_annotations_on_video Modes

## Three Ways to Specify Properties

```
┌─────────────────────────────────────────────────────────────────────┐
│                    INPUT MODES OVERVIEW                              │
└─────────────────────────────────────────────────────────────────────┘

Mode 1: SCALAR (Uniform)
─────────────────────────
All annotations use same properties

    radius = 25
    default_color = "red@0.5"
    default_shape = :circle
    
    Annotation 1 → [25, red, circle]
    Annotation 2 → [25, red, circle]
    Annotation 3 → [25, red, circle]


Mode 2: ARRAY (Per-Set)
────────────────────────
Each annotation set gets unique properties

    radius = [20, 35, 50]
    default_color = ["red@0.5", "blue@0.7", "green@0.6"]
    default_shape = [:circle, :rect, :circle]
    
    Annotation 1 → [20, red, circle]
    Annotation 2 → [35, blue, rect]
    Annotation 3 → [50, green, circle]


Mode 3: :in_annotations (Per-Row)
──────────────────────────────────
Each row in DataFrame has unique properties

    default_color = :in_annotations
    radius = :in_annotations
    default_shape = :in_annotations
    
    DataFrame with columns: frame, px, py, r, g, b, radius, shape
    
    Row 1 → [radius=20, color=(1.0,0.0,0.0), shape=:circle]
    Row 2 → [radius=35, color=(0.0,1.0,0.0), shape=:rect]
    Row 3 → [radius=50, color=(0.0,0.0,1.0), shape=:circle]
```

## DataFrame Structure Examples

### Basic DataFrame (for Mode 1 & 2)
```
┌────────┬──────┬──────┐
│ frame  │  px  │  py  │
├────────┼──────┼──────┤
│   0    │ 100  │ 100  │
│   1    │ 110  │ 105  │
│   2    │ 120  │ 110  │
└────────┴──────┴──────┘
```

### DataFrame with Color Columns (for Mode 3)
```
┌────────┬──────┬──────┬─────┬─────┬─────┬─────┐
│ frame  │  px  │  py  │  r  │  g  │  b  │  a  │
├────────┼──────┼──────┼─────┼─────┼─────┼─────┤
│   0    │ 100  │ 100  │ 1.0 │ 0.0 │ 0.0 │ 0.5 │  ← Red
│   1    │ 110  │ 105  │ 0.0 │ 1.0 │ 0.0 │ 0.7 │  ← Green
│   2    │ 120  │ 110  │ 0.0 │ 0.0 │ 1.0 │ 0.9 │  ← Blue
└────────┴──────┴──────┴─────┴─────┴─────┴─────┘
```

### Full DataFrame with All Properties (for Mode 3)
```
┌────────┬──────┬──────┬─────┬─────┬─────┬─────┬────────┬─────────┐
│ frame  │  px  │  py  │  r  │  g  │  b  │  a  │ radius │  shape  │
├────────┼──────┼──────┼─────┼─────┼─────┼─────┼────────┼─────────┤
│   0    │ 100  │ 100  │ 1.0 │ 0.0 │ 0.0 │ 0.5 │   20   │ :circle │
│   1    │ 110  │ 105  │ 0.0 │ 1.0 │ 0.0 │ 0.7 │   30   │ :rect   │
│   2    │ 120  │ 110  │ 0.0 │ 0.0 │ 1.0 │ 0.9 │   40   │ :circle │
└────────┴──────┴──────┴─────┴─────┴─────┴─────┴────────┴─────────┘
```

## Use Cases

### Use Case 1: Multi-Object Tracking
```
Track multiple objects with distinct colors

┌─────────────┐         ┌─────────────┐
│   Object A  │         │   Object B  │
│  (Red 🔴)   │         │  (Blue 🔵)  │
└─────────────┘         └─────────────┘
      ↓                       ↓
   DataFrame 1            DataFrame 2

annotations = [df_objectA, df_objectB]
default_color = ["red@0.6", "blue@0.6"]
radius = [25, 35]
```

### Use Case 2: Confidence Visualization
```
Show detection confidence using alpha/size

Frame 0: ◯ (low confidence, small, transparent)
Frame 1: ⬤ (medium confidence, medium, semi-opaque)
Frame 2: ⬤ (high confidence, large, opaque)

default_alpha = :in_annotations  ← Read 'a' column
radius = :in_annotations         ← Read 'radius' column
```

### Use Case 3: State Changes
```
Visualize object state changes using color/shape

State 1: Red Circle    → Idle
State 2: Yellow Rect   → Moving
State 3: Green Circle  → Active

default_color = :in_annotations  ← Read 'r,g,b' columns
default_shape = :in_annotations  ← Read 'shape' column
```

## Decision Tree

```
┌─────────────────────────────────────────────────────┐
│  Do you have multiple annotation sets?              │
└────────────────┬────────────────────────────────────┘
                 │
        ┌────────┴────────┐
        │ YES             │ NO
        ↓                 ↓
┌───────────────┐   ┌──────────────────────────────┐
│ Use ARRAY     │   │ Does each ROW need unique    │
│ mode          │   │ properties?                  │
│               │   └──────┬───────────────────────┘
│ radius = [...] │         │
│ color = [...]  │    ┌────┴────┐
└───────────────┘    │ YES     │ NO
                     ↓         ↓
              ┌──────────┐  ┌─────────────┐
              │ Use      │  │ Use SCALAR  │
              │ :in_     │  │ mode        │
              │ annota-  │  │             │
              │ tions    │  │ radius = 25 │
              │ mode     │  │ color = ... │
              └──────────┘  └─────────────┘
```

## Mixing Modes (Hybrid)

You can mix different modes for different properties:

```julia
overlay_annotations_on_video(
    df,
    video_path,
    output_path;
    default_color = :in_annotations,    # ← Per-row from columns
    radius = :in_annotations,           # ← Per-row from columns
    default_alpha = 0.7,                # ← Uniform scalar
    default_shape = :circle             # ← Uniform scalar
)
```

This gives maximum flexibility for your visualization needs!
