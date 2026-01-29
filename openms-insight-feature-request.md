# Feature Request: Heatmap Rendering Order Control

## Summary

Add parameters to control the rendering order (z-order) of points in the Heatmap component, allowing users to specify which points should be drawn on top when overlapping.

## Use Case

When visualizing peptide-spectrum match (PSM) identification results, we display a heatmap with:
- **x-axis**: Retention time (RT)
- **y-axis**: Mass-to-charge ratio (m/z)
- **intensity**: Score (e-value, PEP, or q-value)

For these score types, **lower values indicate better identifications**. Currently, when points overlap, there's no control over which points are rendered on top. Ideally, we want the best identifications (lowest scores) to be visible on top of worse ones.

## Current Workaround

We can reverse the colorscale using `colorscale="Portland_r"` to make low scores appear with high-intensity colors, but this doesn't address the overlapping point visibility issue.

## Proposed Solution

Add one or more of the following parameters to the `Heatmap` component:

### Option 1: `sort_column` and `sort_ascending`

```python
Heatmap(
    data=df,
    x_column="rt",
    y_column="mz",
    intensity_column="score",
    sort_column="score",      # Column to sort by for rendering order
    sort_ascending=False,     # False = high values rendered first (low on top)
)
```

### Option 2: `render_order`

```python
Heatmap(
    data=df,
    x_column="rt",
    y_column="mz",
    intensity_column="score",
    render_order="intensity_asc",  # Options: "intensity_asc", "intensity_desc", "data_order"
)
```

### Option 3: `top_layer`

```python
Heatmap(
    data=df,
    x_column="rt",
    y_column="mz",
    intensity_column="score",
    top_layer="low_intensity",  # Options: "low_intensity", "high_intensity"
)
```

## Additional Context

- This is particularly important for proteomics workflows where score distributions often have many overlapping points
- The feature would complement the existing `colorscale` parameter for full control over visual representation
- Consider how this interacts with the multi-resolution downsampling - perhaps the sort should be applied before binning/downsampling

## Environment

- openms-insight version: >= 0.1.10
- Application: quantms-web (DDA-LFQ) built on streamlit-template
