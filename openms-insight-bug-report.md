# Bug Report: Heatmap-to-Table Row Selection Not Working

## Issue for openms-insight Repository

### Title
Table component does not highlight/scroll to row when Heatmap selection changes

### Description

When using a shared `StateManager` between `Heatmap` and `Table` components with matching `interactivity` identifiers, clicking a point in the Heatmap does not scroll to or highlight the corresponding row in the Table.

### Expected Behavior

1. User clicks a point in the Heatmap
2. Heatmap sets selection via `state_manager.set_selection("identification", id_idx_value)`
3. Table's Vue component (`TabulatorTable.vue`) detects the selection change
4. Table scrolls to the row with matching `id_idx` and highlights it
5. **All rows remain visible** (no filtering)

### Actual Behavior

Clicking a Heatmap point does NOT cause the Table to scroll or highlight the corresponding row. The selection appears to be set (SequenceView and LinePlot update correctly), but the Table does not respond.

### Reproduction Steps

```python
from openms_insight import Table, Heatmap, StateManager

# Initialize components with shared identifier
table = Table(
    cache_id="my_table",
    data=df.lazy(),
    cache_path=str(cache_dir),
    interactivity={"identification": "id_idx"},  # Shared identifier
    column_definitions=[...],
    index_field="id_idx",
)

heatmap = Heatmap(
    cache_id="my_heatmap",
    data=df.lazy(),
    cache_path=str(cache_dir),
    x_column="rt",
    y_column="mz",
    intensity_column="score",
    interactivity={"identification": "id_idx"},  # Same identifier
)

# Render with shared state
state_manager = StateManager()
heatmap(state_manager=state_manager, height=350)
table(state_manager=state_manager, height=533)
```

### Environment

- **openms-insight version**: >=0.1.10
- **Streamlit version**: 1.43.0
- **Browser**: Chrome/Firefox (tested both)
- **Python**: 3.12

### Root Cause Analysis (Deep Investigation)

After thorough code analysis, the root cause has been identified:

#### Primary Issue: Iterator Order Bug in `syncSelectionFromStore()`

**Location:** `TabulatorTable.vue` lines 1005-1034

The `syncSelectionFromStore()` function iterates through interactivity entries and uses the **first identifier with a non-null value**, then breaks:

```javascript
for (const [identifier, column] of Object.entries(interactivity)) {
    const selectedValue = this.selectionStore.$state[identifier]
    if (selectedValue !== undefined && selectedValue !== null) {
        // ... try to find row
        break  // STOPS after first non-null identifier
    }
}
```

**The problem:** When Table has multiple interactivity mappings (common case):
```python
interactivity={"file": "file_index", "spectrum": "scan_id", "identification": "id_idx"}
```

And Heatmap only sets:
```python
interactivity={"identification": "id_idx"}
```

**Failure scenario:**
1. User previously selected a file → "file" has value (e.g., 0)
2. User clicks Heatmap → sets "identification" = 50
3. `syncSelectionFromStore()` iterates
4. "file" has value 0 (non-null) → enters if block
5. Searches for row where `file_index === 0` (wrong identifier!)
6. **Breaks** - never checks "identification"

The function doesn't detect **which identifier changed**, it just uses the first one with a value.

#### Secondary Issue: Two-Render Cycle Timing

**Flow for external selection (Heatmap → Table):**

1. **First render:** Python cache miss → sends `dataChanged: false`
   - Vue selection store IS updated
   - `syncSelectionFromStore()` fires with stale `preparedTableData`
   - Row may not be found → stored as `pendingSelection`

2. **Second render:** Python cache hit → sends `dataChanged: true`
   - Selection store already has correct values (no change detected)
   - Watcher **doesn't fire** (values unchanged)
   - Navigation hints are in the data, but `navigateToPage` watcher may not trigger

This two-phase update cycle can cause the selection highlight to be missed.

### Suggested Fixes

#### Fix 1: Track Which Identifier Changed

Modify `syncSelectionFromStore()` to compare current vs previous selection state and prioritize the changed identifier:

```javascript
syncSelectionFromStore(): void {
    const interactivity = this.args.interactivity || {}

    // Find which identifier actually changed
    for (const [identifier, column] of Object.entries(interactivity)) {
        const selectedValue = this.selectionStore.$state[identifier]
        const previousValue = this.lastSyncedSelections?.[identifier]

        // Prioritize changed identifiers
        if (selectedValue !== previousValue && selectedValue != null) {
            // Handle this identifier first
            this.selectRowByColumn(column, selectedValue)
            break
        }
    }

    // Store current state for next comparison
    this.lastSyncedSelections = {...this.selectionStore.$state}
}
```

#### Fix 2: Ensure Selection Sync After Data Update

In the `currentDataHash` watcher, explicitly call `syncSelectionFromStore()` after data updates:

```javascript
currentDataHash: {
    handler(newHash: string, oldHash: string) {
        // ... existing logic ...

        // After data update, re-sync selection
        this.$nextTick(() => {
            this.syncSelectionFromStore()
        })
    }
}
```

### Verified Working Components

- **SequenceView**: Uses `filters` parameter → actively filters data based on selection (different mechanism)
- **LinePlot**: Linked via SequenceView → inherits working filter-based selection
- **Heatmap**: Correctly sets selection → Python StateManager updates correctly

### Related Code Paths

| File | Lines | Description |
|------|-------|-------------|
| `TabulatorTable.vue` | 989-1038 | `syncSelectionFromStore()` - main issue location |
| `TabulatorTable.vue` | 370-376 | Selection store watcher |
| `TabulatorTable.vue` | 378-410 | `navigateToPage` watcher |
| `table.py` | 711-856 | Python-side navigation hint generation |
| `bridge.py` | 456-797 | Render cycle with two-phase cache handling |
| `streamlit-data.ts` | 39-203 | Vue data store update logic |

### Workaround

Currently, using `filters` instead of `interactivity` works but **filters** data instead of highlighting:

```python
table = Table(
    ...
    filters={"identification": "id_idx"},  # Filters, not highlights
)
```
