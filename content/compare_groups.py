"""Compare Groups Page: which proteins change between the experimental groups."""
import json

import streamlit as st

from src.common.common import page_setup, save_params
from src.common.differential import STEPS, comparison_groups, count_significant
from src.common.guide import (
    ensure_group_comparison,
    needs_results,
    next_steps,
    page_intro,
    page_link,
    show_comparison_settings,
)
from src.common.insight_tables import show_insight_table
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map

params = page_setup()
st.title("Compare groups")

page_intro(
    "Which proteins are more or less abundant in one group of samples than in the other? "
    "The comparison runs automatically with settings that suit most experiments, label-free and TMT alike.",
    background="""
The comparison runs four steps on the protein abundance table:

1. **Filter proteins.** Proteins seen in only a few samples cannot be compared
   reliably, so a protein is kept only if it was measured in at least half of the
   samples of at least one group.
2. **Fill in missing values.** A protein that is missing in a sample was usually
   too faint to detect there. It gets that protein's lowest measured value, so a
   protein switched off in one group still shows up as down-regulated.
3. **Normalize.** Intensities are log2-transformed (so a difference of 1 means
   twice as much) and each sample is shifted so its median matches the others,
   which corrects for unequal amounts of material loaded.
4. **Test.** A moderated t-test (similar to *limma*) compares the groups protein
   by protein. Because thousands of proteins are tested, the p-values are
   adjusted with the Benjamini-Hochberg method: an **adjusted p-value of 0.05**
   means about 5 % of the proteins called significant are expected to be false
   positives.

**log2 fold change** is the difference of the group means on the log2 scale:
+1 means twice as abundant, −1 half as abundant.
""",
    background_title="What happens here, in plain words",
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

result = get_abundance_data(st.session_state["workspace"])
if result is None:
    needs_results()

pivot_df, _, group_map = result
id_col = get_id_column(st.session_state["workspace"], pivot_df)
sample_group_map = get_sample_group_map(st.session_state["workspace"], pivot_df, group_map)
sample_cols = [c for c in pivot_df.columns if c in sample_group_map]

st.subheader("Your groups")
groups: dict[str, list[str]] = {}
for s in sample_cols:
    groups.setdefault(sample_group_map[s], []).append(s)
cols = st.columns(max(len(groups), 1))
for col, (group, members) in zip(cols, sorted(groups.items())):
    col.metric(group, f"{len(members)} sample{'' if len(members) == 1 else 's'}")
if "Unassigned" in groups:
    st.warning(
        "Some samples have no group, so they are left out of the comparison. Give each a group on "
        "the **Describe experiment** page (for TMT, type **skip** for channels you do not want "
        "compared). The comparison updates without running the analysis again."
    )
    page_link("configure")

if st.button(
    "Start over with the recommended settings",
    help="Discard steps changed on the Fine-tune pages and compare the groups again with the defaults.",
):
    for key in [*STEPS, "normalized_transform", "comparison_settings"]:
        st.session_state.pop(key, None)
    with open("default-parameters.json", "r", encoding="utf-8") as f:
        defaults = json.load(f)
    for key, value in defaults.items():
        if key.startswith("postproc-"):
            params[key] = value
            st.session_state.pop(key, None)
    save_params(params)

statistics_df = ensure_group_comparison(params)
if statistics_df is None:
    st.stop()

direction = comparison_groups(sample_group_map, sample_cols)
up, down = count_significant(statistics_df)
st.subheader("Result")
m = st.columns(3)
m[0].metric("Proteins compared", f"{statistics_df['p-adj'].notna().sum():,}")
if direction:
    ref, comp = direction
    m[1].metric(f"Higher in {comp}", f"{up:,}", help="Adjusted p-value ≤ 0.05 and at least twice as abundant.")
    m[2].metric(f"Lower in {comp}", f"{down:,}", help="Adjusted p-value ≤ 0.05 and at most half as abundant.")
    note = f"A positive log2 fold change means higher in **{comp}** than in **{ref}**."
    if len(groups) > 2:
        note += (
            " With more than two groups the p-value tests whether any group differs; the fold "
            f"change still compares **{comp}** with **{ref}**."
        )
    st.caption(note + " Counts use |log2 fold change| ≥ 1 and adjusted p-value ≤ 0.05; explore other cut-offs on the volcano plot.")
show_comparison_settings()

st.subheader("All compared proteins, most significant first")
table = statistics_df.sort_values("p-adj", na_position="last")[[id_col, "log2FC", "p-adj", "p-value"]].copy()
if direction:
    sig = (table["p-adj"] <= 0.05) & (table["log2FC"].abs() >= 1)
    table.insert(1, "Change", "–")
    table.loc[sig & (table["log2FC"] > 0), "Change"] = f"higher in {direction[1]}"
    table.loc[sig & (table["log2FC"] < 0), "Change"] = f"lower in {direction[1]}"
for col in ("p-adj", "p-value"):
    # As text with three significant digits: the table stores floats as float32,
    # which would show 0.155 as 0.1550000011920929. Rows are already in p-value order.
    table[col] = table[col].map(lambda v: "" if v != v else f"{v:.3g}")
table = table.rename(columns={
    id_col: "Protein", "log2FC": "log2 fold change", "p-adj": "adjusted p-value",
})
show_insight_table(
    table, name="compare-groups", layout="fitColumns",
    column_overrides={
        "log2 fold change": {"formatter": "money", "formatterParams": {"thousand": False, "precision": 2}},
    },
)
st.caption("The table's download button saves it as CSV, in the current sort and filter.")

next_steps(
    ("volcano", "See every protein's fold change against its significance and pick cut-offs."),
    ("pca", "Check that samples of the same group look alike."),
    ("enrichment", "Find biological processes over-represented among the changed proteins."),
)
