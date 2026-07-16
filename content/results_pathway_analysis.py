import json
import mygene
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.io as pio
from collections import defaultdict        
from scipy.stats import fisher_exact
from pathlib import Path
from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data

# ================================
# Page setup
# ================================
params = page_setup()
st.title("ProteomicsLFQ Results")

# ================================
# Workspace check
# ================================
if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# ================================
# _run_go_enrichment function
# ================================
def _run_go_enrichment(pivot_df: pd.DataFrame, results_dir: Path):
    p_cutoff = 0.05
    fc_cutoff = 1.0

    analysis_df = pivot_df.dropna(subset=["p-value", "log2FC"]).copy()

    if analysis_df.empty:
        st.error("No valid statistical data found for GO enrichment.")
        st.write("❗ analysis_df is empty")
    else:
        with st.spinner("Fetching GO terms from MyGene.info API..."):
            mg = mygene.MyGeneInfo()

            def get_clean_uniprot(name):
                parts = str(name).split("|")
                return parts[1] if len(parts) >= 2 else parts[0]

            analysis_df["UniProt"] = analysis_df["protein"].apply(get_clean_uniprot)

            bg_ids = analysis_df["UniProt"].dropna().astype(str).unique().tolist()
            fg_ids = analysis_df[
                (analysis_df["p-value"] < p_cutoff) &
                (analysis_df["log2FC"].abs() >= fc_cutoff)
            ]["UniProt"].dropna().astype(str).unique().tolist()
            # st.write("✅ get_clean_uniprot applied")

            if len(fg_ids) < 3:
                st.warning(
                    f"Not enough significant proteins "
                    f"(p < {p_cutoff}, |log2FC| ≥ {fc_cutoff}). "
                    f"Found: {len(fg_ids)}"
                )
                st.write("❗ Not enough significant proteins")
            else:
                res_list = mg.querymany(
                    bg_ids, scopes="uniprot", fields="go", as_dataframe=False
                )
                res_go = pd.DataFrame(res_list)
                if "notfound" in res_go.columns:
                    res_go = res_go[res_go["notfound"] != True]

                def extract_go_terms(go_data, go_type):
                    if not isinstance(go_data, dict) or go_type not in go_data:
                        return []
                    terms = go_data[go_type]
                    if isinstance(terms, dict):
                        terms = [terms]
                    return list({t.get("term") for t in terms if "term" in t})

                for go_type in ["BP", "CC", "MF"]:
                    res_go[f"{go_type}_terms"] = res_go["go"].apply(
                        lambda x: extract_go_terms(x, go_type)
                    )

                annotated_ids = set(res_go["query"].astype(str))
                fg_set = annotated_ids.intersection(fg_ids)
                bg_set = annotated_ids
                # st.write(f"✅ fg_set bg_set are set")

                def run_go(go_type):
                    go2fg = defaultdict(set)
                    go2bg = defaultdict(set)

                    for _, row in res_go.iterrows():
                        uid = str(row["query"])
                        for term in row[f"{go_type}_terms"]:
                            go2bg[term].add(uid)
                            if uid in fg_set:
                                go2fg[term].add(uid)

                    records = []
                    N_fg = len(fg_set)
                    N_bg = len(bg_set)

                    for term, fg_genes in go2fg.items():
                        a = len(fg_genes)
                        if a == 0:
                            continue
                        b = N_fg - a
                        c = len(go2bg[term]) - a
                        d = N_bg - (a + b + c)

                        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
                        records.append({
                            "GO_Term": term,
                            "Count": a,
                            "GeneRatio": f"{a}/{N_fg}",
                            "p_value": p,
                        })

                    df = pd.DataFrame(records)
                    if df.empty:
                        return None, None

                    df["-log10(p)"] = -np.log10(df["p_value"].replace(0, 1e-10))
                    df = df.sort_values("p_value").head(20)

                    # ✅ Plotly Figure
                    fig = px.bar(
                        df,
                        x="-log10(p)",
                        y="GO_Term",
                        orientation="h",
                        title=f"GO Enrichment ({go_type})",
                    )

                    # st.write(f"✅ Plotly Figure generated")

                    fig.update_layout(
                        yaxis=dict(autorange="reversed"),
                        height=500,
                        margin=dict(l=10, r=10, t=40, b=10),
                    )

                    return fig, df

                go_results = {}

                for go_type in ["BP", "CC", "MF"]:
                    fig, df_go = run_go(go_type)
                    if fig is not None:
                        go_results[go_type] = {
                            "fig": fig,
                            "df": df_go
                        }
                # st.write(f"✅ go_type generated")

                go_dir = results_dir / "go-terms"
                go_dir.mkdir(parents=True, exist_ok=True)

                go_data = {}
                
                for go_type in ["BP", "CC", "MF"]:
                    if go_type in go_results:
                        fig = go_results[go_type]["fig"]
                        df = go_results[go_type]["df"]
                        
                        go_data[go_type] = {
                            "fig_json": fig.to_json(),  # Figure → JSON string
                            "df_dict": df.to_dict(orient="records")  # DataFrame → list of dicts
                        }
                
                go_json_file = go_dir / "go_results.json"
                with open(go_json_file, "w") as f:
                    json.dump(go_data, f)
                st.session_state["go_results"] = go_results
                st.session_state["go_ready"] = True if go_data else False
                # st.write("✅ GO enrichment analysis complete")

# ================================
# Load abundance data
# ================================
results_dir = Path(st.session_state["workspace"]) / "topp-workflow" / "results" / "quant_results"
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    st.info("Abundance data not available. Please run the workflow and configure sample groups first.")
    st.page_link("content/results_abundance.py", label="Go to Abundance", icon="📋")
    st.stop()

pivot_df, expr_df, group_map = result

go_json_file = results_dir / "go-terms" / "go_results.json"

go_input_df = pivot_df.copy()
if "ProteinName" in go_input_df.columns:
    go_input_df = go_input_df.rename(columns={"ProteinName": "protein"})

_run_go_enrichment(go_input_df, results_dir)

# ================================
# Tabs
# ================================
protein_tab, = st.tabs(["🧬 Protein Table"])
                
# ================================
# Protein-level results
# ================================
with protein_tab:
    st.markdown("### 🧬 Protein-Level Abundance Table")
    st.info(
        "This protein-level table is generated by grouping all PSMs that map to the "
        "same protein and aggregating their intensities across samples.\n\n"
        "Additionally, log2 fold change and p-values are calculated between sample groups."
    )

    if pivot_df.empty:
        st.info("No protein-level data available.")
    else:
        st.session_state["pivot_df"] = pivot_df
        st.dataframe(pivot_df.sort_values("p-value"), width="stretch")

# ======================================================
# GO Enrichment Results 
# ======================================================
st.markdown("---")
st.subheader("🧬 GO Enrichment Analysis")

if not go_json_file.exists():
    st.info("GO Enrichment results are not available yet. Please run the analysis first.")
else:    
    with open(go_json_file, "r") as f:
        go_data = json.load(f)
    
    bp_tab, cc_tab, mf_tab = st.tabs([
        "🧬 Biological Process",
        "🏠 Cellular Component",
        "⚙️ Molecular Function",
    ])
    
    for tab, go_type in zip([bp_tab, cc_tab, mf_tab], ["BP", "CC", "MF"]):
        with tab:
            if go_type not in go_data:
                st.info(f"No enriched {go_type} terms found.")
                continue
            
            fig_json = go_data[go_type]["fig_json"]
            df_dict = go_data[go_type]["df_dict"]
            
            fig = pio.from_json(fig_json)
            
            df_go = pd.DataFrame(df_dict)
            
            if df_go.empty:
                st.info(f"No enriched {go_type} terms found.")
            else:
                st.plotly_chart(fig, width="stretch")
                
                st.markdown(f"#### {go_type} Enrichment Results")
                st.dataframe(df_go, width="stretch")