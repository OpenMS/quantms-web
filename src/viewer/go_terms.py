"""GO terms of a run's proteins, for the GO enrichment component (Streamlit-free).

The GO terms of the tested proteins are looked up once at MyGene.info and kept
as a term index (``openms_insight.analysis.enrichment.go_term_index``) in the
workspace, so testing another selection only counts and runs the tests.
"""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from openms_insight.analysis.enrichment import (
    fetch_go_annotations,
    go_enrichment_from_index,
    go_term_index,
)

CATEGORIES = {"BP": "Biological process", "CC": "Cellular component", "MF": "Molecular function"}


def clean(accession: str) -> str:
    """``sp|P12345|NAME_HUMAN`` -> ``P12345``, as the term index keys proteins."""
    parts = str(accession).split("|")
    return parts[1] if len(parts) >= 3 else str(accession)


def label(accession: str) -> str:
    """``sp|P12345|NAME_HUMAN`` -> ``NAME_HUMAN``, the name shown for a protein."""
    return str(accession).split("|")[-1]


def _index_path(folder: Path, proteins) -> Path:
    digest = hashlib.sha1("\n".join(sorted({clean(p) for p in proteins})).encode()).hexdigest()[:16]
    return Path(folder, f"go-index-{digest}.parquet")


ALL_PROTEINS = "go-index-all.parquet"


def _from_all(folder: Path, wanted: set) -> pl.DataFrame | None:
    """The rows of ``wanted`` from the index the workflow built, when it holds all of them."""
    path = Path(folder, ALL_PROTEINS)
    if not path.is_file():
        return None
    index = pl.read_parquet(path)
    if not wanted <= set(index["accession"].unique().to_list()):
        return None
    return index.filter(pl.col("accession").is_in(sorted(wanted)))


def prepare(folder: Path, proteins, fetch=fetch_go_annotations) -> str | None:
    """Look up the GO terms of ``proteins`` once, while the workflow runs.

    Adds them to the workspace's index of every protein processed so far
    (``ALL_PROTEINS``), so the GO enrichment component finds them without a
    network call. Returns an error message, or None.
    """
    wanted = {clean(p) for p in proteins if p}
    path = Path(folder, ALL_PROTEINS)
    known = pl.read_parquet(path) if path.is_file() else None
    missing = sorted(wanted - (set(known["accession"].unique().to_list()) if known is not None else set()))
    if not missing:
        return None
    try:
        annotations = fetch(missing)
    except Exception as err:  # the network, MyGene.info  # pylint: disable=broad-except
        return f"The GO terms could not be looked up at MyGene.info ({err})."
    index = go_term_index(annotations, missing)
    if known is not None:
        index = pl.concat([known, index.cast(known.schema)], how="vertical")
    path.parent.mkdir(parents=True, exist_ok=True)
    index.write_parquet(path)
    return None


def index_ready(folder: Path, proteins) -> bool:
    return _index_path(folder, proteins).is_file() or _from_all(folder, {clean(p) for p in proteins}) is not None


def term_index(folder: Path, proteins, fetch=fetch_go_annotations) -> tuple[pl.DataFrame | None, str | None]:
    """(term index of ``proteins``, error message or None), looked up once per protein set."""
    path = _index_path(folder, proteins)
    if path.is_file():
        return pl.read_parquet(path), None
    index = _from_all(folder, {clean(p) for p in proteins})
    if index is not None:
        return index, None
    try:
        annotations = fetch(sorted({clean(p) for p in proteins}))
    except Exception as err:  # the network, MyGene.info  # pylint: disable=broad-except
        return None, f"The GO terms could not be looked up at MyGene.info ({err})."
    index = go_term_index(annotations, sorted({clean(p) for p in proteins}))
    path.parent.mkdir(parents=True, exist_ok=True)
    index.write_parquet(path)
    return index, None


def enrichment(index: pl.DataFrame, selected, background) -> dict:
    """``go_enrichment_from_index`` of the selected proteins against the background."""
    return go_enrichment_from_index(index, [clean(p) for p in selected], [clean(p) for p in background])


def term_counts(index: pl.DataFrame, proteins) -> tuple[int, dict]:
    """(proteins with any GO term, category -> term, count, accessions), most frequent first."""
    wanted = sorted({clean(p) for p in proteins})
    terms = index.filter(pl.col("accession").is_in(wanted)).drop_nulls("term")
    annotated = terms["accession"].n_unique() if terms.height else 0
    tables = {
        category: terms.filter(pl.col("category") == category)
        .group_by("term")
        .agg(
            pl.col("accession").n_unique().alias("count"),
            pl.col("accession").unique().sort().str.join(";").alias("accessions"),
        )
        .sort(["count", "term"], descending=[True, False])
        for category in CATEGORIES
    }
    return annotated, tables


def _names(accessions: str, labels: dict) -> str:
    return ", ".join(labels.get(a, a) for a in str(accessions).split(";") if a)


def term_table(table: pd.DataFrame | None, labels: dict) -> pd.DataFrame:
    """The tested terms of one category, most significant first, proteins by name."""
    columns = ["GO term", "Selected", "Background", "Ratio", "Bg ratio", "Fold", "p-value", "adjusted p-value", "Proteins"]
    if table is None or not len(table):
        return pd.DataFrame(columns=columns)
    table = table.sort_values("p_adj")
    return pd.DataFrame(
        {
            "GO term": table["GO_Term"].astype(str),
            "Selected": table["Count"].astype(int),
            "Background": table["Background"].astype(int),
            "Ratio": table["GeneRatio"].astype(str),
            "Bg ratio": table["BgRatio"].astype(str),
            "Fold": table["FoldEnrichment"].astype(float).round(2),
            "p-value": table["p_value"].map(lambda v: f"{v:.3g}"),
            "adjusted p-value": table["p_adj"].map(lambda v: f"{v:.3g}"),
            "Proteins": [_names(p, labels) for p in table["Proteins"]],
        }
    ).reset_index(drop=True)


def term_bars(table: pd.DataFrame | None, count: int) -> pd.DataFrame:
    """The ``count`` most significant terms, for a bar chart of -log10(adjusted p)."""
    if table is None or not len(table):
        return pd.DataFrame({"term": [], "score": [], "hover": []})
    top = table.sort_values("p_adj").head(count)
    padj = top["p_adj"].astype(float).clip(lower=1e-300)
    return pd.DataFrame(
        {
            "term": top["GO_Term"].astype(str).tolist(),
            "score": (-np.log10(padj)).round(3).tolist(),
            "hover": [
                f"{t}: {c} selected ({r}), adj. p {p:.2e}, fold {f:.2f}"
                for t, c, r, p, f in zip(top["GO_Term"], top["Count"], top["GeneRatio"], padj, top["FoldEnrichment"])
            ],
        }
    )


def count_table(table: pl.DataFrame | None, labels: dict, annotated: int) -> pd.DataFrame:
    columns = ["GO term", "Proteins", "% of annotated", "Annotated proteins"]
    if table is None or not table.height:
        return pd.DataFrame(columns=columns)
    counts = table["count"].to_list()
    return pd.DataFrame(
        {
            "GO term": table["term"].to_list(),
            "Proteins": counts,
            "% of annotated": [round(100.0 * c / annotated, 1) if annotated else 0.0 for c in counts],
            "Annotated proteins": [_names(a, labels) for a in table["accessions"].to_list()],
        }
    )


def count_bars(table: pl.DataFrame | None, count: int) -> pd.DataFrame:
    if table is None or not table.height:
        return pd.DataFrame({"term": [], "proteins": []})
    top = table.head(count)
    return pd.DataFrame({"term": top["term"].to_list(), "proteins": top["count"].to_list()})
