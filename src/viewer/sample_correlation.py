"""Sample-to-sample correlation of protein abundances.

Streamlit-free. The correlation is computed on log2 intensities as measured,
before any missing value is filled in: a filled-in value is the same guess in
many samples and would make them look more alike than they are.
"""
import numpy as np
import pandas as pd

# Proteins two samples must both have measured for their correlation to count
MIN_SHARED_PROTEINS = 10


def sample_correlation(log2_mat: pd.DataFrame, min_shared: int = MIN_SHARED_PROTEINS) -> pd.DataFrame:
    """Pearson correlation between the columns (samples) of a protein x sample log2 matrix.

    Each pair uses the proteins measured in both samples; a pair sharing fewer
    than ``min_shared`` proteins gets NaN.
    """
    mat = log2_mat.replace([np.inf, -np.inf], np.nan)
    corr = mat.corr(method="pearson", min_periods=min_shared)
    corr.index.name = None
    corr.columns.name = None
    return corr


def shared_proteins(log2_mat: pd.DataFrame) -> pd.DataFrame:
    """Number of proteins measured in both samples, for every pair of samples."""
    present = log2_mat.replace([np.inf, -np.inf], np.nan).notna().astype(int)
    return present.T @ present


def median_shared_proteins(shared: pd.DataFrame) -> int:
    """Median number of proteins shared by two different samples."""
    values = shared.to_numpy(dtype=float, copy=True)
    np.fill_diagonal(values, np.nan)
    finite = values[~np.isnan(values)]
    return int(np.median(finite)) if finite.size else 0


def median_correlation(corr: pd.DataFrame) -> pd.Series:
    """Each sample's median correlation with the other samples, lowest first."""
    values = corr.to_numpy(dtype=float, copy=True)
    np.fill_diagonal(values, np.nan)
    medians = pd.Series(np.nanmedian(values, axis=1), index=corr.index) if len(corr) > 1 else pd.Series(dtype=float)
    return medians.sort_values()


def least_similar(corr: pd.DataFrame, gap: float = 0.05) -> list[str]:
    """Samples whose median correlation is ``gap`` or more below the median of all samples.

    A hint for samples that may not belong with the rest (a failed run, a
    swapped sample), not a test.
    """
    medians = median_correlation(corr).dropna()
    if len(medians) < 3:
        return []
    overall = float(medians.median())
    return [s for s, m in medians.items() if overall - m >= gap]
