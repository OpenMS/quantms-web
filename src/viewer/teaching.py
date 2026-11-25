"""What each view and component of Explore results shows and how to read it.

The Viewer keeps its components free of explanations: each one's title carries
a short tooltip linking here, and the **Understand your Results** page
(``content/understand_results.py``) shows the full texts, one section per view
and per component. Component names are the catalog's (``src/viewer/catalog.py``).
"""

PAGE = "understand_results"

VIEWS = {
    "quantification": (
        "Which proteins differ in abundance between your groups. The summary counts the quantified "
        "proteins, the samples and groups, and how many proteins are higher or lower in the compared "
        "group. Pick a protein in the volcano plot or the group comparison table: its abundance in "
        "every sample appears beside the volcano plot, and it stays picked in every other view. "
        "PCA and the heatmap show whether the samples of one group resemble each other."
    ),
    "pathways": (
        "What the proteins that differ between your groups have in common. The volcano plot and the "
        "group comparison table show the differing proteins; GO enrichment tests which biological "
        "processes, molecular functions and cellular components are over-represented among them."
    ),
    "identification": (
        "Which proteins and peptides were identified, down to the spectrum behind each match. Read "
        "from left to right: picking a protein narrows the peptides to that protein, picking a peptide "
        "narrows the peptide-spectrum matches (PSMs), and picking a PSM shows its annotated spectrum "
        "below. Mass shifts show modifications the search did not look for."
    ),
    "quality_control": (
        "Whether the identifications and the group comparison can be trusted. The first rows check the "
        "files: how many spectra each identified, and which fall short. The next checks the search: "
        "target and decoy scores should separate, and precursor errors should centre near zero. The "
        "last rows check the steps of the group comparison: missing values, normalization, the values "
        "filled in for missing ones and the p-values."
    ),
}

COMPONENTS = {
    "dda_summary": (
        "Key figures of the identification: files, target peptide-spectrum matches (PSMs), distinct "
        "peptides and proteins after the false discovery rate (FDR) filters of the run, the share of MS2 "
        "spectra that gave a PSM, and the median precursor mass error. The grey line names the run's "
        "filters. A share of PSMs per MS2 spectrum far below 20 % usually points at a database or "
        "tolerance that does not fit the sample."
    ),
    "quant_summary": (
        "Key figures of the quantification: proteins with an abundance, samples, groups, the proteins "
        "significantly higher (up) or lower (down) in the compared group, and the share of abundances "
        "that were not measured. Up and down count proteins past the cut-offs of Edit comparison, by "
        "default |log2 fold change| ≥ 1 (at least twice or half as abundant) and an adjusted p-value ≤ 0.05."
    ),
    "sample_counts": (
        "PSMs, peptides and proteins per sample and per group. Samples of one group should reach "
        "similar counts; one far below the rest may have failed."
    ),
    "id_diagnostics": (
        "Every file's identification rate next to the rate expected for its search settings, the "
        "files that fall short first, each with its likely causes (wrong database or species, a "
        "precursor tolerance that cuts off true matches, spectra recorded in profile mode, few MS2 "
        "spectra)."
    ),
    "qc_scores": (
        "How many target and decoy PSMs reach each score, before the FDR filter. Decoys are matches to "
        "reversed or shuffled proteins and show what wrong matches score. A good search shows a second "
        "peak of target PSMs at high scores, apart from the decoys. The grey line names the score the "
        "PSM FDR lets pass."
    ),
    "qc_mass_error": (
        "The difference between measured and expected precursor mass of the target PSMs, in ppm, after "
        "removing picks of a 13C isotope peak. A peak away from zero means the instrument was "
        "miscalibrated; a flat spread up to the edges means the precursor tolerance is too wide or too "
        "narrow."
    ),
    "protein_table": (
        "Every identified protein with its peptides, PSMs, sequence coverage and the files it was found "
        "in. Contaminants (common lab proteins such as keratins and trypsin) are hidden; untick the "
        "filter to show them. Picking a protein picks it in every component."
    ),
    "protein_peptides": (
        "The sequence coverage of the picked protein (identified residues highlighted) and its "
        "peptides. Picking a peptide narrows the PSMs to it."
    ),
    "protein_psms": (
        "The peptide-spectrum matches of the picked protein and peptide: file, scan, charge, retention "
        "time, m/z and score. Picking one shows its annotated spectrum."
    ),
    "psm_spectrum": (
        "The MS2 spectrum of the picked PSM with its fragment ions annotated on the peptide sequence. "
        "Many matched b and y ions along the sequence make an identification convincing."
    ),
    "psms": "Every PSM after the FDR filter, with its spectrum, fragment ions and their mass errors.",
    "psms_rescored": "The PSMs as rescored by Percolator, which learns from targets and decoys to separate them.",
    "psms_search": "The raw hits of the database search, before rescoring and the FDR filter.",
    "mass_shifts": (
        "How far the measured precursor masses lie from the identified peptides. Shifts shared by many "
        "PSMs are modifications or artefacts the search did not look for; picking a shift lists its PSMs."
    ),
    "protein_abundance": "The abundance of every protein and peptide in every sample, as the quantification reported it.",
    "abundance_profile": (
        "The abundance of the picked protein in each sample (log2), coloured by group. It shows whether a "
        "difference between groups is consistent across the samples or comes from one outlier."
    ),
    "group_comparison": (
        "Every compared protein, most significant first: its log2 fold change (the difference of the "
        "groups' mean log2 abundance; 1 means twice as abundant), adjusted p-value and whether it counts "
        "as higher or lower. The adjusted p-value corrects for testing thousands of proteins at once. "
        "How the groups are compared (filter, missing values, normalization, test) is set with ⚖️ "
        "Edit comparison."
    ),
    "volcano": (
        "Each protein's log2 fold change (x) against its significance, -log10 p-value (y). Proteins far "
        "to the right or left and high up differ most reliably. The colours follow the cut-offs set in "
        "Edit comparison (what counts as changed). Clicking a point picks the protein everywhere; the "
        "buttons above the plot, or box and lasso select, pick a set of proteins for GO enrichment."
    ),
    "pca": (
        "Principal component analysis places each sample by its protein abundances, so that similar "
        "samples lie close together. Samples of one group should cluster; one far from its group may be "
        "an outlier."
    ),
    "heatmap": "The abundance of the most variable proteins in every sample, scaled per protein.",
    "clustered_heatmap": "The heatmap with proteins and samples ordered by similarity, so groups of proteins that behave alike stand out.",
    "go_enrichment": (
        "Gene Ontology (GO) terms over-represented among the proteins you pick, compared with all "
        "compared proteins: biological processes, molecular functions and cellular components. "
        "Pick proteins to test which GO terms are over-represented in them: use All changed, Higher "
        "or Lower above the volcano plot, box or lasso select in its mode bar (Shift adds), or the "
        "check boxes of the group comparison table (its check mark button turns them on). A term "
        "with a small adjusted p-value is shared by more picked proteins than chance would give. "
        "With nothing picked, it counts the terms of all compared proteins. The GO category and "
        "the number of terms in the chart are in the chart's settings (gear)."
    ),
    "comparison_missing_values": (
        "How many abundances are missing per sample and protein. Missing values are common in "
        "data-dependent acquisition and mostly belong to low-abundant proteins."
    ),
    "comparison_imputed_values": (
        "The values filled in for missing abundances, against the measured ones. Filled-in values should "
        "sit at the low end, where the missing values come from."
    ),
    "comparison_normalization": (
        "The abundance distribution of each sample before and after normalization. After it, the "
        "samples' medians should line up, so differences between samples come from biology rather than "
        "from how much was loaded."
    ),
    "comparison_pvalues": (
        "The distribution of the comparison's p-values. Flat with a peak near zero means some proteins "
        "differ; flat only means none do; a peak elsewhere hints at a test that does not fit the data."
    ),
    "file_chromatogram": (
        "The run's total ion chromatogram (TIC, the summed intensity of each MS1 spectrum) and base "
        "peak chromatogram (its highest peak) over retention time. Click a point to pick that spectrum. "
        "Gaps or a sudden drop point at a spray or chromatography problem."
    ),
    "file_scans": (
        "Every spectrum of the run: its scan number, MS level, retention time, precursor mass (MS2) "
        "and number of peaks. Pick a row to show the spectrum beside it."
    ),
    "file_spectrum": (
        "The picked spectrum as recorded, m/z against intensity, with its highest peaks labelled. "
        "Zoom in by dragging."
    ),
    "file_peak_map": (
        "The MS1 peaks of the run, retention time against mass, coloured by intensity. Peptides show "
        "as short streaks; zoom in to see their isotope patterns. The picked spectrum is marked."
    ),
    "file_psms": (
        "The peptide-spectrum matches (PSMs) of this run after the false discovery rate filter, with "
        "the annotated spectrum of the picked one."
    ),
}


def tooltip(name: str) -> str | None:
    """A component title's tooltip: the first sentence of its text and a link to the full one."""
    text = COMPONENTS.get(name)
    if not text:
        return None
    first = text.split(". ")[0].rstrip(".") + "."
    return f"{first}\n\n[More in Understand your Results](/{PAGE}#{anchor(name)})"


def help_text(name: str) -> str | None:
    """The full explanation of a component, for its info icon."""
    return COMPONENTS.get(name)


def anchor(name: str) -> str:
    return name.replace("_", "-")
