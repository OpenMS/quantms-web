# PXD044981 mini example: Yeast + UPS2 (1x vs 10x), DDA, 3-min RT slice

Source: PRIDE PXD044981 (Q Exactive HF, HCD, DDA; doi:10.1021/acs.jproteome.4c00048).
Raw files converted with ThermoRawFileParser (centroided, indexed mzML), then cut with
OpenMS FileFilter `-rt` to one 3-min window per run. The window is the same peptide
region in each run (reference B000272, 99–102 min), shifted per run by its RT offset:

| File | Condition (UPS2) | RT kept (s) |
|---|---|---|
| B000248_..._ratio01_DDA   | 1x  | 6300–6480 |
| B000272_..._ratio01_DDA_2 | 1x  | 5940–6120 |
| B000296_..._ratio01_DDA_3 | 1x  | 6066–6246 |
| B000264_..._ratio10_DDA   | 10x | 5988–6168 |
| B000288_..._ratio10_DDA_2 | 10x | 6006–6186 |
| B000312_..._ratio10_DDA_3 | 10x | 6060–6240 |

Note: `ratio10` = 10x UPS2 (about 6x more UPS2 PSMs in the full runs). The community
SDRF on PRIDE labels these the other way round.

Database: `db/yeast_ups2.fasta` = UniProt reviewed S. cerevisiae (UP000002311, 6067) + 48 UPS2
proteins (`_HUMAN`). Targets only; let the workflow add decoys.

Expected result with quantms-web LFQ defaults (Comet -> Percolator -> IDFilter q<=0.1 ->
ProteomicsLFQ 1% PSM/protein FDR): ~900–980 PSMs/run, 594 protein groups,
13 UPS2 proteins; UPS2 log2(10x/1x) ~ +3.4 to +3.9 (expected 3.32) vs yeast median ~ +0.4.
