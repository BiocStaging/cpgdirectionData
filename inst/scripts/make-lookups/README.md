# make-lookups: the three catalogue lookups

`lookup_blood_hg19`, `lookup_nasal_epithelium_hg19` and
`lookup_solid_tissue_hg19` pair every EPIC v2 CpG with candidate target genes
and give, for each pair, the direction predicted by a tissue-specific model
of measured eQTMs (a `PREDICTED` call or an `ABSTAIN`). Pairs measured in a
training cohort are `DIRECT_eQTM`, with the measured sign. The scripts here
rebuilt them in cpgdirectionData 0.99.5 (8 November 2026).

## Why they were rebuilt

The build that made the 0.99.0-0.99.4 tables was never committed. Its gene
table had two coordinate errors.

1. **Wrong chromosome.** Each gene had its correct TSS position but, for 46%
   of the blood lookup's Atlas genes (54% solid, 7% nasal), a wrong
   chromosome.
   - The Atlas's own `gene.chr` was not used. It is consistent for every gene
     and agrees with the HELIX transcript-cluster chromosome for 99.92% of
     genes.
   - Result: the shipped tables paired CpGs with a gene on another chromosome
     in 450,553 of 1,408,361 checkable blood rows, 288,290 of 829,224 solid
     rows and 52,327 of 1,395,873 nasal rows.
2. **Wrong genome build.** The TSS of 934 genes came from the Framingham and
   GTEx Atlas exports, which are not hg19.
   - Their CpG positions sit a median 342 kb and 542 kb from the hg19
     manifest.
   - The TSSs were a median 326 kb off.
   - Result: 8.0% of blood rows and 1.1% of solid rows were placed at those
     TSSs.

The training pools themselves are clean: no Atlas cis row and no TCGA pair
crosses chromosomes. The saved models therefore stand unchanged.

## What the rebuild keeps and what it changes

Everything is the original's except the gene placement.

**Kept as originally built:**
- the pools (`analysis/02` and `04` of eqtm-direction; 44,791 blood, 29,338
  nasal and 11,987 solid pairs, equal to the saved models' `n`);
- the saved `per_tissue_models.pkl` (models and isotonic calibrators,
  scikit-learn 1.7.2);
- the features (`lib_eqtm.full_features`);
- the candidate rule;
- the tiers and call thresholds;
- the `DIRECT_eQTM` rule.

All of these were recovered from the shipped tables. `build_lookups.py
--check` proves it: with the shipped gene placement, it regenerates each
shipped table identically in every column. The exceptions are 492 / 44 / 205
rows, where two genes share an identical TSS and the original tie order
cannot be recovered.

**Changed:** each gene's chromosome and hg19 TSS. The rules are in the header
of `build_lookups.py`.
- **Verification:** every gene used has its chromosome confirmed by the EPIC v2
  probe annotation and HGNC (aliases through HGNChelper). Where neither has
  the gene, it is confirmed by HELIX, the measured eQTMs and the gene's own cis
  training partners.
- **Agreement rates:** HELIX 99.92%, EPIC v2 99.86%, measured eQTMs 99.91%,
  training partners 99.80%.
- **Excluded:** genes whose sources conflict (5, all HELIX gene-label errors)
  and genes with no usable hg19 TSS (197).
- **TSS check:** the hg19 Atlas TSS agrees with the EPIC v2 Zhou GENCODE v41 TSS
  within 1 kb for 87-89% of genes, and within 10 kb for 94%.

## Result

| lookup | rows (was) | CpGs (was) | genes (was) | unchanged rows | rows on another chromosome |
|---|---|---|---|---|---|
| blood | 2,051,017 (1,663,799) | 672,298 (641,325) | 7,482 (6,572) | 937,029 | 0 |
| nasal epithelium | 1,523,103 (1,472,435) | 555,440 (553,145) | 4,712 (4,634) | 1,386,905 | 0 |
| solid tissue | 1,188,710 (849,229) | 542,032 (454,821) | 3,404 (2,973) | 510,098 | 0 |

The "was" figures are those of the 0.99.0-0.99.4 tables. "Unchanged" means
identical in every column. The last column is checked against the EPIC v2
annotation and HGNC.

Blood `DIRECT_eQTM` rows rise from 14,983 to 25,752. This is because the
measured pairs of genes that used to be misplaced can now be candidates.

md5 of the 0.99.5 resources:
`lookup_blood_hg19.rds 24361aebd5546db994c9971eeb8fc4e1`,
`lookup_nasal_epithelium_hg19.rds 18d6531edc83a23c0e4b1e97e7d0586f`,
`lookup_solid_tissue_hg19.rds 194e2d607bc74ee0d0be8f04c09c701b`.

## Reproduce

1. **Python environment.** Python 3 with scikit-learn 1.7.2 (the version the
   models were pickled with; the script refuses any other), pandas 2.3.3 and
   numpy 2.2.6:
   `python3 -m venv venv && venv/bin/pip install scikit-learn==1.7.2 pandas==2.3.3 numpy==2.2.6`
2. **R environment.** cpgdirection with its data layers, and HGNChelper
   (0.8.16 was used).
3. **Sources (RAW).**
   - The five eQTM Atlas exports, under their original file names:
     `EVA-PR_ All Subjects_data.csv`,
     `EVA-PR_ Atopic asthma + controls_data.csv`,
     `EVA-PR_ White Blood Cell_data.csv`,
     `Framingham Heart Study_data (1).csv` and `GTEx_ WholeBlood_data.csv`.
     They come from <https://shiny.crc.pitt.edu/eqtm_browser/> (Downloads
     tab).
   - The HELIX catalogue `eQTM_autosome_adj.cells_SIG.txt.gz` (Ruiz-Arenas et
     al. 2022).
   - The five `EWAS_OpenEWAS_EPICv1_block_*` folders (TCGA via OpenEWAS).
   - `ann_hg19.bed`, the EPIC v2 manifest lifted to hg19.
4. **Commands:**

```
Rscript export_inputs.R inputs <cpgdirection data dir>
python3 build_lookups.py --raw RAW --models per_tissue_models.pkl --inputs inputs --out out --check
Rscript write_rds.R out hub_upload <old hub_upload>
```

`lib_eqtm.py` is a verbatim copy of `eqtm-direction/analysis/lib_eqtm.py`
(md5 0defbd144a6b13995d1fdd296d7638b6). The run writes `out/build_report.csv`
(regression and comparison counts) and `out/gene_placement.csv.gz` (each
gene's placement, its sources, its independent checks and its verdict). A
full run takes about four minutes.
