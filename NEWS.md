# cpgdirectionData 0.99.5

* The three catalogue lookups (`lookup_blood_hg19`,
  `lookup_nasal_epithelium_hg19`, `lookup_solid_tissue_hg19`) are rebuilt.
  The build behind 0.99.0-0.99.4 had two errors in its gene table:
  - a wrong chromosome for 46% of the blood lookup's Atlas genes (54% solid,
    7% nasal), so a CpG was paired with a gene on another chromosome in
    450,553 of 1,408,361 checkable blood rows, 288,290 of 829,224 solid rows
    and 52,327 of 1,395,873 nasal rows;
  - for 934 genes, a TSS taken from Atlas exports that are not hg19 (a median
    326 kb off).
* The pools, the saved models, the features, the candidate rule, the tiers
  and the `DIRECT_eQTM` rule are unchanged. With the old gene placement the
  new pipeline regenerates the old tables identically, except 492 / 44 / 205
  rows where two genes share an identical TSS.
* Every gene used now has its chromosome confirmed by the EPIC v2 annotation
  and HGNC, or, where neither has the gene, by HELIX, the measured eQTMs and
  its own training partners. Rows on another chromosome: 0.
* New sizes: blood 2,051,017 pairs (672,298 CpGs; was 1,663,799 / 641,325),
  nasal 1,523,103 (555,440; was 1,472,435 / 553,145), solid 1,188,710
  (542,032; was 849,229 / 454,821). Blood `DIRECT_eQTM` rows: 25,752 (was
  14,983).
* The pipeline that builds them is `inst/scripts/make-lookups/` (README,
  `build_lookups.py`, `export_inputs.R`, `write_rds.R`); `make-data.R` no
  longer converts the lookups.
* The other seven resources are byte-identical to 0.99.4.

# cpgdirectionData 0.99.4

* Zenodo deposit moved to version 2.99.8 (record 23257261; concept DOI
  10.5281/zenodo.22225163 unchanged). `distal_links` was rebuilt and now
  carries scE2G single-cell enhancer-gene predictions (Sheth et al. 2026)
  for nine PBMC cell types at E2G.Score.qnorm >= 0.177 beside the
  ENCODE-rE2G sets: 10,205,122 rows, 1,289,742 CpG-gene pairs, 525,222 CpGs.
  The other nine resources are byte-identical to v2.99.6. Every RDataPath
  in metadata.csv now points at the new record; SourceVersion is 2.99.8.

# cpgdirectionData 0.99.3

* ExperimentHub metadata in the form requested by the Hubs team: `Location_Prefix`
  is `https://zenodo.org/` and `RDataPath` the complete
  `records/23122949/files/<name>.rds` path, without `?download=1`.
* `cpgdData()` resolves a resource by its exact Hub title within
  `query(eh, "cpgdirectionData")` and, when the Hub is reachable, says how many
  cpgdirectionData resources it lists.
* The Hub chunk of the vignette is evaluated.

# cpgdirectionData 0.99.2

* Two new layers, both opt-in on the cpgdirection side and neither entering
  `best_direction`:
  - `distal_links` (hg38): ENCODE-rE2G thresholded enhancer-gene predictions
    for 26 resting blood cell types and brain regions/cell types (47
    released ENCODE annotation sets), intersected with EPIC v2 probe
    positions. 8,036,254 rows, 1,149,157 CpG-gene pairs, 514,549 CpGs.
    Built by `tools/build_distal_links.R` in cpgdirection; the ABC-model
    tables ENCODE ships beside the rE2G tables are excluded so one score
    scale applies.
  - `onco_eqtm_consensus` (hg19, 450K): Onco-eQTM cross-tumour cis-eQTM
    consensus (Korra et al. 2026; 27 TCGA cancer types). 2,464,990 CpG-gene
    pairs with consensus direction, agreement tier (O1 113,560; O2 566,351;
    O3 549,150; O4 1,235,929), number of cancer types and correlation
    summaries recovered from the t statistics. Sanity check: promoter pairs
    81-86% negative, gene-body pairs 39% negative.
* `epicv2_probe_gene_annotation` rebuilt on the Exeter re-annotated manifest
  v3.0 (GENCODE 49): 1,006,212 pairs (793,233 CpGs; previously 996,697), 9,647
  pairs gained, 132 lost, 14,214 TSS200/TSS1500 relabels; the Exeter
  Promoter_2000bp/Enhancer_5000bp labels now live in `exeter_regulatory`.
* `getDistalLinks()` and `getOncoEqtmConsensus()` accessors.

# cpgdirectionData 0.99.0

* First submission. Serves the eight cpgdirection data layers that exceed
  Bioconductor's per-file size cap: the three tissue catalogue lookups
  (blood, nasal epithelium, solid tissue; hg19), the SMR causal-direction
  table, the EPIC v2 hg19 probe positions, the brain SMR directions, the
  saliva bridge scores, and the multi-source EPIC v2 probe-to-gene
  annotation (five cross-validated tracks).
* `cpgdData(name)` plus one named accessor per layer;
  `options(cpgdirectionData.local=)` serves the same resources from a local
  directory for air-gapped machines and pinned analyses.
