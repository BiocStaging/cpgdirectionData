# cpgdirectionData 0.99.8

* **Two new resources; the ten existing ones are untouched.**
  - `smr_directions_v2` (hg19): the second SMR layer, read by cpgdirection
    1.1.0 and later. GoDMC mQTL instruments against three expression
    outcomes (eQTLGen cis-eQTL, FinnGen immune-multiome PBMC and cell types,
    INTERVAL blood RNA-seq): 768,862 CpG-gene pairs over 133,498 CpGs
    (413,982 with an eQTLGen record, 135,334 seen in FinnGen alone, 219,546
    in INTERVAL alone). Each pair carries its per-outcome records, a
    replication rung (`smr_rung`: replicated 108,135; high 14,284; moderate
    528,204; weak 52,878; reported 65,361) and a per-pair accuracy
    (`smr_accuracy`) from a leave-one-chromosome-out model against measured
    blood eQTMs (25,875 validated pairs; AUC 0.850; mean calibration gap
    0.009). Built by `tools/build_smr_layer_interval.py` in cpgdirection.
    Accessor: `getSmrDirectionsV2()`.
  - `cpg_coordinates` (hg19 and hg38): array membership and coordinates for
    every probe ID of the Illumina HM27, 450K, EPIC v1 and EPIC v2 arrays
    (1,084,297 IDs), read by `cpgd_cpg_annotation()` in cpgdirection. hg38
    from the Zhou-lab InfiniumAnnotation (SeSAMe) manifests, hg19 from the
    Illumina 450K manifest or the EPIC v1 coordinates; a position is stored
    only where `epicv2_probe_gene_annotation` (hg38) or `cpg_positions_hg19`
    does not already give it (352,454 hg38, 141,945 hg19); 1,094 probes
    that two sources place differently in hg38 are flagged. Built by
    `tools/build_cpg_coordinates.R` in cpgdirection. Accessor:
    `getCpgCoordinates()`.
* **Why the second SMR layer has its own name.** cpgdirection 0.99.x and
  1.0.x take every row of `smr_directions` at its tier. Served under that
  name, the second layer would hand those versions 354,880 pairs they were
  never validated on, 50,798 of them pairs the second layer itself marks
  `reported` (not to be used for a direction). `smr_directions` (413,982
  pairs) therefore stays exactly as deposited, and `getSmrDirections()`
  still returns it.
* `metadata.csv`: the ten existing rows are unchanged and keep Zenodo
  record 24091940 (deposit 2.99.10); the two new rows point to the record of
  deposit 2.99.11. `inst/scripts/make-metadata.R` takes resource names after
  the record and the deposit version and then rewrites only those rows.
* `inst/scripts/make-data.R` converts the two new resources and no longer
  rebuilds `smr_directions`.

# cpgdirectionData 0.99.7

* **The catalogue's training chain is rebuilt and committed.** 0.99.5 kept
  the original training pools and the saved models, whose build was never
  committed. Checked against Ensembl GRCh37, those pools had:
  - wrong-end distances: 2,654 Framingham minus-strand blood rows (5.9%),
    measured from a point a median 15.8 kb from the Ensembl TSS;
  - non-hg19 TSSs in 156 TCGA solid rows;
  - one gene under several names (GTEx ENSG IDs, Framingham clone names,
    EVA-PR symbols, renamed HGNC symbols), each treated as its own gene by
    the pools, the gene prior and the lookups;
  - a pseudo-gene, "NAN": 3,935 HELIX rows from 1,461 unnamed transcript
    clusters.

  The saved calibrators and the solid model could not be regenerated. The
  original procedure also fitted each calibrator to models unlike the final
  one: with scikit-learn's early stopping (on above 10,000 rows), the solid
  fold models ran 400 iterations, the final solid model stopped early, and
  every final model left 10% of its pool out.
* **`inst/scripts/make-lookups/` now holds the whole chain**, and it is
  deterministic: `export_inputs.R`, `build_pools.py`, `train_models.py`,
  `build_lookups.py` and `write_rds.R`.
  - **Identity.** Each identified name resolves to one Ensembl GRCh37 gene
    (annotables `grch37`, HGNC previous and alias symbols through
    HGNChelper). The match is checked against the cohort's own coordinate
    wherever the cohort gives one. 11,005 of 11,478 pool names give 10,725
    genes.
  - **TSS.** One 1-based, strand-aware TSS per gene: the cohort's own hg19 TSS
    (EVA-PR, else HELIX) where it agrees with the identity, else the Ensembl
    gene TSS. Training distances and lookups use the same TSS.
  - **Models.** The three models are retrained on the harmonized pools by the
    original procedure, with early stopping off: 400 iterations on all rows,
    so each calibrator matches its model. Pools: blood 43,982 rows, nasal
    28,978, solid 11,935. Out-of-fold AUCs are 0.848 / 0.865 / 0.662. The
    blood rise from 0.793 comes from removing the "NAN" pseudo-gene; on
    identified genes the AUC is 0.848 before and 0.851 after.
  - **Lookups.** The candidate rule, tiers, call thresholds and
    `DIRECT_eQTM` rule are unchanged; a gene's names are one candidate.
  - **Solid tier B abstains** (PI ruling). The solid model's tier-B calls
    agree 0.41 with the solid yardsticks (123 pairs, 95% CI 0.32-0.49), below
    chance.
* **New sizes.**

  | lookup | pairs (0.99.5) | CpGs (0.99.5) |
  |---|---|---|
  | blood | 1,970,053 (2,051,017) | 668,711 (672,298) |
  | nasal | 1,500,942 (1,523,103) | 552,167 (555,440) |
  | solid | 1,184,630 (1,188,710) | 541,536 (542,032) |

* **Validation.** Against independent yardsticks (MESA, HELIX Dryad, EVA-PR
  nasal, FUSION, GTEx solid), with genes matched by Ensembl ID:
  - blood and nasal hold, with differences the size refitting alone
    produces;
  - where a tissue model calls a yardstick pair in both builds, 7 of 14,670
    signs differ;
  - solid is lower: GTEx solid −0.021 at ladder level. The solid model's
    tier-A calls agree 0.644 with FUSION and 0.618 with GTEx solid (was
    0.677 and 0.674). The solid model is unstable: over 11 retrainings it
    ranged 0.60-0.69 with FUSION.

  This is a correction, not a gain in accuracy.
* The other seven resources are byte-identical to 0.99.5. The build outputs
  (pools, models, reports) are deposited with Zenodo deposit 2.99.10.
* 0.99.6 was prepared and withdrawn before release.

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
