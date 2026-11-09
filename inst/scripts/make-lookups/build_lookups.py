#!/usr/bin/env python3
"""Step 3 of the v2 catalogue build: the three lookups (lookup_<tissue>_hg19).

Inputs are the outputs of steps 1 and 2: the harmonized pools and gene table of
build_pools.py and the models of train_models.py. Nothing is placed or
identified here; every gene is the Ensembl GRCh37 gene of the gene table, on
its chromosome at the table's TSS, the same TSS its training distances used.

Candidates: the Ensembl genes of the tissue's training pool (HELIX transcript
clusters with no identified gene train the model but are never candidates).
On the CpG's chromosome, genes sorted by TSS (ties: the old map order, then
gene ID); the 3 genes before and the 4 after the CpG's insertion point, kept
when |pos - TSS| <= 250 kb. tss_dist = pos - TSS (genomic, not strand-signed,
as in every shipped lookup). target_gene = the gene table's display name.

Calls: tier A = a catalogued CpG within 50 kb and the gene catalogued; B = one
of the two (every candidate is a pool gene, so in practice B = no catalogued
CpG within 50 kb, and tier C does not occur); PREDICTED when tier A and
confidence >= 0.30 or tier B and confidence >= 0.50, else ABSTAIN; a (CpG,
gene) pair of the training pool is DIRECT_eQTM with the pool's majority sign
(ties +1), pairs matched by gene ID. The gene prior counts the CpG's own pool
pairs, as the shipped build did. One tissue-specific rule (PI ruling,
2026-11-09): the solid-tissue model's tier-B calls ABSTAIN. Against the solid
yardsticks they agree 0.41 (123 pairs, 95% CI 0.32-0.49; FUSION 0.38), below
chance: far from any catalogued CpG the weak solid model extrapolates.

Usage:
  python build_lookups.py --pools OUT_OF_BUILD_POOLS --models MODELS.pkl --raw RAW --out OUT
      [--lib DIR] [--tissues blood,nasal_epithelium,solid_tissue]
RAW holds ann_hg19.bed. Writes OUT/lookup_<tissue>.csv.gz (the shipped columns),
OUT/lookup_<tissue>_gene_id.csv.gz (target_gene -> gene_id, for checks) and
OUT/build_report.csv. write_rds.R turns the lookups into the Hub resources.
Requires scikit-learn 1.7.2.
"""
import argparse, os, pickle, sys, time
import numpy as np, pandas as pd
ap = argparse.ArgumentParser()
ap.add_argument("--pools", required=True); ap.add_argument("--models", required=True)
ap.add_argument("--raw", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--lib", default=None); ap.add_argument("--tissues", default="blood,nasal_epithelium,solid_tissue")
a = ap.parse_args()
if a.lib: sys.path.insert(0, a.lib)
import lib_eqtm as L, sklearn
assert sklearn.__version__ == "1.7.2", f"scikit-learn 1.7.2 required, this is {sklearn.__version__}"
os.makedirs(a.out, exist_ok=True)
T0 = time.time()
def log(*x): print(f"[{time.time() - T0:6.0f}s]", *x, flush=True)
W, NL, NR = 250000, 3, 4
TIER_B_CALL = {"blood": 0.50, "nasal_epithelium": 0.50, "solid_tissue": np.inf}   # solid tier B never calls

ann = pd.read_csv(os.path.join(a.raw, "ann_hg19.bed"), sep="\t", header=None, comment="t", names=["chr", "s", "e", "name", "sc", "st"])
ann["cpg_id"] = ann.name.astype(str).str.split("_").str[0].str.lower(); ann["pos"] = ann.s + 2
MAN = ann.drop_duplicates("cpg_id")[["cpg_id", "chr", "pos"]].reset_index(drop=True)
GT = pd.read_csv(os.path.join(a.pools, "gene_table.csv.gz"))
assert GT.gene_id.is_unique and GT.display.is_unique and GT.tss.notna().all() and GT.chr.notna().all()
MODELS = pickle.load(open(a.models, "rb"))

def build(t):
    P = pd.read_csv(os.path.join(a.pools, f"pool_{t}.csv.gz"))
    P = P.assign(gene=P.gene_id)[["cpg_id", "gene", "dist", "y", "chr", "pos", "cohort", "tissue"]].reset_index(drop=True)
    MOD = MODELS[t]; assert MOD["n"] == len(P), f"{t}: the model was trained on {MOD['n']} rows, the pool has {len(P)}"
    G = GT[GT.gene_id.isin(set(P.gene))]
    parts = []
    for ch, m in MAN.groupby("chr"):
        g = G[G.chr == ch].sort_values(["tss", "map_order", "gene_id"], kind="stable")
        if not len(g): continue
        Tt = g.tss.values.astype(float); gid = g.gene_id.values; nm = g.display.values
        Pp = m.pos.values.astype(float); ids = m.cpg_id.values
        ins = np.searchsorted(Tt, Pp, side="left")
        for off in range(-NL, NR):
            j = ins + off; ok = (j >= 0) & (j < len(Tt)); jj = j[ok]; d = Pp[ok] - Tt[jj]; keep = np.abs(d) <= W
            parts.append(pd.DataFrame({"cpg_id": ids[ok][keep], "gene_id": gid[jj][keep], "target_gene": nm[jj][keep],
                                       "tss_dist": d[keep], "chr": ch, "pos": Pp[ok][keep]}))
    Q = pd.concat(parts, ignore_index=True)
    q = pd.DataFrame({"cpg_id": Q.cpg_id.values, "gene": Q.gene_id.values, "dist": Q.tss_dist.values,
                      "chr": Q.chr.values, "pos": Q.pos.values})
    qg = q.copy(); qg["cpg_id"] = qg.cpg_id + "|q"
    X = pd.concat([L.full_features(q, P)[["abs", "log", "p500", "p2k", "p10k", "p100k"]],
                   L.neighbour_prior(q, P), L.gene_prior(qg, P)], axis=1)
    p = MOD["calibrator"].predict(MOD["model"].predict_proba(X[MOD["features"]].values)[:, 1])
    near = X.nn50000.values > 0; gene = X.gn.values > 0
    tier = np.where(near & gene, "A", np.where(near | gene, "B", "C"))
    conf = np.abs(p - 0.5) * 2
    call = ((tier == "A") & (conf >= 0.30)) | ((tier == "B") & (conf >= TIER_B_CALL[t]))
    ym = P.groupby(["cpg_id", "gene"]).y.mean().reindex(pd.MultiIndex.from_arrays([q.cpg_id, q.gene])).values
    direct = ~np.isnan(ym)
    status = np.where(direct, "DIRECT_eQTM", np.where(call, "PREDICTED", "ABSTAIN"))
    direction = np.where(direct, np.where(ym >= 0.5, 1.0, -1.0), np.where(call, np.where(p >= 0.5, 1.0, -1.0), np.nan))
    R = pd.DataFrame({"cpg_id": Q.cpg_id.values, "target_gene": Q.target_gene.values, "tss_dist": Q.tss_dist.astype(int).values,
                      "status": status, "direction": direction, "probability_plus1": np.round(p, 4),
                      "confidence": np.round(conf, 4), "evidence_tier": tier, "gene_id": Q.gene_id.values})
    R = R.sort_values(["cpg_id", "tss_dist", "gene_id"], kind="stable").reset_index(drop=True)
    assert not R.duplicated(["cpg_id", "gene_id"]).any()
    pool_pairs = P[~P.gene.str.startswith("HELIX:")].drop_duplicates(["cpg_id", "gene"])
    in_lookup = pd.MultiIndex.from_arrays([pool_pairs.cpg_id, pool_pairs.gene]).isin(pd.MultiIndex.from_arrays([R.cpg_id, R.gene_id]))
    rep = dict(tissue=t, rows=len(R), cpgs=int(R.cpg_id.nunique()), genes=int(R.gene_id.nunique()), candidate_genes=len(G),
               pool_gene_pairs=len(pool_pairs), pool_pairs_in_lookup=int(in_lookup.sum()),
               **{f"n_{k}": int(v) for k, v in R.status.value_counts().items()},
               **{f"tier_{k}": int(v) for k, v in R.evidence_tier.value_counts().items()})
    return R, rep

report = []
for t in a.tissues.split(","):
    R, rep = build(t)
    R.drop(columns="gene_id").to_csv(os.path.join(a.out, f"lookup_{t}.csv.gz"), index=False)
    R[["target_gene", "gene_id"]].drop_duplicates().to_csv(os.path.join(a.out, f"lookup_{t}_gene_id.csv.gz"), index=False)
    report.append(rep); log(t, rep)
pd.DataFrame(report).to_csv(os.path.join(a.out, "build_report.csv"), index=False)
log("done")
