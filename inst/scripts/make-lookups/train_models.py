#!/usr/bin/env python3
"""Step 2 of the v2 catalogue build: fit, calibrate and save the three per-tissue models.

The procedure of eqtm-direction/analysis/04_per_tissue_models.py: lib_eqtm features
(|distance| terms, neighbour prior, gene prior), lib_eqtm.gbm() (HistGradientBoosting, seed 0),
5-fold GroupKFold by CpG for out-of-fold predictions, an isotonic calibrator fitted on them, the
model refitted on the whole pool, and leave-one-cohort-out AUCs. One change: early stopping is
off, so every fit runs lib_eqtm.gbm()'s 400 iterations on all its rows. With scikit-learn's
default (early stopping on above 10,000 rows) the solid pool's fold models (about 9,500 rows)
ran 400 iterations while its final model (11,935) stopped early on 90% of the rows, so the
calibrator was fitted to a different model than the one shipped; every final model also left
10% of its pool out. What differs is the input: the
harmonized pools of build_pools.py, whose gene prior is keyed by gene identity (gene_id) and
whose distances come from the TSS the lookups use. Deterministic: the same pools give the same
models and calibrators.

Usage:
  python train_models.py --pools OUT_OF_BUILD_POOLS --out MODELS.pkl [--lib DIR] [--report CSV]
Requires scikit-learn 1.7.2.
"""
import argparse, os, pickle, sys
import numpy as np, pandas as pd
ap = argparse.ArgumentParser()
ap.add_argument("--pools", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--lib", default=None); ap.add_argument("--report", default=None)
a = ap.parse_args()
if a.lib: sys.path.insert(0, a.lib)
import lib_eqtm as L, sklearn
def gbm():
    return L.gbm().set_params(early_stopping=False)
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
assert sklearn.__version__ == "1.7.2", f"scikit-learn 1.7.2 required, this is {sklearn.__version__}"
out, rep = {}, []
for t in ["blood", "nasal_epithelium", "solid_tissue"]:
    P = pd.read_csv(os.path.join(a.pools, f"pool_{t}.csv.gz"))
    P = P.assign(gene=P.gene_id)[["cpg_id", "gene", "dist", "y", "chr", "pos", "cohort", "tissue"]].reset_index(drop=True)
    oof = np.zeros(len(P))
    for tr, te in GroupKFold(5).split(P, P.y, P.cpg_id.values):
        oof[te] = gbm().fit(L.full_features(P.iloc[tr], P.iloc[tr]).values, P.y.values[tr]) \
                         .predict_proba(L.full_features(P.iloc[te], P.iloc[tr]).values)[:, 1]
    iso = IsotonicRegression(out_of_bounds="clip").fit(oof, P.y.values)
    X = L.full_features(P, P); m = gbm().fit(X.values, P.y.values)
    out[t] = dict(model=m, calibrator=iso, features=list(X.columns), n=len(P), base_rate=float(P.y.mean()),
                  auc_oof=float(roc_auc_score(P.y.values, oof)))
    rep.append(dict(tissue=t, held_out="(grouped 5-fold)", n=len(P), AUC=out[t]["auc_oof"]))
    print(f"{t:<18} n={len(P):>6}  grouped-CV AUC {out[t]['auc_oof']:.4f}  base rate {out[t]['base_rate']:.3f}", flush=True)
    for co in sorted(P.cohort.unique()):
        te_, tr_ = P[P.cohort == co], P[P.cohort != co]
        p = gbm().fit(L.full_features(tr_, tr_).values, tr_.y.values).predict_proba(L.full_features(te_, tr_).values)[:, 1]
        rep.append(dict(tissue=t, held_out=co, n=len(te_), AUC=roc_auc_score(te_.y.values, p)))
pickle.dump(out, open(a.out, "wb"))
if a.report: pd.DataFrame(rep).to_csv(a.report, index=False)
print("wrote", a.out)
