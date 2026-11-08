#!/usr/bin/env python3
"""Rebuild the three cpgdirection catalogue lookups (lookup_<tissue>_hg19).

The original build script was never committed. This script reconstructs it
from what survives (eqtm-direction/analysis/lib_eqtm.py and
04_per_tissue_models.py, the saved per_tissue_models.pkl, the Atlas, HELIX
and OpenEWAS inputs) and proves the reconstruction by regenerating the
shipped tables (--check). It then rebuilds them with each gene placed on its
own chromosome at an hg19 TSS (the default), changing nothing else: same
training pools, same saved models and calibrators, same features, same
candidate rule, same tiers and call thresholds.

Candidate rule (recovered from the shipped tables): on the CpG's chromosome,
genes sorted by TSS; the 3 genes before and the 4 after the CpG's insertion
point, kept when |pos - TSS| <= 250 kb. tss_dist = pos - TSS. Only genes of
the tissue's training pool are candidates.

Calls: tier A = a catalogued CpG within 50 kb and the gene catalogued; B =
one of the two; PREDICTED when tier A and confidence >= 0.30 or tier B and
confidence >= 0.50, else ABSTAIN; a (CpG, gene) pair of the training pool is
DIRECT_eQTM with the pool's majority sign (ties +1). The gene-prior feature
counts the CpG's own pool pairs (the shipped DIRECT rows are reproduced only
that way; for every other row it makes no difference).

Gene placement (default build):
  chromosome  Atlas gene.chr; else the HELIX transcript cluster's chromosome;
              else the chromosome on which the EPIC v2 annotation and the
              gene's cis training partners agree.
  TSS         the Atlas first-row TSS when that row's file is hg19 (the three
              EVA-PR exports); the Framingham and GTEx exports are not hg19,
              so genes first seen there take the EPIC v2 Zhou GENCODE v41 TSS
              carried to hg19 through the manifest, else the HELIX
              TC_gene_TSS, else the Atlas TSS carried to hg19 through the
              gene's own cis rows (rows agreeing within 1 kb); HELIX-only
              genes keep the HELIX TSS.
  verified    the chromosome agrees with every independent annotation that
              has the gene (EPIC v2 probe annotation, HGNC with aliases), or,
              where neither has it, with every eQTM source that does (HELIX,
              the measured eQTMs, the gene's cis training partners). Only
              verified genes with an in-bounds hg19 TSS are candidates.

Usage:
  python build_lookups.py --raw RAW --models per_tissue_models.pkl \\
      --inputs INPUTS --out OUT [--check] [--tissues blood,nasal_epithelium,solid_tissue]
RAW holds the five Atlas exports (original file names), the HELIX
eQTM_autosome_adj.cells_SIG.txt.gz, ann_hg19.bed and the five
EWAS_OpenEWAS_EPICv1_block_* folders; INPUTS is the directory written by
export_inputs.R. Requires lib_eqtm.py (eqtm-direction/analysis) on the path
and scikit-learn 1.7.2 (the version the models were pickled with).
"""
import argparse, glob, os, pickle, sys, time
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--raw", required=True); ap.add_argument("--models", required=True)
ap.add_argument("--inputs", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--lib", default=None, help="directory holding lib_eqtm.py")
ap.add_argument("--check", action="store_true", help="also regenerate the shipped tables and compare")
ap.add_argument("--tissues", default="blood,nasal_epithelium,solid_tissue")
a = ap.parse_args()
if a.lib: sys.path.insert(0, a.lib)
import lib_eqtm as L
import sklearn
assert sklearn.__version__ == "1.7.2", f"the models were pickled with scikit-learn 1.7.2, this is {sklearn.__version__}"
os.makedirs(a.out, exist_ok=True)
T0 = time.time()
def log(*x): print(f"[{time.time() - T0:6.0f}s]", *x, flush=True)

W, NL, NR = 250000, 3, 4
ATLAS = [("allsubjects", "EVA-PR_ All Subjects_data.csv", "nasal_epithelium", "NasalEpi_EVAPR_All", True),
         ("asthmactrl", "EVA-PR_ Atopic asthma + controls_data.csv", "nasal_epithelium", "NasalEpi_EVAPR_AsthmaCtrl", True),
         ("wbc", "EVA-PR_ White Blood Cell_data.csv", "blood", "Blood_EVAPR_WBC", True),
         ("framingham", "Framingham Heart Study_data (1).csv", "blood", "Blood_Framingham", False),
         ("gtex", "GTEx_ WholeBlood_data.csv", "blood", "Blood_GTEx", False)]
HG19 = {"chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
        "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
        "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753, "chr17": 81195210, "chr18": 78077248,
        "chr19": 59128983, "chr20": 63025520, "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566}
nchr = lambda s: "chr" + s.astype(str).str.replace("^chr", "", regex=True).str.replace(r"\.0$", "", regex=True)

# ---- the Atlas exports: first row per gene, cis rows, every gene.chr ------------------
cols = ["probeid", "gene", "TSS", "dist", "beta", "pos", "meth.chr", "gene.chr", "eQTM_label"]
AT = {}
for key, fn, tis, cohort, hg19 in ATLAS:
    first, cis, chrs, seen = [], [], [], set()
    for ch in pd.read_csv(os.path.join(a.raw, fn), usecols=cols, chunksize=1_000_000, low_memory=False):
        ch = ch[cols]
        g = ch.dropna(subset=["gene", "TSS", "pos", "dist"])
        g = g[~g.gene.astype(str).isin(seen)].drop_duplicates("gene"); seen.update(g.gene.astype(str)); first.append(g)
        cis.append(ch[(ch.eQTM_label == "cis") & ch.dist.abs().lt(1_000_000) & ch.beta.notna()])
        chrs.append(ch[["gene", "gene.chr"]].dropna().drop_duplicates())
    AT[key] = dict(first=pd.concat(first), cis=pd.concat(cis), chrs=pd.concat(chrs).drop_duplicates(), tissue=tis, cohort=cohort, hg19=hg19)
    log(f"Atlas {key}: first rows {len(AT[key]['first'])}, cis rows {len(AT[key]['cis'])}")

# ---- training pools, exactly as analysis/02 + 04 build them ------------------------------
def atlas_pool(c):
    return pd.DataFrame({"cpg_id": c.probeid.str.lower(), "gene": c.gene.astype(str).str.upper(),
                         "dist": c.dist.astype(float), "y": (c.beta > 0).astype(int)}).drop_duplicates(["cpg_id", "gene"])
GT = pd.concat([AT[k]["first"].assign(src=k) for k, *_ in ATLAS], ignore_index=True)
GT["g"] = GT.gene.astype(str).str.upper()
GT = GT.drop_duplicates("g").reset_index(drop=True); GT["map_order"] = np.arange(len(GT))
tss_map = GT.set_index("g").TSS

def openewas():
    bed = pd.concat([pd.read_csv(f, sep="\t", header=None, names=["chr", "s", "e", "cpg_id", "t", "st"])
                     for f in sorted(glob.glob(os.path.join(a.raw, "EWAS_*/motif_input.bed")))])
    bed["cpg_id"] = bed.cpg_id.str.lower(); bed["pos"] = (bed.s + bed.e) // 2
    ec = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(a.raw, "EWAS_*/expression_correlation.tsv")))])
    ec = ec[np.isfinite(ec.cor) & (ec.cor != 0)]
    ec["cpg_id"] = ec.probe.str.lower(); ec["g"] = ec.gene_symbol.astype(str).str.upper()
    ec = ec.merge(bed.drop_duplicates("cpg_id")[["cpg_id", "pos"]], on="cpg_id")
    ec["tss"] = ec.g.map(tss_map); ec = ec[ec.tss.notna()].copy()
    ec["dist"] = ec.pos - ec.tss                      # every feature uses |dist|; the sign is irrelevant
    ec = ec[ec.dist.abs() < 1_000_000]
    return {f"TCGA_{t}": pd.DataFrame({"cpg_id": g.cpg_id, "gene": g.g, "dist": g.dist.astype(float),
            "y": (g.cor > 0).astype(int)}).drop_duplicates(["cpg_id", "gene"]) for t, g in ec.groupby("tissue")}

HX = pd.read_csv(os.path.join(a.raw, "eQTM_autosome_adj.cells_SIG.txt.gz"), sep="\t")
C = {AT[k]["cohort"]: atlas_pool(AT[k]["cis"]) for k, *_ in ATLAS}
C["Blood_HELIX"] = L.load_helix(os.path.join(a.raw, "eQTM_autosome_adj.cells_SIG.txt.gz"))
C.update(openewas())
C = L.add_positions(C, os.path.join(a.raw, "ann_hg19.bed"))
MODELS = pickle.load(open(a.models, "rb"))
POOL = {}
for t in ["blood", "nasal_epithelium", "solid_tissue"]:
    P = pd.concat([v for k, v in C.items() if L.tissue_class(k) == t], ignore_index=True)
    POOL[t] = P[P.pos.notna()].reset_index(drop=True)
    assert len(POOL[t]) == MODELS[t]["n"], f"{t}: pool {len(POOL[t])} differs from the model's n {MODELS[t]['n']}"
    log(f"pool {t}: {len(POOL[t])} pairs (= the saved model's n)")

ann = pd.read_csv(os.path.join(a.raw, "ann_hg19.bed"), sep="\t", header=None, comment="t", names=["chr", "s", "e", "name", "sc", "st"])
ann["cpg_id"] = ann.name.astype(str).str.split("_").str[0].str.lower(); ann["pos"] = ann.s + 2
MAN = ann.drop_duplicates("cpg_id")[["cpg_id", "chr", "pos"]].reset_index(drop=True)

# ---- gene sources ------------------------------------------------------------------------
R = pd.concat([AT[k]["chrs"] for k, *_ in ATLAS], ignore_index=True)
R["g"] = R.gene.astype(str).str.upper(); R["c"] = nchr(R["gene.chr"])
atl = R.groupby("g").c.agg(chr_atlas=lambda x: x.mode().iloc[0], n_chr="nunique")
HX["g"] = HX.TC_gene.astype(str).str.upper()
HX["chr_tc"] = "chr" + HX.TC.str.slice(2, 4).astype(int).astype(str).replace({"23": "X", "24": "Y"})
hx = HX.groupby("g").agg(tss_helix=("TC_gene_TSS", "first"), chr_helix=("chr_tc", lambda x: x.mode().iloc[0]))
G0 = GT.set_index("g")[["TSS", "src", "map_order"]].rename(columns={"TSS": "tss_atlas"}).join(atl, how="outer").join(hx, how="outer")
log(f"genes {len(G0)}; Atlas gene.chr inconsistent across rows: {(G0.n_chr > 1).sum()}")

# independent checks
e = pd.read_csv(os.path.join(a.inputs, "epic_gene_chr.csv.gz")).sort_values(["g", "N"], ascending=[True, False])
ep = e.groupby("g").agg(c=("chr", "first"), top=("N", "first"), tot=("N", "sum"))
G0["chr_epic"] = ep[ep.top / ep.tot >= 0.9].c
Mx = pd.read_csv(os.path.join(a.inputs, "measured_pairs.csv.gz"))
G0["chr_meas"] = Mx.groupby("g").chr.agg(lambda x: x.mode().iloc[0])
pp = pd.concat([atlas_pool(AT[k]["cis"])[["cpg_id", "gene"]] for k, *_ in ATLAS] +
               [L.load_helix(os.path.join(a.raw, "eQTM_autosome_adj.cells_SIG.txt.gz"))[["cpg_id", "gene"]]])
pp = pp.merge(MAN[["cpg_id", "chr"]], on="cpg_id")
G0["chr_pool"] = pp.groupby("gene").chr.agg(lambda x: x.mode().iloc[0])
H = pd.read_csv(os.path.join(a.inputs, "hgnc_table.csv.gz"))
H["c"] = H.chromosome.astype(str).str.extract(r"^(\d+|X|Y)\b")[0]; H = H.dropna(subset=["c"]); H["c"] = "chr" + H.c
hg = H.assign(s=H["Approved.Symbol"].str.upper()).groupby("s").c.agg(lambda x: ";".join(sorted(set(x)))).combine_first(
     H.assign(s=H.Symbol.str.upper()).groupby("s").c.agg(lambda x: ";".join(sorted(set(x)))))
def hgnc_chr(g):
    cs = set()
    for part in g.split(";"):
        if part in hg.index: cs.update(hg[part].split(";"))
    return ";".join(sorted(cs)) if cs else np.nan
G0["chr_hgnc"] = [hgnc_chr(g) for g in G0.index]
Z = pd.read_csv(os.path.join(a.inputs, "epic_zhou.csv.gz")).merge(MAN.rename(columns={"chr": "chr19", "pos": "pos19"}), on="cpg_id")
Z = Z[Z.chr == Z.chr19]
Z["a"] = ((Z.pos19 + Z.zhou_dist_tss) / 100).round() * 100; Z["b"] = ((Z.pos19 - Z.zhou_dist_tss) / 100).round() * 100
def zt(d):
    ma = d.a.mode().iloc[0]; fa = (d.a == ma).mean(); mb = d.b.mode().iloc[0]; fb = (d.b == mb).mean()
    m, f = (ma, fa) if fa >= fb else (mb, fb)
    return pd.Series({"tss_zhou": m, "zhou_share": f, "zhou_n": len(d), "chr_zhou": d.chr.iloc[0]})
zz = Z.groupby("g").apply(zt, include_groups=False)
G0 = G0.join(zz[(zz.zhou_n >= 2) & (zz.zhou_share >= 0.5)])

# Atlas TSS carried to hg19 for the non-hg19 exports, through each gene's own cis rows
carried = {}
for k, fn, tis, cohort, hg19 in ATLAS:
    if hg19: continue
    d = AT[k]["cis"].dropna(subset=["TSS", "pos"])
    d = d[nchr(d["meth.chr"]) == nchr(d["gene.chr"])].copy()
    d["cpg_id"] = d.probeid.str.lower(); d = d.merge(MAN.rename(columns={"pos": "pos19"}), on="cpg_id")
    d = d[d.chr == nchr(d["gene.chr"])]
    d["t19"] = d.pos19 - (d.pos - d.TSS); d["g"] = d.gene.astype(str).str.upper()
    s = d.groupby("g").t19.agg(["median", "min", "max"])
    carried[k] = s[(s["max"] - s["min"]) <= 1000]["median"]

# ---- placement ----------------------------------------------------------------------------
pl = []
for g, r in G0.iterrows():
    if pd.isna(r.tss_atlas) and pd.isna(r.tss_helix): continue
    if pd.notna(r.chr_atlas) and r.chr_atlas != "chrnan": c, cs = r.chr_atlas, "Atlas gene.chr"
    elif pd.notna(r.chr_helix): c, cs = r.chr_helix, "HELIX TC"
    elif pd.notna(r.chr_epic) and pd.notna(r.chr_pool) and r.chr_epic == r.chr_pool: c, cs = r.chr_epic, "EPIC v2 + training partners"
    else: c, cs = np.nan, None
    if pd.notna(r.tss_atlas) and r.src in ("allsubjects", "asthmactrl", "wbc"): t, ts = r.tss_atlas, f"Atlas {r.src} (hg19)"
    elif pd.notna(r.tss_atlas):
        if pd.notna(r.tss_zhou) and r.chr_zhou == c: t, ts = r.tss_zhou, "EPIC v2 Zhou GENCODE v41 TSS, carried to hg19"
        elif pd.notna(r.tss_helix): t, ts = r.tss_helix, "HELIX (hg19)"
        elif g in carried.get(r.src, {}): t, ts = carried[r.src][g], f"Atlas {r.src}, carried to hg19"
        else: t, ts = np.nan, None
    else: t, ts = r.tss_helix, "HELIX (hg19)"
    pl.append(dict(g=g, chr=c, tss=t, chr_source=cs, tss_source=ts, first_row=r.src,
                   tss_shipped_source=r.tss_atlas if pd.notna(r.tss_atlas) else r.tss_helix))
PL = pd.DataFrame(pl).set_index("g").join(G0[["chr_atlas", "chr_helix", "chr_epic", "chr_hgnc", "chr_meas", "chr_pool", "tss_zhou", "chr_zhou"]])
PL["in_bounds"] = [(c in HG19) and pd.notna(t) and 0 < t <= HG19[c] for c, t in zip(PL.chr, PL.tss)]
def verdict(r):
    if pd.isna(r.chr) or pd.isna(r.tss) or not r.in_bounds: return "unplaceable"
    ann_ = []
    if pd.notna(r.chr_epic) and r.chr_source != "EPIC v2 + training partners": ann_.append(r.chr_epic == r.chr)
    if pd.notna(r.chr_hgnc): ann_.append(r.chr in r.chr_hgnc.split(";"))
    if ann_: return "confirmed" if all(ann_) else ("conflict" if not any(ann_) else "mixed")
    oth = [r[k] == r.chr for k in ["chr_meas", "chr_pool"] if pd.notna(r[k]) and r.chr_source != "EPIC v2 + training partners"]
    if r.chr_source != "HELIX TC" and pd.notna(r.chr_helix): oth.append(r.chr_helix == r.chr)
    if not oth: return "single_source"
    return "confirmed" if all(oth) else "conflict"
PL["verdict"] = PL.apply(verdict, axis=1)
PL["use"] = PL.verdict.eq("confirmed") | (PL.verdict.eq("mixed") & (PL.chr_epic == PL.chr))
PL.to_csv(os.path.join(a.out, "gene_placement.csv.gz"))
log("placement verdicts:", PL.verdict.value_counts().to_dict(), "; used:", int(PL.use.sum()))

# ---- one lookup ----------------------------------------------------------------------------
def build(t, G):
    P, MOD = POOL[t], MODELS[t]
    G = G.copy(); G["mo"] = G.g.map(GT.set_index("g").map_order).fillna(1e9)
    parts = []
    for ch, m in MAN.groupby("chr"):
        g = G[G.chr == ch].sort_values(["tss", "mo", "g"], kind="stable")
        if not len(g): continue
        Tt = g.tss.values.astype(float); nm = g.g.values; Pp = m.pos.values.astype(float); ids = m.cpg_id.values
        ins = np.searchsorted(Tt, Pp, side="left")
        for off in range(-NL, NR):
            j = ins + off; ok = (j >= 0) & (j < len(Tt)); jj = j[ok]; d = Pp[ok] - Tt[jj]; keep = np.abs(d) <= W
            parts.append(pd.DataFrame({"cpg_id": ids[ok][keep], "target_gene": nm[jj][keep], "tss_dist": d[keep],
                                       "chr": ch, "pos": Pp[ok][keep]}))
    Q = pd.concat(parts, ignore_index=True)
    q = pd.DataFrame({"cpg_id": Q.cpg_id.values, "gene": Q.target_gene.str.upper().values, "dist": Q.tss_dist.values,
                      "chr": Q.chr.values, "pos": Q.pos.values})
    qg = q.copy(); qg["cpg_id"] = qg.cpg_id + "|q"
    X = pd.concat([L.full_features(q, P)[["abs", "log", "p500", "p2k", "p10k", "p100k"]],
                   L.neighbour_prior(q, P), L.gene_prior(qg, P)], axis=1)
    p = MOD["calibrator"].predict(MOD["model"].predict_proba(X[MOD["features"]].values)[:, 1])
    near = X.nn50000.values > 0; gene = X.gn.values > 0
    tier = np.where(near & gene, "A", np.where(near | gene, "B", "C"))
    conf = np.abs(p - 0.5) * 2
    call = ((tier == "A") & (conf >= 0.30)) | ((tier == "B") & (conf >= 0.50))
    ym = P.groupby(["cpg_id", "gene"]).y.mean().reindex(pd.MultiIndex.from_arrays([q.cpg_id, q.gene])).values
    direct = ~np.isnan(ym)
    status = np.where(direct, "DIRECT_eQTM", np.where(call, "PREDICTED", "ABSTAIN"))
    direction = np.where(direct, np.where(ym >= 0.5, 1.0, -1.0), np.where(call, np.where(p >= 0.5, 1.0, -1.0), np.nan))
    return pd.DataFrame({"cpg_id": Q.cpg_id.values, "target_gene": Q.target_gene.values, "tss_dist": Q.tss_dist.astype(int).values,
                         "status": status, "direction": direction, "probability_plus1": np.round(p, 4),
                         "confidence": np.round(conf, 4), "evidence_tier": tier})

def compare(Rb, S):
    m = Rb.merge(S, on=["cpg_id", "target_gene"], how="outer", suffixes=("", "_s"), indicator=True)
    b = m[m._merge == "both"]
    same = lambda c: float(np.mean(b[c].astype(str).values == b[c + "_s"].astype(str).values))
    num = lambda c: float(np.mean(np.abs(b[c] - b[c + "_s"]) <= 1e-4))
    allsame = ((b.tss_dist == b.tss_dist_s) & (b.status == b.status_s) & (b.evidence_tier == b.evidence_tier_s) &
               (np.abs(b.probability_plus1 - b.probability_plus1_s) <= 1e-4) & (np.abs(b.confidence - b.confidence_s) <= 1e-4) &
               (b.direction.fillna(0) == b.direction_s.fillna(0)))
    return dict(rows=len(Rb), shipped=len(S), both=len(b), identical=int(allsame.sum()), new=int((m._merge == "left_only").sum()),
                dropped=int((m._merge == "right_only").sum()), tss_dist=same("tss_dist"), status=same("status"),
                tier=same("evidence_tier"), probability=num("probability_plus1"), confidence=num("confidence"),
                direction=float(np.mean(b.direction.fillna(0).values == b.direction_s.fillna(0).values)))

report = []
for t in a.tissues.split(","):
    pool_genes = set(POOL[t].gene.unique())
    S = pd.read_csv(os.path.join(a.inputs, f"lookup_{t}.csv.gz"))
    if a.check:
        s = S.merge(MAN, on="cpg_id", how="left")
        Gs = s.assign(g=s.target_gene.astype(str), tss=s.pos - s.tss_dist).groupby("g").agg(chr=("chr", "first"), tss=("tss", "first")).reset_index()
        rc = compare(build(t, Gs), S); report.append(("check", t, rc)); log(f"{t} regression (shipped placement):", rc)
    G = PL[PL.use & PL.index.isin(pool_genes)].reset_index()[["g", "chr", "tss"]]
    Rb = build(t, G)
    Rb.to_csv(os.path.join(a.out, f"lookup_{t}.csv.gz"), index=False)
    rf = compare(Rb, S)
    rf.update(pool_genes=len(pool_genes), placed=len(G), cpgs=int(Rb.cpg_id.nunique()), genes=int(Rb.target_gene.nunique()),
              status_counts=Rb.status.value_counts().to_dict())
    report.append(("rebuilt", t, rf)); log(f"{t} rebuilt:", rf)
pd.DataFrame([dict(kind=k, tissue=t, **{x: y for x, y in r.items() if not isinstance(y, dict)}) for k, t, r in report]).to_csv(
    os.path.join(a.out, "build_report.csv"), index=False)
log("done")
