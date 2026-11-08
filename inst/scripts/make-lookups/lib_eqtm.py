"""
lib_eqtm.py — shared functions for the eQTM direction analysis.

Coordinate handling is the part to read carefully. The twelve source catalogues
are in different genome builds, use opposite and non-strand-aware sign
conventions for their distance field, and disagree about where genes start.
Everything here therefore uses |distance to TSS|, which is invariant to all
three. See docs/AUDIT_REPORT_coordinates.md.
"""
from __future__ import annotations
import glob, os
import numpy as np, pandas as pd

RAW = os.environ.get("EQTM_RAW", "../data/raw")

ATLAS = {
    "EVA-PR_ All Subjects_data.csv":            "NasalEpi_EVAPR_All",
    "EVA-PR_ Atopic asthma + controls_data.csv":"NasalEpi_EVAPR_AsthmaCtrl",
    "EVA-PR_ White Blood Cell_data.csv":        "Blood_EVAPR_WBC",
    "Framingham Heart Study_data (1).csv":      "Blood_Framingham",
    "GTEx_ WholeBlood_data.csv":                "Blood_GTEx",
}

def tissue_class(cohort: str) -> str:
    if cohort.startswith("Blood"):    return "blood"
    if cohort.startswith("NasalEpi"): return "nasal_epithelium"
    return "solid_tissue"

# ---------------------------------------------------------------- cohorts
def load_atlas(raw=RAW) -> dict:
    """eQTM Atlas cohorts. Uses each cohort's own distance field, absolute value only."""
    out = {}
    for fn, key in ATLAS.items():
        p = os.path.join(raw, fn)
        if not os.path.exists(p):
            print(f"[skip] {fn} not found"); continue
        d = pd.read_csv(p, usecols=["probeid","gene","dist","beta","eQTM_label"])
        d = d[(d.eQTM_label == "cis") & d.dist.abs().lt(1_000_000) & d.beta.notna()]
        out[key] = pd.DataFrame({
            "cpg_id": d.probeid.str.lower(),
            "gene":   d.gene.astype(str).str.upper(),
            "dist":   d.dist.astype(float),
            "y":      (d.beta > 0).astype(int),
        }).drop_duplicates(["cpg_id","gene"])
    return out

def load_helix(path=None) -> pd.DataFrame:
    """HELIX from the original release; the Atlas copy has beta = NA and an unsigned dist."""
    path = path or os.path.join(RAW, "eQTM_autosome_adj.cells_SIG.txt.gz")
    h = pd.read_csv(path, sep="\t")
    strand = np.where(h.TC_gene_TSS == h.TC_gene_start, 1, -1)
    return pd.DataFrame({
        "cpg_id": h.CpG.str.lower(),
        "gene":   h.TC_gene.astype(str).str.upper(),
        "dist":   (h.CpG_pos - h.TC_gene_TSS) * strand,
        "y":      (h.log2FC > 0).astype(int),
    }).drop_duplicates(["cpg_id","gene"])

def gene_tss_map(raw=RAW):
    """Gene -> TSS and strand, built from the Atlas cohorts, which ship a genuine TSS
    column. Required for the TCGA sets, which do not. An earlier attempt used an SMR
    .epi file whose fourth column is the probe midpoint, not the TSS; see the audit."""
    G, S = [], []
    for fn in ATLAS:
        p = os.path.join(raw, fn)
        if not os.path.exists(p): continue
        d = pd.read_csv(p, usecols=["gene","TSS","pos","dist"]).dropna()
        g = d.gene.astype(str).str.upper()
        G.append(pd.DataFrame({"g": g, "tss": d.TSS.astype(float)}))
        S.append(pd.DataFrame({"g": g, "s": np.sign(d.dist) * np.sign(d.pos - d.TSS)}))
    return (pd.concat(G).drop_duplicates("g").set_index("g").tss,
            pd.concat(S).groupby("g").s.median())

def load_openewas(raw=RAW) -> dict:
    """TCGA analyses via OpenEWAS. Probe positions come from the block motif_input.bed
    files; gene TSS from gene_tss_map(). Validate with validate_tcga_mapping()."""
    tss, strand = gene_tss_map(raw)
    bed = pd.concat([pd.read_csv(f, sep="\t", header=None,
                     names=["chr","s","e","cpg_id","t","st"])
                     for f in glob.glob(os.path.join(raw, "EWAS_*/motif_input.bed"))])
    bed["cpg_id"] = bed.cpg_id.str.lower(); bed["pos"] = (bed.s + bed.e) // 2
    ec = pd.concat([pd.read_csv(f, sep="\t")
                    for f in glob.glob(os.path.join(raw, "EWAS_*/expression_correlation.tsv"))])
    ec = ec[np.isfinite(ec.cor) & (ec.cor != 0)]
    ec["cpg_id"] = ec.probe.str.lower(); ec["g"] = ec.gene_symbol.astype(str).str.upper()
    ec = ec.merge(bed.drop_duplicates("cpg_id")[["cpg_id","pos"]], on="cpg_id")
    ec["tss"] = ec.g.map(tss); ec["strand"] = ec.g.map(strand).fillna(1)
    ec = ec[ec.tss.notna()].copy()
    ec["dist"] = (ec.pos - ec.tss) * np.where(ec.strand >= 0, 1, -1)
    ec = ec[ec.dist.abs() < 1_000_000]
    return {f"TCGA_{t}": pd.DataFrame({"cpg_id": g.cpg_id, "gene": g.g,
            "dist": g.dist.astype(float), "y": (g.cor > 0).astype(int)}
            ).drop_duplicates(["cpg_id","gene"]) for t, g in ec.groupby("tissue")}, ec

def validate_tcga_mapping(ec: pd.DataFrame) -> pd.DataFrame:
    """Gate on this. Promoter-labelled rows must sit near zero. Published run:
    promoter median +10 bp, 88.1% within 2 kb; body 30.7% within 2 kb."""
    return ec.groupby("location").dist.agg(
        n="size", median="median", within_2kb=lambda s: (s.abs() < 2000).mean())

def load_all(raw=RAW) -> dict:
    C = load_atlas(raw)
    C["Blood_HELIX"] = load_helix()
    tcga, ec = load_openewas(raw)
    v = validate_tcga_mapping(ec)
    prom = float(v.loc["promoter","within_2kb"]) if "promoter" in v.index else np.nan
    if not (prom > 0.80):
        raise RuntimeError(f"TCGA coordinate validation failed: only {prom:.1%} of "
                           "promoter-labelled rows within 2 kb of the computed TSS.")
    C.update(tcga)
    return C

# ---------------------------------------------------------------- features
def distance_features(d: pd.DataFrame) -> np.ndarray:
    """Strand-, sign-convention- and build-invariant. This is deliberate."""
    a = d.dist.abs()
    return pd.DataFrame({
        "abs":  a, "log": np.log10(a + 1),
        "p500": (a < 500).astype(int),   "p2k":  (a < 2000).astype(int),
        "p10k": (a < 10000).astype(int), "p100k":(a < 100000).astype(int),
    }).values

def neighbour_prior(q, ref, wins=(2000, 10000, 50000)) -> pd.DataFrame:
    """Mean direction of catalogued CpGs within each window, target CpG excluded.
    Requires chr/pos columns; take them from one manifest for all cohorts."""
    out = {}
    for w in wins:
        out[f"nb{w}"] = np.zeros(len(q)); out[f"nn{w}"] = np.zeros(len(q))
    lab = (ref.dropna(subset=["pos"]).groupby(["chr","cpg_id"])
              .agg(pos=("pos","first"), m=("y","mean")).reset_index())
    lab["s"] = 2 * lab.m - 1
    for ch, g in lab.groupby("chr"):
        g = g.sort_values("pos")
        P, S, ID = g.pos.values.astype(float), g.s.values, g.cpg_id.values
        cs = np.concatenate([[0.], np.cumsum(S)]); idp = {x: i for i, x in enumerate(ID)}
        sel = np.where(q.chr.values == ch)[0]
        if not len(sel): continue
        qp = q.pos.values[sel].astype(float)
        si = np.array([idp.get(x, -1) for x in q.cpg_id.values[sel]])
        for w in wins:
            lo = np.searchsorted(P, qp - w); hi = np.searchsorted(P, qp + w)
            ss = cs[hi] - cs[lo]; n = (hi - lo).astype(float)
            has = (si >= 0) & (si >= lo) & (si < hi)
            ss -= np.where(has, S[np.clip(si, 0, len(S) - 1)], 0.); n -= has.astype(float)
            out[f"nb{w}"][sel] = np.where(n > 0, ss / np.maximum(n, 1), 0.)
            out[f"nn{w}"][sel] = n
    return pd.DataFrame(out, index=q.index)

def gene_prior(q, ref) -> pd.DataFrame:
    """Mean direction of other CpGs for the same gene, all pairs of the target CpG excluded."""
    g  = ref.groupby("gene").y.agg(["sum","size"])
    cg = ref.groupby(["gene","cpg_id"]).y.agg(["sum","size"])
    tot = g.reindex(q.gene.values).fillna(0)
    own = cg.reindex(pd.MultiIndex.from_arrays([q.gene, q.cpg_id])).fillna(0)
    s = tot["sum"].values - own["sum"].values
    n = tot["size"].values - own["size"].values
    return pd.DataFrame({"gd": np.where(n > 0, 2 * (s / np.maximum(n, 1)) - 1, 0.),
                         "gn": n}, index=q.index)

def full_features(q, ref) -> pd.DataFrame:
    a = q.dist.abs()
    X = pd.DataFrame({"abs": a, "log": np.log10(a + 1),
                      "p500": (a < 500).astype(int),   "p2k":  (a < 2000).astype(int),
                      "p10k": (a < 10000).astype(int), "p100k":(a < 100000).astype(int)},
                     index=q.index)
    return pd.concat([X, neighbour_prior(q, ref), gene_prior(q, ref)], axis=1)

def add_positions(C: dict, manifest_bed: str) -> dict:
    """Attach hg19 probe positions from ONE manifest for every cohort, which removes
    the build heterogeneity from feature construction."""
    ann = pd.read_csv(manifest_bed, sep="\t", header=None, comment="t",
                      names=["chr","s","e","name","sc","st"])
    ann["cpg_id"] = ann.name.astype(str).str.split("_").str[0].str.lower()
    ann["pos"] = ann.s + 2
    A = ann.drop_duplicates("cpg_id").set_index("cpg_id")[["chr","pos"]]
    for k in C:
        C[k] = C[k].join(A, on="cpg_id")
        C[k]["cohort"] = k; C[k]["tissue"] = tissue_class(k)
    return C

# ---------------------------------------------------------------- metrics
def wilson(s, n, z=1.96):
    s = np.asarray(s, float); n = np.maximum(np.asarray(n, float), 1)
    p = s / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h

def excess_concordance(A: pd.DataFrame, B: pd.DataFrame, min_pairs=30):
    """Observed agreement minus agreement expected from each cohort's own base rates.
    Two catalogues each 60% negative agree 52% of the time knowing nothing."""
    from scipy import stats
    m = A.merge(B, on=["cpg_id","gene"], suffixes=("_a","_b"))
    if len(m) < min_pairs: return None
    ya = 2 * m.y_a - 1; yb = 2 * m.y_b - 1
    agree = (ya == yb)
    pa = (ya > 0).mean(); pb = (yb > 0).mean()
    chance = pa * pb + (1 - pa) * (1 - pb)
    k, n = int(agree.sum()), len(m)
    return dict(n=n, agree=k / n, chance=chance, excess=k / n - chance,
                p=stats.binomtest(k, n, chance, alternative="greater").pvalue)

def gbm(seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(max_iter=400, learning_rate=0.05, max_depth=5,
                                          l2_regularization=1.0, random_state=seed)
