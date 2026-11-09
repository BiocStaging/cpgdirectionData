#!/usr/bin/env python3
"""Step 1 of the v2 catalogue build: one gene identity, one TSS, harmonized training pools.

The training pools are the original ones (analysis/02 + 04 of eqtm-direction): the same
cohort x CpG x gene rows, the same labels. What changes:

  identity   every gene name of every cohort is resolved to one Ensembl GRCh37 gene
             (annotables grch37, primary chromosomes), on the chromosome its own pool
             CpGs sit on:
               by ENSG ID, by symbol, or by HGNC previous/alias symbol (HGNChelper),
               accepted when the cohort's own hg19 coordinate for the gene (EVA-PR or
               HELIX TSS; Framingham's leftmost coordinate or GTEx's TSS, carried to hg19
               through the gene's own cis rows) lies within the gene +- 10 kb (several such
               genes: the one whose TSS lies within 10 bp of that coordinate);
               else, for an exact symbol carried by several Ensembl IDs there (duplicate
               annotations), the lowest ID;
               else an exact, non-clone symbol, the only gene of that name on the
               chromosome, whose coordinate lies within 100 kb of the gene (an upstream
               first exon Ensembl GRCh37 lacks);
               else by position: the one Ensembl gene whose TSS (Framingham: start) lies
               within 10 bp of that coordinate;
               else the name is not identified and its rows are dropped.
             HELIX rows whose transcript cluster names no single identified gene (blank
             or multi-gene TC_gene) keep the cluster as their identity, with HELIX's own
             strand-aware TSS. Names of different cohorts that resolve to one gene are
             one gene.
  TSS        one 1-based, strand-aware TSS per gene (the strand is the Ensembl gene's):
               --tss cohort (default; the shipped build): the cohort's own hg19 TSS where
                 it has the gene and the coordinate agrees with the identity (on the gene's
                 chromosome, within the gene +- 10 kb, or +- 100 kb for a name identified at
                 up to 100 kb): the EVA-PR first-row TSS (0-based on the plus strand there,
                 so +1), else the HELIX TC TSS; else the Ensembl GRCh37 gene TSS (start on
                 the plus strand, end on the minus strand);
               --tss ensembl: the Ensembl TSS for every gene (validated worse on the
                 solid yardsticks; not shipped).
             HELIX clusters take the HELIX TC TSS.
  distance   every training row's dist is recomputed as (hg19 CpG position - TSS) x
             strand, from the same TSS the lookups use. The cohorts' own distance
             fields are not used: Framingham's runs to the gene's leftmost coordinate,
             and TCGA's came from the first Atlas row's TSS, not hg19 for genes first
             seen in Framingham or GTEx.
  duplicates within a cohort, rows of one CpG and one gene (two names of one gene) are
             one row: kept when their labels agree, dropped when they disagree.

Usage:
  python build_pools.py --raw RAW --inputs INPUTS --out OUT [--lib DIR] [--tss cohort|ensembl]
RAW holds the five Atlas exports (original file names), the HELIX
eQTM_autosome_adj.cells_SIG.txt.gz, ann_hg19.bed and the five EWAS_OpenEWAS_EPICv1_block_*
folders; INPUTS is the directory written by export_inputs.R (ens_grch37.csv.gz,
hgnc_table.csv.gz). Writes OUT/pool_<tissue>.csv.gz (cpg_id, gene_id, gene = display name,
dist, y, chr, pos, cohort, tissue, name = the cohort's own name), OUT/gene_table.csv.gz (one
row per Ensembl gene: chromosome, TSS and its source, strand, display name, every name
(separated by " | "), the tie order of the candidate rule), OUT/name_identity.csv.gz (every pool name and how it was
resolved) and OUT/pools_report.csv. Deterministic.
Display name: a cohort's own name, preferring EVA-PR's (for continuity with earlier lookups),
never a name identified by position that is, or contains, another Ensembl gene's symbol (with no
name left, the Ensembl symbol).
"""
import argparse, glob, os, re, sys, time
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--raw", required=True); ap.add_argument("--inputs", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--lib", default=None, help="directory holding lib_eqtm.py")
ap.add_argument("--tss", choices=["ensembl", "cohort"], default="cohort", help="TSS source (see header)")
a = ap.parse_args()
if a.lib: sys.path.insert(0, a.lib)
import lib_eqtm as L
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
POOL = {}
for t in ["blood", "nasal_epithelium", "solid_tissue"]:
    P = pd.concat([v for k, v in C.items() if L.tissue_class(k) == t], ignore_index=True)
    POOL[t] = P[P.pos.notna()].reset_index(drop=True)
    log(f"pool {t}: {len(POOL[t])} pairs")

ann = pd.read_csv(os.path.join(a.raw, "ann_hg19.bed"), sep="\t", header=None, comment="t", names=["chr", "s", "e", "name", "sc", "st"])
ann["cpg_id"] = ann.name.astype(str).str.split("_").str[0].str.lower(); ann["pos"] = ann.s + 2
MAN = ann.drop_duplicates("cpg_id")[["cpg_id", "chr", "pos"]].reset_index(drop=True)


# ---- reference tables ------------------------------------------------------------------
PRIM = set(HG19)
EN = pd.read_csv(os.path.join(a.inputs, "ens_grch37.csv.gz"))
EN = EN[EN.chr.isin(PRIM)].drop_duplicates("ensgene").reset_index(drop=True)
EN["sym"] = EN.symbol.fillna("").astype(str).str.upper()
HG = pd.read_csv(os.path.join(a.inputs, "hgnc_table.csv.gz"))
HG["c"] = "chr" + HG.chromosome.astype(str).str.extract(r"^(\d+|X|Y)")[0]
HG["s"] = HG.Symbol.astype(str).str.upper(); HG["a"] = HG["Approved.Symbol"].astype(str).str.upper()
APPROVED = HG.dropna(subset=["c"]).groupby(["s", "c"]).a.agg(lambda x: set(x)).to_dict()
def name_keys(n, c): return {n} | APPROVED.get((n, c), set())
KEY = {}
for i, (s, c) in enumerate(zip(EN.sym, EN.chr)):
    if s:
        for k in name_keys(s, c): KEY.setdefault((k, c), set()).add(i)
BYCHR = {c: x for c, x in EN.groupby("chr")}
CLONE = re.compile(r"^(A[A-Z]\d{6}\.\d+|RP\d+-|CT[ABCD]-|CH\d+-|XX[A-Za-z]*-)")

# ---- each name's own hg19 coordinates ----------------------------------------------------
GT_ = GT.set_index("g")
HXg = HX.assign(g=HX.TC_gene.astype(str).str.upper()).groupby("g").TC_gene_TSS.first()
carried = {}
for k, fn, tis, cohort, hg19 in ATLAS:
    if hg19: continue
    d = AT[k]["cis"].dropna(subset=["TSS", "pos"]).copy()
    mc = nchr(d["meth.chr"]); gc = nchr(d["gene.chr"]).where(d["gene.chr"].notna(), mc)
    d = d.assign(gchr=gc)[mc == gc]
    d["cpg_id"] = d.probeid.str.lower(); d = d.merge(MAN.rename(columns={"pos": "pos19"}), on="cpg_id")
    d = d[d.chr == d.gchr]
    d["f"] = d.pos19 - (d.pos - d.TSS); d["r"] = d.pos19 + (d.pos - d.TSS); d["g"] = d.gene.astype(str).str.upper()
    s = d.groupby("g").agg(f=("f", "median"), f_lo=("f", "min"), f_hi=("f", "max"), r=("r", "median"), r_lo=("r", "min"), r_hi=("r", "max"))
    s["f_ok"] = (s.f_hi - s.f_lo) <= 1000; s["r_ok"] = (s.r_hi - s.r_lo) <= 1000
    carried[k] = s
def own_coords(g):
    """(coordinate, kind) pairs: kind 'tss' compares with an Ensembl TSS, 'start' with a gene start."""
    out = []
    src = GT_.src.get(g) if g in GT_.index else None
    if src in ("allsubjects", "asthmactrl", "wbc"): out.append((float(GT_.at[g, "TSS"]), "tss"))
    if g in HXg.index and pd.notna(HXg[g]): out.append((float(HXg[g]), "tss"))
    if src == "framingham" and g in carried["framingham"].index and carried["framingham"].at[g, "f_ok"]:
        out.append((float(carried["framingham"].at[g, "f"]), "start"))
    if src == "gtex" and g in carried["gtex"].index:
        cr = carried["gtex"].loc[g]
        if cr.f_ok: out.append((float(cr.f), "tss"))
        if cr.r_ok: out.append((float(cr.r), "tss"))
    return out

# ---- identity of every pool name -----------------------------------------------------
allp = pd.concat(POOL.values(), ignore_index=True)
gchr = allp.groupby("gene").chr.agg(lambda x: x.mode().iloc[0])
rows = []
for g, c in gchr.items():
    co = own_coords(g); e = BYCHR.get(c)
    gb = re.sub(r"\.\d+$", "", g)
    if gb.startswith("ENSG"): cand, how = set(EN.index[(EN.ensgene == gb) & (EN.chr == c)]), "ENSG ID"
    else:
        cand = {i for i in KEY.get((g, c), set()) if EN.at[i, "sym"] == g}; how = "symbol"
        if not cand: cand, how = set().union(*[KEY.get((k, c), set()) for k in name_keys(g, c)]), "HGNC previous/alias symbol"
    def consistent(i):
        return not co or any(EN.at[i, "start"] - 10000 <= x <= EN.at[i, "end"] + 10000 for x, _ in co)
    ok = {i for i in cand if consistent(i)}
    pos = set()
    if e is not None:
        for x, kind in co:
            col = "tss" if kind == "tss" else "start"
            pos |= set(e.index[(e[col] - x).abs() <= 10])
    # an exact symbol carried by several Ensembl GRCh37 IDs on the chromosome (duplicate
    # annotations, e.g. CFB, KMT2B, HOXC4): the lowest ID, the long-standing one
    dup = {i for i in ok if EN.at[i, "sym"] == g} if how == "symbol" else set()
    # an exact, non-clone symbol whose cohort coordinate lies upstream of Ensembl's gene span
    # (an alternative first exon Ensembl GRCh37 does not have, e.g. TLR2, NRXN3, TP63)
    def near(i, w=100000):
        return any(EN.at[i, "start"] - w <= x <= EN.at[i, "end"] + w for x, _ in co)
    far = [i for i in cand if EN.at[i, "sym"] == g and near(i)] if (how == "symbol" and not ok and not CLONE.match(g)) else []
    if len(ok) == 1: idx, rule = next(iter(ok)), how
    elif len(ok) > 1 and len(ok & pos) == 1: idx, rule = next(iter(ok & pos)), how + " + position"
    elif len(dup) > 1 and dup == ok: idx, rule = min(dup, key=lambda i: EN.at[i, "ensgene"]), "symbol (lowest of several Ensembl IDs)"
    elif len(far) == 1 and len(cand) == 1: idx, rule = far[0], "symbol (cohort TSS up to 100 kb outside the gene)"
    elif len(pos) == 1: idx, rule = next(iter(pos)), ("position (name points elsewhere)" if cand else "position")
    else: idx, rule = None, ("ambiguous" if (len(ok) > 1 or len(pos) > 1) else ("name inconsistent with own coordinates" if cand else "not identified"))
    rows.append(dict(name=g, chr=c, ensgene=EN.at[idx, "ensgene"] if idx is not None else None, rule=rule,
                     n_name=len(cand), n_consistent=len(ok), n_position=len(pos), n_own_coords=len(co), clone=bool(CLONE.match(g))))
NI = pd.DataFrame(rows)
log("name identity:", NI.rule.value_counts().to_dict())

# ---- HELIX rows: the transcript cluster of each pool row ------------------------------
h = HX.copy(); h["cpg_id"] = h.CpG.str.lower(); h["gene"] = h.TC_gene.astype(str).str.upper()
h = h.drop_duplicates(["cpg_id", "gene"])        # the rows load_helix keeps, in its order
h["hx_strand"] = np.where(h.TC_gene_TSS == h.TC_gene_start, 1, -1)
h["hx_chr"] = "chr" + h.TC.str.slice(2, 4).astype(int).astype(str).replace({"23": "X", "24": "Y"})
HXrow = h.set_index(["cpg_id", "gene"])[["TC", "TC_gene_TSS", "hx_strand", "hx_chr"]]

# ---- genes: identity, TSS, display name ---------------------------------------------
idmap = NI.set_index("name").ensgene
G = EN.set_index("ensgene").loc[sorted(set(NI.ensgene.dropna()))][["chr", "tss", "strand", "symbol", "start", "end"]]
names = NI.dropna(subset=["ensgene"]).groupby("ensgene").name.agg(list)
order = {g: i for i, g in enumerate(GT.g)}                      # the old map order: EVA-PR, then Framingham, GTEx
RULE = NI.dropna(subset=["ensgene"]).set_index(["name", "ensgene"]).rule.to_dict()
OWNER = EN.groupby("sym").ensgene.agg(set).to_dict()
def display(eid):
    # a name identified by position that is another Ensembl gene's symbol (a mislabeled or a
    # reused clone name, e.g. HELIX "IFI6" at TRBV6-1) is not shown; with none left, the Ensembl symbol
    other = lambda n: any(OWNER.get(x, set()) - {eid} for x in n.split(";"))
    ns = [n for n in names[eid] if not (RULE.get((n, eid), "").startswith("position") and other(n))]
    if not ns: return str(G.at[eid, "symbol"]).upper()
    pref = sorted(ns, key=lambda n: (bool(CLONE.match(n)) or n.startswith("ENSG"), order.get(n, 10**9), n))
    return pref[0]
G["display"] = [display(i) for i in G.index]
# one TSS per gene. --tss ensembl: the Ensembl GRCh37 gene TSS. --tss cohort: the cohorts' own hg19
# TSS where they have the gene (EVA-PR first row, 0-based on the plus strand there, so +1; else the
# HELIX TC TSS), else the Ensembl TSS. Both strand-aware and 1-based. A cohort TSS is used only
# where it agrees with the identity: on the gene's chromosome, within the gene +- 10 kb (+- 100 kb
# for a name identified at up to 100 kb).
G = G.rename(columns={"tss": "tss_ensembl"})
HXc = HX.assign(g=HX.TC_gene.astype(str).str.upper()).groupby("g").TC.first().str.slice(2, 4).astype(int).astype(str).replace({"23": "X", "24": "Y"})
def agrees(eid, n, x, c):
    """the cohort coordinate x (chromosome c) of name n lies where n's identity accepts it"""
    w = 100000 if "100 kb" in RULE.get((n, eid), "") else 10000
    return c == G.at[eid, "chr"] and G.at[eid, "start"] - w <= x <= G.at[eid, "end"] + w
def cohort_tss(eid):
    st = G.at[eid, "strand"]
    eva = [n for n in names[eid] if n in GT_.index and GT_.at[n, "src"] in ("allsubjects", "asthmactrl", "wbc")
           and agrees(eid, n, float(GT_.at[n, "TSS"]), "chr" + str(GT_.at[n, "gene.chr"]).replace("chr", "").replace(".0", ""))]
    if eva:
        n = min(eva, key=lambda n: order.get(n, 10**9))
        return float(GT_.at[n, "TSS"]) + (1 if st == 1 else 0), f"EVA-PR {GT_.at[n, 'src']}"
    hx = [n for n in names[eid] if n in HXg.index and pd.notna(HXg[n]) and agrees(eid, n, float(HXg[n]), "chr" + HXc.get(n, ""))]
    if hx: return float(HXg[sorted(hx)[0]]), "HELIX"
    return float(G.at[eid, "tss_ensembl"]), "Ensembl GRCh37"
if a.tss == "cohort":
    ct = [cohort_tss(i) for i in G.index]; G["tss"] = [x for x, _ in ct]; G["tss_source"] = [y for _, y in ct]
else:
    G["tss"] = G.tss_ensembl; G["tss_source"] = "Ensembl GRCh37"
log("TSS source:", G.tss_source.value_counts().to_dict())
G["names"] = [" | ".join(sorted(names[i])) for i in G.index]   # some cohort names contain ";"
G["map_order"] = [min(order.get(n, 10**9) for n in names[i]) for i in G.index]   # tie order of the candidate rule
assert G.display.is_unique, G[G.display.duplicated(keep=False)]
log(f"genes: {len(G)} Ensembl genes from {int(NI.ensgene.notna().sum())} names; {int((names.map(len) > 1).sum())} genes carry more than one name")

# ---- harmonized pools ----------------------------------------------------------------------
rep = []
OUTP = {}
for t, P in POOL.items():
    P = P.copy(); n0 = len(P); P["gene_name"] = P.gene
    P["gene_id"] = P.gene.map(idmap)
    hx = P.cohort.eq("Blood_HELIX")
    if hx.any():
        k = pd.MultiIndex.from_arrays([P.cpg_id[hx], P.gene[hx]])
        hr = HXrow.reindex(k); assert hr.TC.notna().all(), "HELIX pool rows without a cluster"
        cl = P.loc[hx, "gene_id"].isna().values
        P.loc[hx & P.gene_id.isna(), "gene_id"] = ("HELIX:" + hr.TC.values)[cl]
        P["tss"] = P.gene_id.map(G.tss); P["strand"] = P.gene_id.map(G.strand); P["gchr"] = P.gene_id.map(G.chr)
        hcl = hx & P.gene_id.astype(str).str.startswith("HELIX:")
        P.loc[hcl, "tss"] = hr.TC_gene_TSS.values[cl]; P.loc[hcl, "strand"] = hr.hx_strand.values[cl]; P.loc[hcl, "gchr"] = hr.hx_chr.values[cl]
    else:
        P["tss"] = P.gene_id.map(G.tss); P["strand"] = P.gene_id.map(G.strand); P["gchr"] = P.gene_id.map(G.chr)
    unid = P.gene_id.isna()
    xchr = ~unid & (P.gchr != P.chr)
    P = P[~unid & ~xchr].copy()
    P["dist"] = (P.pos - P.tss) * P.strand
    far = P.dist.abs() >= 1_000_000
    P = P[~far]
    g = P.groupby(["cohort", "cpg_id", "gene_id"]).y
    nlab = g.transform("nunique"); n1 = len(P)
    P = P[nlab == 1].drop_duplicates(["cohort", "cpg_id", "gene_id"]); ndup = n1 - len(P)
    P["gene"] = P.gene_id.map(G.display).fillna(P.gene_id)
    P["name"] = P.pop("gene_name") if "gene_name" in P else P.gene
    P = P[["cpg_id", "gene_id", "gene", "dist", "y", "chr", "pos", "cohort", "tissue", "name"]].reset_index(drop=True)
    OUTP[t] = P
    rep.append(dict(tissue=t, rows_original=n0, not_identified=int(unid.sum()), other_chromosome=int(xchr.sum()), beyond_1Mb=int(far.sum()),
                    duplicate_rows_dropped=int(ndup), rows=len(P), genes=P.gene_id.nunique(),
                    helix_clusters=int(P.gene_id.str.startswith("HELIX:").sum())))
    P.to_csv(os.path.join(a.out, f"pool_{t}.csv.gz"), index=False)
    log(f"pool {t}:", rep[-1])
pd.DataFrame(rep).to_csv(os.path.join(a.out, "pools_report.csv"), index=False)
G.reset_index().rename(columns={"index": "gene_id", "ensgene": "gene_id"}).to_csv(os.path.join(a.out, "gene_table.csv.gz"), index=False)
NI.to_csv(os.path.join(a.out, "name_identity.csv.gz"), index=False)
log("done")
