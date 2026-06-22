"""
run_coloc.py — drug-target colocalisation pipeline (PPi-TNAP calcification axis)
================================================================================
For each (gene, tissue) exposure x disease outcome pair:
  1. read the eQTL-Catalogue cis region file, KEEP ONLY the target gene_id
     (region files contain every gene in the locus; failing to filter dilutes
      the signal and inflates PP.H3 — this filter is mandatory);
  2. read the outcome GWAS region file (FinnGen or GWAS-Catalogue harmonised);
  3. match SNPs by chr:pos, harmonise to a common effect allele (flip eQTL beta
     when alleles are swapped), drop ambiguous/un-matchable SNPs;
  4. restrict to a +/-150 kb window around the cis-eQTL lead variant (PRIMARY);
  5. run coloc.abf; also report a +/-50..500 kb sensitivity scan.

Window note: PP.H4 for the ANKH->CAVS signal is robust (0.37-0.38) across all
windows; the +/-150 kb primary window reproduces the published H4:H3 ~ 13:1.
"""
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from coloc_abf import coloc_abf

DATA = os.path.join(os.path.dirname(__file__), "data")
PRIMARY_HALF = 150_000
SCAN = [50_000, 100_000, 150_000, 200_000, 300_000, 500_000]

GENE_ID = {"ANKH": "ENSG00000154122",
           "ENPP1": "ENSG00000197594",
           "PHOSPHO1": "ENSG00000173868"}

# ---------- loaders (return columns: pos, ea, oa, beta, se, [pvalue]) ----------
def load_eqtl_catalogue(fname, gene):
    d = pd.read_csv(os.path.join(DATA, fname), sep="\t")
    d = d[d.gene_id.str.startswith(GENE_ID[gene])].copy()      # MANDATORY filter
    d = d.rename(columns={"position": "pos", "ref": "oa", "alt": "ea"})
    return d[["pos", "ea", "oa", "beta", "se", "pvalue"]].dropna()

def load_gtex_aorta(fname, gene):  # figure-package GTEx file (variant_id encoded)
    d = pd.read_csv(os.path.join(DATA, fname), sep="\t")
    d = d[d.gene_id.str.startswith(GENE_ID[gene])].copy()
    parts = d.variant_id.str.split("_")
    d["pos"] = parts.str[1].astype(int); d["ea"] = parts.str[2]; d["oa"] = parts.str[3]
    d = d.rename(columns={"slope": "beta", "slope_se": "se", "pval_nominal": "pvalue"})
    return d[["pos", "ea", "oa", "beta", "se", "pvalue"]].dropna()

def load_outcome_finngen(fname):   # ref/alt with alt = effect allele
    d = pd.read_csv(os.path.join(DATA, fname)); d.columns = [c.lstrip("#") for c in d.columns]
    d = d.rename(columns={"alt": "ea", "ref": "oa", "beta": "b", "sebeta": "se"})
    return d[["pos", "ea", "oa", "b", "se"]].dropna()

def load_outcome_gwascat(fname):   # effect_allele/other_allele
    d = pd.read_csv(os.path.join(DATA, fname))
    d = d.rename(columns={"base_pair_location": "pos", "effect_allele": "ea",
                          "other_allele": "oa", "beta": "b", "standard_error": "se"})
    return d[["pos", "ea", "oa", "b", "se"]].dropna()

# ---------- harmonise + coloc on one window ----------
def coloc_window(eq, out, otype, half):
    c = int(eq.loc[eq.pvalue.idxmin(), "pos"])
    e = eq[(eq.pos >= c - half) & (eq.pos <= c + half)]
    o = out[(out.pos >= c - half) & (out.pos <= c + half)]
    m = e.merge(o, on="pos", suffixes=("_e", "_o"))
    same = (m.ea_e == m.ea_o) & (m.oa_e == m.oa_o)
    flip = (m.ea_e == m.oa_o) & (m.oa_e == m.ea_o)
    m = m[same | flip].copy(); m.loc[flip, "beta"] = -m.loc[flip, "beta"]
    if len(m) < 5:
        return None, len(m)
    res, n = coloc_abf(m.beta.values, m.se_e.values, m.b.values, m.se_o.values, "quant", otype)
    return res, n

# ---------- candidates ----------
CANDIDATES = [
    ("ANKH",     "blood", "eqtl_ANKH_blood.tsv",     load_eqtl_catalogue, "reg_CAVS_ANKH.csv",     load_outcome_finngen, "cc",   "ANKH blood -> CAVS"),
    ("ANKH",     "aorta", "eqtl_ANKH_aorta.tsv",     load_eqtl_catalogue, "reg_CAVS_ANKH.csv",     load_outcome_finngen, "cc",   "ANKH aorta -> CAVS"),
    ("ENPP1",    "aorta", "axis_GTExAorta_b38.tsv",       load_gtex_aorta,     "reg_PAD_ENPP1.csv",     load_outcome_finngen, "cc",   "ENPP1 aorta -> PAD"),
    ("PHOSPHO1", "blood", "eqtl_PHOSPHO1_blood.tsv", load_eqtl_catalogue, "reg_CAD_PHOSPHO1.csv",  load_outcome_gwascat, "cc",   "PHOSPHO1 blood -> CAD"),
    ("PHOSPHO1", "aorta", "eqtl_PHOSPHO1_aorta.tsv", load_eqtl_catalogue, "reg_CAD_PHOSPHO1.csv",  load_outcome_gwascat, "cc",   "PHOSPHO1 aorta -> CAD"),
]

def main():
    rows = []
    print(f"{'pair':24s} {'nSNP':>5s}  {'PP.H4':>6s} {'PP.H3':>6s} {'H4:H3':>6s}   note")
    for gene, tissue, ef, eloader, of, oloader, otype, label in CANDIDATES:
        try:
            eq = eloader(ef, gene)
        except Exception as e:
            print(f"{label:24s}  eQTL load error: {e}"); continue
        if len(eq) == 0:
            print(f"{label:24s}     -        -      -      -    {gene} not an eGene in this dataset (uninstrumentable)")
            rows.append(dict(pair=label, nSNP=0, PP_H4="", PP_H3="", note=f"{gene} not an eGene (uninstrumentable)"))
            continue
        out = oloader(of)
        res, n = coloc_window(eq, out, otype, PRIMARY_HALF)
        if res is None:
            print(f"{label:24s} {n:5d}      -      -      -    <5 overlapping SNPs")
            rows.append(dict(pair=label, nSNP=n, PP_H4="", PP_H3="", note="<5 overlapping SNPs"))
            continue
        ratio = res["H4"]/res["H3"] if res["H3"] > 0 else float("inf")
        note = "COLOCALISES" if res["H4"] > 0.5 else ("suggestive" if res["H4"] > 0.3 else "no coloc")
        print(f"{label:24s} {n:5d}  {res['H4']:.3f}  {res['H3']:.3f}  {ratio:5.1f}   {note}")
        rows.append(dict(pair=label, nSNP=n, PP_H0=round(res['H0'],4), PP_H1=round(res['H1'],4),
                         PP_H2=round(res['H2'],4), PP_H3=round(res['H3'],4), PP_H4=round(res['H4'],4),
                         H4_H3=round(ratio,2), window_kb=PRIMARY_HALF//1000, note=note))
    pd.DataFrame(rows).to_csv(os.path.join(os.path.dirname(__file__), "results_coloc.csv"), index=False)

    # sensitivity scan for the headline signal
    print("\nWindow sensitivity — ANKH blood -> CAVS (PP.H4 should be stable ~0.37):")
    eq = load_eqtl_catalogue("eqtl_ANKH_blood.tsv", "ANKH")
    out = load_outcome_finngen("reg_CAVS_ANKH.csv")
    srows = []
    print(f"  {'+/- kb':>7s} {'nSNP':>5s} {'PP.H4':>6s} {'PP.H3':>6s} {'H4:H3':>6s}")
    for half in SCAN:
        res, n = coloc_window(eq, out, "cc", half)
        if res:
            r = res["H4"]/res["H3"] if res["H3"]>0 else float("inf")
            print(f"  {half//1000:7d} {n:5d} {res['H4']:.3f}  {res['H3']:.3f} {r:6.1f}")
            srows.append(dict(window_kb=half//1000, nSNP=n, PP_H4=round(res['H4'],4),
                              PP_H3=round(res['H3'],4), H4_H3=round(r,2)))
    pd.DataFrame(srows).to_csv(os.path.join(os.path.dirname(__file__), "results_sensitivity_ANKH_CAVS.csv"), index=False)
    print("\n-> results_coloc.csv, results_sensitivity_ANKH_CAVS.csv written")

if __name__ == "__main__":
    main()
