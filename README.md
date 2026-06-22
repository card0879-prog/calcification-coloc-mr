# Calcification drug-target coloc + MR engine

Reproducible colocalisation and Mendelian-randomisation analysis for the
PPi-TNAP calcification axis (ALPL / ENPP1 / ANKH / PHOSPHO1). Self-contained:
code + input data + result tables. Running the scripts reproduces every number
in the paper.

**Quick start**
```
pip install numpy pandas scipy matplotlib
python run_coloc.py      # colocalisation + sensitivity scan
python run_mr.py          # Mendelian randomisation estimates
python validate_mr.py     # proves MR output matches TwoSampleMR
```

**Repository layout**
```
coloc_abf.py  run_coloc.py      colocalisation engine + pipeline
mr_engine.py  run_mr.py  validate_mr.py   MR estimators + validation
data/                             input region files (see DATA.md)
results_*.csv                     generated outputs
```

**How to cite** — see `CITATION.cff` (GitHub shows a "Cite this repository" button).

**Licences** — code: MIT (`LICENSE`); data: CC-BY 4.0 with source attribution (`DATA.md`).

---

# Colocalisation engine — PPi-TNAP calcification drug-target axis

Bayesian colocalisation (coloc.abf, Wakefield ABF, single-causal-variant) for the
ALPL / ENPP1 / ANKH / PHOSPHO1 axis against vascular-calcification and
cardiovascular outcomes. This package is self-contained: code + the exact input
data + generated result tables. Running `run_coloc.py` reproduces the numbers below.

## Files
- `coloc_abf.py` — the ABF colocalisation engine (Wakefield 2009; Giambartolomei 2014).
- `run_coloc.py` — pipeline: gene filtering, allele harmonisation, windowing, coloc.
- `data/` — input region files (eQTL-Catalogue cis regions; FinnGen / GWAS-Catalogue outcomes).
- `results_coloc.csv`, `results_sensitivity_ANKH_CAVS.csv` — generated outputs.

## Run
```
pip install numpy pandas scipy
python run_coloc.py
```

## Method (fixed parameters)
- Engine: Wakefield approximate Bayes factors, single-causal-variant assumption.
- Priors: p1 = 1e-4, p2 = 1e-4, p12 = 1e-5.
- Prior effect SD: quantitative W = 0.15^2; case-control W = 0.20^2.
- **Gene filtering (mandatory):** eQTL-Catalogue region files contain every gene in
  the locus. Each exposure is filtered to its target `gene_id`
  (ANKH ENSG00000154122, ENPP1 ENSG00000197594, PHOSPHO1 ENSG00000173868)
  before colocalisation. Omitting this filter dilutes the signal and inflates PP.H3.
- SNP matching: by chr:position; alleles harmonised to a common effect allele,
  eQTL beta sign-flipped on allele swap; non-matching alleles dropped.
- Window: +/- 150 kb around the cis-eQTL lead variant (primary). A +/- 50..500 kb
  sensitivity scan is reported for the ANKH->CAVS signal.

## Results (this package, reproduced)
| Exposure (tissue) -> Outcome | nSNP | PP.H4 | PP.H3 | H4:H3 | verdict |
|---|---|---|---|---|---|
| ANKH (blood) -> CAVS  | 976 | 0.372 | 0.031 | 12.0 | sole coloc-supported signal |
| ANKH (aorta) -> CAVS  | 778 | 0.012 | 0.010 |  1.2 | uninstrumentable in aorta |
| ENPP1 (aorta) -> PAD  |  80 | 0.008 | 0.001 | 10.5 | no colocalisation |
| PHOSPHO1 (blood) -> CAD | 866 | 0.027 | 0.172 | 0.2 | no colocalisation (H3-dominant) |
| PHOSPHO1 (aorta) -> CAD | 843 | 0.065 | 0.077 | 0.8 | no colocalisation (H3-dominant) |

ANKH->CAVS PP.H4 is stable at 0.366-0.379 across +/- 50-500 kb windows; the
+/- 150 kb window gives H4:H3 = 12.0. This is the project's single positive
colocalisation; the remaining axis genes are eliminated or uninstrumentable.

## Note on PHOSPHO1
PHOSPHO1 (Ensembl gene ENSG00000173868) is a tested cis-eGene in both whole blood
and aorta. Colocalisation with CAD is not supported (PP.H4 = 0.027 blood, 0.065
aorta; PP.H3 dominant), i.e. the eQTL and CAD association peaks are driven by
distinct variants. This eliminates PHOSPHO1 as a colocalising target.


## Data provenance
- eQTL: eQTL-Catalogue cis region summary statistics (whole blood; GTEx aorta).
- Outcomes: FinnGen (CAVS, PAD; b38, alt = effect allele); GWAS-Catalogue harmonised
  CARDIoGRAMplusC4D (CAD).

---

# Mendelian randomisation engine

Two-sample MR estimators re-implemented to match TwoSampleMR (Hemani et al. 2018).

## Files
- `mr_engine.py` — Wald ratio, IVW, MR-Egger, weighted median, single-SNP, leave-one-out.
- `run_mr.py` — runs all estimators on `data/scatter_dat.csv`, writes `results_mr*.csv`.
- `validate_mr.py` — reproduces the archived TwoSampleMR output to machine precision.

## Run
```
python run_mr.py        # produce MR estimates + single-SNP + leave-one-out
python validate_mr.py   # prove the engine matches TwoSampleMR (deltas ~1e-16)
```

## Method
- Per-SNP Wald ratio: b = beta_out/beta_exp, se = se_out/|beta_exp|.
- IVW: weighted regression of outcome on exposure through the origin, weights 1/se_out^2;
  standard error keeps over-dispersion (se = se_raw * max(1, residual sigma)),
  matching TwoSampleMR's `mr_ivw`.
- MR-Egger: weighted regression with intercept (directional-pleiotropy term).
- Weighted median: 1000-iteration parametric bootstrap SE.

## Validation (machine precision vs TwoSampleMR)
| estimate | Δb | Δse |
|---|---|---|
| IVW   | 7e-18  | 2e-17 |
| MR-Egger | 1e-16 | 0 |
| per-SNP Wald (n=492) | 5e-15 | 5e-15 |
| leave-one-out (n=492) | 1e-16 | 7e-17 |

Reference IVW = 0.038617 (se 0.026171); MR-Egger = 0.015991 (se 0.045464),
matching `data/singlesnp.csv` exactly.

## Dependencies
Python 3, numpy, pandas, scipy.
