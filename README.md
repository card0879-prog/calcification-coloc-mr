# PPi–TNAP calcification axis: exploratory colocalisation and MR analyses

> **Important: preliminary-analysis archive, not the reproducibility package for the submitted manuscript.**
>
> This repository preserves an **earlier pilot analysis** of ALPL, ENPP1, ANKH and PHOSPHO1. It **does not reproduce** the final analyses, tables, or figures in the manuscript *“Colocalization-gated validity audit of drug-target Mendelian randomization in the PPi–TNAP axis”* (submitted to *BMC Cardiovascular Disorders*). In particular, the pilot uses different locus windows, outcomes, data resources and analysis versions. Results here **must not be cited as the results of the submitted manuscript**.
>
> The final manuscript's analysis code is **not currently contained in this repository**. The authors are auditing and organising the final-version code and provenance separately. No claim of end-to-end manuscript reproducibility is made for this repository.

## Scope and provenance

This repository is retained to document exploratory calculations and sensitivity checks performed during the development of a colocalisation-based drug-target MR audit. It includes archived inputs and outputs for the **pilot analysis only**. The pilot used an approximately **±150 kb** lead-variant-centred window for its primary colocalisation runs and includes **CAVS, PAD and CAD** outcomes; these settings should not be assumed to apply to the final paper.

## Repository contents

- `coloc_abf.py`, `run_coloc.py`: pilot colocalisation engine and pipeline
- `mr_engine.py`, `run_mr.py`, `validate_mr.py`: pilot MR implementations and validation
- `data/`, `DATA.md`: inputs and source documentation for the pilot
- `results_coloc.csv`, `results_mr*.csv`, `results_sensitivity_ANKH_CAVS.csv`: pilot outputs
- `LICENSE`, `CITATION.cff`: licensing and citation metadata

## Running the archived pilot

```bash
pip install numpy pandas scipy matplotlib
python run_coloc.py
python run_mr.py
python validate_mr.py
```

Running these scripts is intended to reproduce **the archived pilot results** where all listed inputs and software dependencies remain available. **It is not a replication of the final submitted study.**

## Pilot-specific methodological details

- Colocalisation engine: Wakefield approximate Bayes factors (`coloc.abf`), with a single-causal-variant approximation
- Priors: `p1 = 1e-4`, `p2 = 1e-4`, `p12 = 1e-5`
- Prior effect standard deviations: quantitative `0.15`, case-control `0.20`
- Pilot primary window: **±150 kb** around the lead cis-eQTL variant
- Pilot input resources included eQTL Catalogue cis-region summary statistics and FinnGen/GWAS Catalog outcomes

For example, one pilot run yielded `ANKH (blood) → CAVS: PP.H4 = 0.372`. **This is a historical pilot result, not the final manuscript estimate** (reported as `PP.H4 = 0.27`). Do not interchange these results.

## Reuse and citation

Code and data remain subject to the licences and source-attribution requirements in `LICENSE` and `DATA.md`. If you use this repository, explicitly identify it as a *preliminary exploratory code archive*. Do **not** describe it as the published/submitted manuscript's fully reproducible analysis package.
