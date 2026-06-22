"""
run_mr.py — two-sample MR on harmonised SNP-level data
================================================================================
Input : data/scatter_dat.csv  (SNP, beta.exposure, se.exposure, beta.outcome, se.outcome)
Output: results_mr.csv            (IVW / MR-Egger / Weighted median summary)
        results_mr_singlesnp.csv  (per-SNP Wald ratios)
        results_mr_loo.csv        (leave-one-out IVW)
"""
import os, numpy as np, pandas as pd
from mr_engine import mr_ivw, mr_egger, mr_weighted_median, single_snp, leave_one_out

HERE = os.path.dirname(__file__)
d = pd.read_csv(os.path.join(HERE, "data", "scatter_dat.csv"))
bx, sx = d["beta.exposure"].values, d["se.exposure"].values
by, sy = d["beta.outcome"].values,  d["se.outcome"].values
snp = d["SNP"].values

summary = pd.DataFrame([mr_ivw(bx,by,sx,sy), mr_egger(bx,by,sx,sy),
                        mr_weighted_median(bx,by,sx,sy)])
summary.to_csv(os.path.join(HERE,"results_mr.csv"), index=False)
single_snp(bx,by,sx,sy,snp).to_csv(os.path.join(HERE,"results_mr_singlesnp.csv"), index=False)
leave_one_out(bx,by,sx,sy,snp).to_csv(os.path.join(HERE,"results_mr_loo.csv"), index=False)

print(summary[["method","b","se","pval","nsnp"]].to_string(index=False))
print("\n-> results_mr.csv, results_mr_singlesnp.csv, results_mr_loo.csv written")
