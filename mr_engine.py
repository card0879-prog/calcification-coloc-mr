"""
mr_engine.py — two-sample Mendelian randomisation estimators
================================================================================
Re-implements the TwoSampleMR estimators used in this project, on harmonised
SNP-level data (columns: beta.exposure, se.exposure, beta.outcome, se.outcome).

Estimators: per-SNP Wald ratio, inverse-variance weighted (IVW), MR-Egger,
weighted median; plus single-SNP and leave-one-out diagnostics.
Conventions follow TwoSampleMR (Hemani et al. 2018) so results match its output.
"""
import numpy as np, pandas as pd
from scipy import stats

def wald_ratio(bx, by, sx, sy):
    b = by / bx
    se = sy / np.abs(bx)                       # first-order (TwoSampleMR mr_wald_ratio)
    return b, se

def _wls_through_origin(bx, by, w):
    """weighted least squares  by ~ -1 + bx ; returns (slope, se, sigma)."""
    b = np.sum(w * bx * by) / np.sum(w * bx**2)
    resid = by - b * bx
    n = len(bx)
    sigma = np.sqrt(np.sum(w * resid**2) / (n - 1))         # residual SE
    se_raw = np.sqrt(1.0 / np.sum(w * bx**2))
    return b, se_raw, sigma

def mr_ivw(bx, by, sx, sy):
    w = 1.0 / sy**2
    b, se_raw, sigma = _wls_through_origin(bx, by, w)
    # TwoSampleMR: se = lm Std.Error / min(1, sigma) = (sigma * se_raw)/min(1, sigma)
    #            = se_raw * max(1, sigma)  -> keeps over-dispersion, floors under-dispersion
    se = se_raw * max(1.0, sigma)
    p = 2 * stats.norm.sf(abs(b / se))
    return dict(method="IVW", b=b, se=se, pval=p, nsnp=len(bx))

def mr_egger(bx, by, sx, sy):
    # orient exposure effects positive (TwoSampleMR convention)
    sign = np.sign(bx); sign[sign == 0] = 1
    bx_o, by_o = bx * sign, by * sign
    w = 1.0 / sy**2
    X = np.column_stack([np.ones_like(bx_o), bx_o])
    W = np.diag(w)
    XtW = X.T @ W
    beta = np.linalg.solve(XtW @ X, XtW @ by_o)
    resid = by_o - X @ beta
    n = len(bx_o)
    sigma = np.sqrt((resid @ (w * resid)) / (n - 2))
    cov = np.linalg.inv(XtW @ X)
    se_slope = np.sqrt(cov[1, 1]) * max(sigma, 1e-12)
    se_int   = np.sqrt(cov[0, 0]) * max(sigma, 1e-12)
    b = beta[1]; se = se_slope
    p = 2 * stats.t.sf(abs(b / se), df=n - 2)
    pint = 2 * stats.t.sf(abs(beta[0] / se_int), df=n - 2)
    return dict(method="MR-Egger", b=b, se=se, pval=p,
                egger_intercept=beta[0], intercept_pval=pint, nsnp=n)

def mr_weighted_median(bx, by, sx, sy, n_boot=1000, seed=1):
    b_iv = by / bx
    w = (bx**2) / (sy**2)                       # weights ~ 1/se(wald)^2
    order = np.argsort(b_iv); b_s, w_s = b_iv[order], w[order]
    cw = np.cumsum(w_s) - 0.5 * w_s; cw /= np.sum(w_s)
    below = np.where(cw < 0.5)[0][-1]
    b = b_s[below] + (b_s[below+1]-b_s[below]) * (0.5-cw[below])/(cw[below+1]-cw[below])
    rng = np.random.default_rng(seed); boots = np.empty(n_boot)
    for i in range(n_boot):
        by_b = rng.normal(by, sy); bx_b = rng.normal(bx, sx)
        biv = by_b / bx_b
        o = np.argsort(biv); bs, ws = biv[o], ((bx_b**2)/(sy**2))[o]
        c = np.cumsum(ws) - 0.5*ws; c /= np.sum(ws)
        k = np.where(c < 0.5)[0][-1]
        boots[i] = bs[k] + (bs[k+1]-bs[k])*(0.5-c[k])/(c[k+1]-c[k])
    se = np.std(boots); p = 2*stats.norm.sf(abs(b/se))
    return dict(method="Weighted median", b=b, se=se, pval=p, nsnp=len(bx))

def single_snp(bx, by, sx, sy, snp):
    b, se = wald_ratio(bx, by, sx, sy)
    p = 2*stats.norm.sf(np.abs(b/se))
    return pd.DataFrame(dict(SNP=snp, b=b, se=se, p=p))

def leave_one_out(bx, by, sx, sy, snp):
    rows = []
    for i in range(len(bx)):
        m = np.arange(len(bx)) != i
        r = mr_ivw(bx[m], by[m], sx[m], sy[m])
        rows.append(dict(SNP=snp[i], b=r["b"], se=r["se"], p=r["pval"]))
    return pd.DataFrame(rows)
