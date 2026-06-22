"""
coloc_abf.py — Bayesian colocalisation via Wakefield Approximate Bayes Factors
================================================================================
Single-causal-variant assumption (Giambartolomei et al. 2014, PLoS Genet),
ABF per Wakefield (2009, Genet Epidemiol).

Posterior over 5 hypotheses for a genomic region:
  H0: no causal variant for either trait
  H1: causal variant for trait 1 only (e.g. eQTL exposure)
  H2: causal variant for trait 2 only (e.g. disease outcome)
  H3: distinct causal variants for the two traits
  H4: a single shared causal variant (colocalisation)

Priors (coloc defaults, fixed for this project):
  p1 = 1e-4   p2 = 1e-4   p12 = 1e-5
Prior effect-size SD (sd.prior):
  quantitative trait  W = 0.15**2
  case-control trait  W = 0.20**2
"""
import numpy as np
from scipy.special import logsumexp

P1, P2, P12 = 1e-4, 1e-4, 1e-5
W_QUANT, W_CC = 0.15**2, 0.20**2

def _lABF(beta, se, W):
    """log Wakefield ABF per SNP for a single trait."""
    beta = np.asarray(beta, float); se = np.asarray(se, float)
    V = se**2
    z2 = (beta / se)**2
    r = W / (W + V)
    return 0.5 * np.log1p(-r) + 0.5 * z2 * r

def coloc_abf(beta1, se1, beta2, se2, type1, type2,
              p1=P1, p2=P2, p12=P12):
    """
    Colocalise two aligned summary-stat vectors (same SNPs, same effect allele).
    type1/type2 in {"quant","cc"}. Returns (dict PP.H0..H4, n_snps).
    """
    W1 = W_QUANT if type1 == "quant" else W_CC
    W2 = W_QUANT if type2 == "quant" else W_CC
    l1, l2 = _lABF(beta1, se1, W1), _lABF(beta2, se2, W2)
    sabf1, sabf2, sabf12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    # H3: product of marginals minus the shared (H4) term, in log space
    logdiff = lambda a, b: a + np.log1p(-np.exp(b - a))
    lH = np.array([
        0.0,                                          # H0
        np.log(p1) + sabf1,                           # H1
        np.log(p2) + sabf2,                           # H2
        np.log(p1) + np.log(p2) + logdiff(sabf1 + sabf2, sabf12),  # H3
        np.log(p12) + sabf12,                         # H4
    ])
    pp = np.exp(lH - logsumexp(lH))
    return dict(zip(["H0", "H1", "H2", "H3", "H4"], pp)), int(np.size(l1))
