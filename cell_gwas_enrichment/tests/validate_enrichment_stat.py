"""
validate_enrichment_stat.py

Standalone (R-free) validation of the statistical method implemented in
R/enrichment_engine.R, using a Python re-implementation of the same
closed-form Frisch-Waugh-Lovell (FWL) regression. This is a numerical/
statistical sanity check of the *method*, run once during development
(no R installation was available in the environment this pipeline was
authored in) -- it does not exercise the R code itself, so re-run the R
unit tests (if you add any, e.g. with testthat) against real Signac
objects before trusting results on real data.

Checks:
  1. The FWL closed-form coefficient matches a naive multivariate OLS fit
     (np.linalg.lstsq) exactly, with and without covariates.
  2. Under the global null (no true association between annotation and
     effect magnitude), empirical permutation p-values are ~Uniform(0,1)
     across many simulated cells (Kolmogorov-Smirnov test).
  3. Under an injected true enrichment signal, the test has high power to
     detect it at a nominal alpha of 0.05.

Run: python3 validate_enrichment_stat.py
"""
import numpy as np
from scipy import stats as st

rng = np.random.default_rng(42)


def fwl_beta1_matrix(X, y, C):
    """Closed-form OLS coefficient(s) on X (n x k, k annotation columns),
    controlling for covariates C (n x p, include an intercept column),
    via the Frisch-Waugh-Lovell theorem. Matches the R implementation in
    enrichment_engine.R (which uses stats::lm.fit instead of an explicit
    projector, for memory efficiency at scale -- these are numerically
    equivalent)."""
    CtC_inv = np.linalg.pinv(C.T @ C)
    Hc = C @ CtC_inv @ C.T
    y_resid = y - Hc @ y
    X_resid = X - Hc @ X
    num = (X_resid * y_resid[:, None]).sum(axis=0)
    den = (X_resid ** 2).sum(axis=0)
    return num / den


def naive_ols_multi(x, y, C):
    design = np.column_stack([C, x])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return coef[-1]


def test_fwl_matches_naive_ols():
    n = 2000
    maf = rng.uniform(0.01, 0.5, n)
    ld = rng.exponential(1.0, n)
    C = np.column_stack([np.ones(n), maf, ld])
    chi2 = rng.chisquare(df=1, size=n)
    x_obs = rng.binomial(1, 0.1, n).astype(float)

    beta_fwl = fwl_beta1_matrix(x_obs[:, None], chi2, C)[0]
    beta_naive = naive_ols_multi(x_obs, chi2, C)
    assert np.isclose(beta_fwl, beta_naive), (beta_fwl, beta_naive)
    print(f"[PASS] FWL closed form matches naive OLS: {beta_fwl:.6f} vs {beta_naive:.6f}")


def test_null_calibration(n=2000, nperm=500, ncells_sim=300):
    maf = rng.uniform(0.01, 0.5, n)
    ld = rng.exponential(1.0, n)
    C = np.column_stack([np.ones(n), maf, ld])
    chi2 = rng.chisquare(df=1, size=n)  # no true signal anywhere

    pvals = []
    for _ in range(ncells_sim):
        n_frag_variants = rng.integers(20, 80)
        x_obs = np.zeros(n)
        x_obs[rng.choice(n, n_frag_variants, replace=False)] = 1
        X_perm = np.zeros((n, nperm))
        for k in range(nperm):
            X_perm[rng.choice(n, n_frag_variants, replace=False), k] = 1
        beta_obs = fwl_beta1_matrix(x_obs[:, None], chi2, C)[0]
        beta_perm = fwl_beta1_matrix(X_perm, chi2, C)
        pvals.append((1 + np.sum(beta_perm >= beta_obs)) / (nperm + 1))

    pvals = np.array(pvals)
    ks = st.kstest(pvals, "uniform")
    print(f"[INFO] null calibration: mean p = {pvals.mean():.3f} (expect ~0.5), "
          f"KS p-value = {ks.pvalue:.3f} (expect > 0.05, i.e. consistent with uniform)")
    assert ks.pvalue > 0.01, "Null p-values deviate from Uniform(0,1) -- method miscalibrated"
    print("[PASS] Null empirical p-values are consistent with Uniform(0,1)")


def test_power_under_true_signal(n=2000, nperm=500, ncells_sim=200):
    maf = rng.uniform(0.01, 0.5, n)
    ld = rng.exponential(1.0, n)
    C = np.column_stack([np.ones(n), maf, ld])
    chi2 = rng.chisquare(df=1, size=n)

    causal_idx = rng.choice(n, 100, replace=False)
    chi2_signal = chi2.copy()
    chi2_signal[causal_idx] += rng.chisquare(df=1, size=100) * 3

    noncausal_pool = np.setdiff1d(np.arange(n), causal_idx)
    pvals = []
    for _ in range(ncells_sim):
        n_frag_variants = 50
        n_causal_in_frag = min(rng.binomial(n_frag_variants, 0.6), len(causal_idx))
        chosen = np.concatenate([
            rng.choice(causal_idx, n_causal_in_frag, replace=False),
            rng.choice(noncausal_pool, n_frag_variants - n_causal_in_frag, replace=False),
        ])
        x_obs = np.zeros(n); x_obs[chosen] = 1
        X_perm = np.zeros((n, nperm))
        for k in range(nperm):
            X_perm[rng.choice(n, n_frag_variants, replace=False), k] = 1
        beta_obs = fwl_beta1_matrix(x_obs[:, None], chi2_signal, C)[0]
        beta_perm = fwl_beta1_matrix(X_perm, chi2_signal, C)
        pvals.append((1 + np.sum(beta_perm >= beta_obs)) / (nperm + 1))

    power = np.mean(np.array(pvals) < 0.05)
    print(f"[INFO] power at alpha=0.05 under injected true signal: {power:.2f}")
    assert power > 0.7, "Power unexpectedly low for a strong injected signal"
    print("[PASS] Method detects a true enrichment signal with high power")


if __name__ == "__main__":
    test_fwl_matches_naive_ols()
    test_null_calibration()
    test_power_under_true_signal()
    print("\nAll validation checks passed.")
