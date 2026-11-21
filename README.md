# pTFCE: Probabilistic Threshold-Free Cluster Enhancement (HallquistLab Fork)

High-performance fork of [`spisakt/pTFCE`](https://github.com/spisakt/pTFCE) implementing **Probabilistic Threshold-Free Cluster Enhancement** (Spisák et al., 2019) with a fast C++ (`Rcpp`) engine, exact log-space posterior quadrature, and native bi-sided (two-tailed) fMRI inference.

> **Detailed Technical Report on Fork Modifications:** See [`docs/fork_modifications_report.qmd`](docs/fork_modifications_report.qmd) for mathematical derivations, benchmarks, and a full record of all statistical, numerical, and API departures from upstream `spisakt/pTFCE` (`v0.2.2.1`).

## Installation

```r
remotes::install_github("HallquistLab/pTFCE")
```

## Key Enhancements in This Fork (`v0.4.1`)

- **Fast `Rcpp` Engine (~400–900x Faster, Minimal Memory):**
  - Replaces repeated 3D connected-component labeling (`mmand::components` + `table`) and 4D intermediate voxel-by-threshold arrays with a **single-pass descending Union-Find (Disjoint-Set)** sweep in C++ (~0.33 s for a whole-brain bi-sided contrast vs. ~2.5–5 min previously).
  - Replaces slow triple-nested R loops in `smoothest()` with a cache-friendly, zero-copy C++ estimator for both 3D Z-maps (`0.006 s`) and 4D GLM residual images (`res4d`, `0.08 s`).
  - Drops `mmand` dependency; supports `RNifti` (`niftiImage`), `oro.nifti` (`nifti`, optional), plain R `array`s, and NIfTI file paths directly.
- **Native Bi-Sided (`sided = "two.sided"`) Enhancement:**
  - Enhances positive (`Z > 0`) and negative (`Z < 0`) excursion sets symmetrically using an identical threshold resolution, returning signed enhanced Z-scores (`$Z`), two-sided uncorrected p-values (`$p`), `-log(p)` (`$logp`), and two-sided FWER thresholds (`$fwer0.05.Z`, `$fwe.Z`).
- **Optional Non-Stationary Smoothness Adjustment (`nonstationary = TRUE`):**
  - Implements Worsley et al. (1999) / Hayasaka et al. (2004) / Salimi-Khorshidi et al. (2011) non-stationary cluster extent adjustment using spatially regularized voxelwise roughness weights `w(v) = dLh(v) / mean(dLh)` from 4D residuals (`res4d`), 4D `copes`, or precomputed `local_smoothness` / `RPV` maps—equalizing regional false-positive rates across smooth cortex vs. rougher subcortex. Defaults to `nonstationary = FALSE` for backward compatibility.
- **Empirical FWER Threshold Calibration (`fwe_method = c("grf", "grf_sim", "sign_flip")`):**
  - Because pTFCE smooths supra-threshold peaks, classical raw-field Euler characteristic thresholds (`fwe_method = "grf"`, default) are conservative (`FWER ~ 0.005–0.015` at nominal `alpha = 0.05`).
  - **`"sign_flip"`**: Fast Rademacher (`+/- 1`) sign-flipping in C++ across 4D subject contrasts (`copes`) or 4D GLM residuals (`residual`), achieving exact nominal FWER (`~0.050`) in ~1–7 seconds whole-brain while respecting empirical spatial covariance, non-stationarity, and heavy-tailed noise.
  - **`"grf_sim"`**: Parametric 2D/3D GRF Monte Carlo simulation (stationary anisotropic or multi-scale non-stationary) requiring no 4D subject/residual stacks.
- **Grid-Invariant Integration (`dh = 0.05` + Analytical Sub-`ZestThr`):**
  - Evaluates the unenhanced interval `[0, -log(Phi(-ZestThr))]` analytically and integrates above `ZestThr` with a fixed step `dh = 0.05`.
  - Eliminates cross-cluster interference from distant extreme Z-peaks (`Z_max` no longer coarsens step size `dh`) and removes the ~20% zero dead-zone around `Z = 0` as well as the `-8.21` out-of-mask floor. Legacy `Nh`-grid behavior remains available via `Nh = 100` or `legacy = TRUE`.
- **Exact Log-Space Bayesian Quadrature (Zero `NaN`s or Underflow):**
  - Replaces nested `stats::integrate(..., -Inf, Inf)` calls with mode-shifted log-space adaptive Gauss-Kronrod quadrature and per-cluster-size caching, completely eliminating `p > 1` cluster-border `NaN`s, silent score deflation inside large clusters, and `Z = 37.5` underflow ceilings.
- **Additional Utilities:**
  - `t2z(img, dof)` for exact log-space Student's *t*-to-*Z* conversion without overflow at high `|t|`.
  - `connectivity = c(6, 18, 26)` (matching AFNI `NN = 1, 2, 3`).
  - Vectorized log-space Newton-Raphson `fwe.p2z()` and `fwe.z2p()`.

## Quick Start

```r
library(pTFCE)

# Bi-sided non-stationary pTFCE with fast Rademacher sign-flipping FWER calibration:
res <- ptfce(
  img = "zstat1.nii.gz",
  mask = "mask.nii.gz",
  residual = "res4d.nii.gz",
  dof = 42,
  sided = "two.sided",
  nonstationary = TRUE,     # FALSE (default) for stationary cluster sizes
  fwe_method = "sign_flip", # "grf" (default), "sign_flip", or "grf_sim"
  n_perm = 500L,
  fwe_p = c(0.05, 0.01)
)

# Signed enhanced Z-map (RNifti niftiImage when given file paths)
RNifti::writeNifti(res$Z, "zstat1_ptfce.nii.gz")

# Threshold at calibrated two-sided FWER p < 0.05
z_thr <- res$Z
z_thr[abs(z_thr) < res$fwer0.05.Z] <- 0
RNifti::writeNifti(z_thr, "zstat1_ptfce_fwep0.05.nii.gz")
```

## Reference

> Tamás Spisák, Zsófia Spisák, Matthias Zunhammer, Ulrike Bingel, Stephen Smith, Thomas Nichols, Tamás Kincses (2019). *Probabilistic TFCE: a generalized combination of cluster size and voxel intensity to increase statistical power.* NeuroImage, 185, 12–26. [DOI: 10.1016/j.neuroimage.2018.09.078](https://doi.org/10.1016/j.neuroimage.2018.09.078)
