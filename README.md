# pTFCE: Probabilistic Threshold-Free Cluster Enhancement (HallquistLab Fork)

High-performance fork of [`spisakt/pTFCE`](https://github.com/spisakt/pTFCE) implementing **Probabilistic Threshold-Free Cluster Enhancement** (Spisák et al., 2019) with a fast C++ (`Rcpp`) engine, exact log-space posterior quadrature, and native bi-sided (two-tailed) fMRI inference.

> **Documentation Website & Simulation Studies (`https://hallquistlab.github.io/pTFCE/`):**
> - **Practitioner's Guide & Simulation Walkthrough:** [Rendered HTML Vignette](https://hallquistlab.github.io/pTFCE/articles/ptfce-fmri-guide.html) ([source `.Rmd`](vignettes/ptfce-fmri-guide.Rmd)) — conceptual overview, visual examples, Monte Carlo simulation results, and copy-paste recipes for FSL, SPM, and AFNI.
> - **Technical Departure Report:** [Rendered Quarto Report](https://hallquistlab.github.io/pTFCE/fork_modifications_report.html) ([source `.qmd`](docs/fork_modifications_report.qmd)) — mathematical derivations, C++ architecture, and complete departure log from upstream `spisakt/pTFCE` (`v0.2.2.1`).
> - **Function Reference & Changelog:** [API Reference](https://hallquistlab.github.io/pTFCE/reference/index.html) | [Changelog](https://hallquistlab.github.io/pTFCE/news/index.html) | [Reproducible Benchmark Suite (`run_simulation_studies.R`)](docs/benchmarks/run_simulation_studies.R).

## Installation

```r
remotes::install_github("HallquistLab/pTFCE", build_vignettes = TRUE)

# Open the user-friendly fMRI guide & simulation study results in R:
vignette("ptfce-fmri-guide", package = "pTFCE")
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

## Which Settings Should I Use?

| Available fMRI Inputs | Recommended `ptfce()` Call | Why |
| :--- | :--- | :--- |
| **4D GLM residuals** (`res4d` + `dof`) | `sided = "two.sided", nonstationary = TRUE, fwe_method = "sign_flip"` | **Recommended default for group fMRI.** Equalizes regional FPR, calibrates exact FWER (`~0.050`) in ~1–7 s whole-brain, and handles heavy-tailed outliers. |
| **4D subject contrasts** (`copes`, $N \ge 6$) | `copes = copes4d, sided = "two.sided", nonstationary = TRUE` | Auto-estimates 4D smoothness (`dof = N - 1`) and runs exact one-sample Rademacher sign-flipping. |
| **3D `Z`/`t` map + `RPV.nii`** (no 4D data) | `local_smoothness = "RPV.nii", sided = "two.sided", fwe_method = "grf_sim"` | Non-stationary cluster scaling + exact multi-scale GRF FWER calibration without 4D stacks. |
| **Single 3D `Z`-map + mask only** | `sided = "two.sided", nonstationary = FALSE, fwe_method = "grf_sim"` | Recovers exact stationary pTFCE FWER power without noisy single-map local smoothness estimates. |

## Empirical Simulation Highlights

In realistic 3D group-fMRI Monte Carlo studies spanning rough subcortex (`FWHM = 3.0` vox), lateral cortex (`FWHM = 4.8` vox), and smooth medial DMN (`FWHM = 7.2` vox; $N = 25$–$30$ subjects, nominal two-sided $\alpha = 0.05$):

- **Regional false-positive parity (`nonstationary = TRUE`):** Stationary pTFCE produces a **17.8$\times$** false-positive disparity between smooth DMN cortex and rough subcortex (`1.89` vs. `0.11` per $10^3$ voxels at $p < 0.001$); `nonstationary = TRUE` collapses this disparity to **1.1$\times$** (`0.34` vs. `0.31` per $10^3$).
- **Exact FWER & power recovery (`sign_flip` & `grf_sim`):** Analytical raw-field GRF (`fwe_method = "grf"`) yields conservative empirical FWER of `0.010–0.015`. Both `sign_flip` (`FWER = 0.050`) and `grf_sim` (`FWER = 0.045`) restore exact nominal Type I error control, boosting subcortical cluster power from `18.0%` (`NS_GRF`) and `24.0%` (`Raw Z`) to **`32.7%`** (doubling voxel sensitivity from `3.7%` to **`7.4%`**) and tripling lateral cortical power (`8.0%` $\to$ **`24.7%`**) with **`0.0%`** false-positive leakage in null smooth DMN cortex.

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
