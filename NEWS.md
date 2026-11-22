# pTFCE 0.4.2

## Bug Fixes & 2D / Non-Stationary Smoothness Enhancements

* Fixed external 2D/3D `local_smoothness` handling (e.g., SPM `RPV.nii` or `dLh` maps passed without manual `V`/`Rd`/`resels`, as in Vignette Recipe C): `ptfce()` now derives global Riemannian `V`, `Rd`, and `resels` directly from `RPV`/`dLh` maps via `local_smoothness_units = c("auto", "rpv", "dlh", "weights")` instead of falling back to 3D `Z`-map `smoothest()`.
* Added exact 2D (`D = 2`) Gaussian Random Field cluster-extent and Expected Euler Characteristic formulas across `Es()`, `ptfce_tail_cpp()`, `NullPtfceEngine`, `fwe.z2p()`, and `fwe.p2z()`.
* Added small-ROI rank-cap and `dof`-mismatch diagnostics in `calibrate_fwer_signflip_cpp()` when `n_act < nt` or supplied `dof` diverges from empirical residual rank.
* Winsorized extreme supra-threshold excursions (`|Z| > 3`) during single-volume 3D (`nt == 1`) local roughness estimation in `smoothest_cpp()` so strong focal activations do not suppress their own local roughness weights, and made `interpolate(v < 5)` monotone.

# pTFCE 0.4.1

## Empirical FWER Threshold Calibration

* Added `fwe_method = c("grf", "grf_sim", "sign_flip")`, `n_perm = 500L`, and `copes = NULL` to `ptfce()` (default `"grf"` preserves exact backward compatibility) to overcome the conservativeness of raw-field Euler characteristic FWER thresholds on pTFCE-enhanced maps (`FWER ~ 0.010–0.015` at nominal `alpha = 0.05`).
* Implemented `calibrate_fwer_signflip_cpp()`: fast SIMD-friendly Rademacher (`+/- 1`) sign-flipping calibration across either 4D subject contrast stacks (`copes`, auto-estimating 4D residual smoothness when `residual` is omitted) or 4D GLM residuals (`residual`, using two-pass Modified Gram-Schmidt projection onto the empirical residual subspace).
* Implemented `calibrate_fwer_grf_cpp()`: parametric 2D/3D Gaussian Random Field Monte Carlo simulator supporting both stationary anisotropic smoothness and multi-scale non-stationary smoothness (`nonstationary = TRUE`) with exact analytical cross-scale variance normalization.
* Added lazy 2D posterior lookup-table (`NullPtfceEngine`) and supra-`ZestThr` candidate pre-filtering in C++ to evaluate `max(Z)` / `max(|Z|)` across hundreds of null replicates with zero per-replicate quadrature or heap churn.
* Added practitioner vignette (`vignette("ptfce-fmri-guide", package = "pTFCE")`) and reproducible Monte Carlo simulation suite (`docs/benchmarks/run_simulation_studies.R`).

# pTFCE 0.4.0

## Non-Stationary Smoothness Adjustment

* Added optional non-stationary smoothness adjustment (`nonstationary = TRUE` in `ptfce()` and `smoothest()`, default `FALSE`) based on Worsley et al. (1999), Hayasaka et al. (2004), and Salimi-Khorshidi et al. (2011), collapsing smooth-vs-rough regional false-positive disparity from `17.8x` to `1.1x` in 3D group-fMRI simulations.
* Extended `smoothest()` C++ engine with normalized within-mask separable 3D Gaussian smoothing (`smooth_fwhm = 3.0`) and empirical-Bayes shrinkage of voxelwise derivative sufficient statistics, outputting `dLh_map`, `rpv_map` (RESELs per voxel), `dLh_nonstationary` (Riemannian integral `mean(dLh_map)`), `resels_nonstationary`, and mean-1 `voxel_weights` (winsorized to `[0.02, 50.0]`).
* Extended `ptfce_tail_cpp()` with continuous resel-weighted cluster mass accumulation and fast 4-point cubic Lagrange log-space interpolation tables.
* Added `local_smoothness` parameter to `ptfce()` allowing reuse of precomputed `smoothest(..., nonstationary = TRUE)` objects or external RPV/roughness maps (e.g., SPM `RPV.nii`) across multiple contrasts.

# pTFCE 0.3.0

## HallquistLab Fork Release (C++ Engine & Native Bi-Sided Inference)

* Replaced connected-components loop (`mmand::components` + `table`) and nested `stats::integrate` calls with a C++ (`Rcpp`) single-pass descending Union-Find engine and log-space adaptive Gauss-Kronrod posterior quadrature (~400–900x faster, ~50x lower RAM, drops `mmand` dependency).
* Eliminated numerical underflow and quadrature spike-miss bugs (`pvox.clust > 1`, cluster-border `NaN`s, silent score deflation inside large clusters, and artificial `745` / `Z = 37.53` ceiling).
* Added native bi-sided enhancement (`sided = c("positive", "two.sided", "negative")`) with shared symmetric threshold grids, signed enhanced Z-scores, two-sided p-values, and two-sided FWER thresholds.
* Switched default threshold integration (`dh = 0.05`, `legacy = FALSE`) to exact analytical sub-`ZestThr` integration plus fixed-`dh` midpoint integration above `ZestThr`, eliminating cross-cluster interference from distant extreme Z-peaks and the ~20% zero dead-zone around `Z = 0` (legacy `Nh`-grid mode preserved via `Nh = ...` or `legacy = TRUE`).
* Fixed `if (!is.na(residual))` error on R >= 4.2, 3D trailing-1 dimension crash in `smoothest()`, 2D `zstart` loop bug, DOF lookup interpolation bug in `interpolate()`, out-of-mask `-8.2095` floor, and spurious warnings from out-of-mask `NA`/`Inf` values.
* Rewrote `smoothest()` backend in C++ for fast, zero-copy 2D/3D/4D residual smoothness estimation.
* Added log-space `t2z()` helper, `connectivity = c(6, 18, 26)` support, vectorized log-space Newton-Raphson `fwe.p2z()` / `fwe.z2p()`, and seamless support for `RNifti`, `oro.nifti` (optional in `Suggests`), plain R `array`s, and NIfTI file paths.

# pTFCE 0.2.2.1

* Upstream `spisakt/pTFCE` baseline (`30.04.2020`): bugfix Issue #8 tolerance adjustments and prior `0.0.3`–`0.2.2` releases.
