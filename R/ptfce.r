#' @useDynLib pTFCE, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Probabilistic Threshold-Free Cluster Enhancement (pTFCE)
#'
#' Enhances 2D or 3D statistical Z-score neuroimages by integrating multi-threshold
#' cluster extent information via Bayes' rule and Gaussian Random Field (GRF) theory,
#' producing enhanced p-values and Z-scores directly without permutation testing.
#'
#' @references
#' T. Spisák, Z. Spisák, M. Zunhammer, U. Bingel, S. Smith, T. Nichols, T. Kincses (2019).
#' Probabilistic TFCE: a generalized combination of cluster size and voxel intensity
#' to increase statistical power. \emph{NeuroImage}, 185, 12-26.
#' \doi{10.1016/j.neuroimage.2018.09.078}
#'
#' @param img A 2D or 3D Z-score image to enhance: an \code{RNifti} \code{niftiImage},
#'   an \code{oro.nifti} \code{nifti} object, a numeric \code{array}, or a file path
#'   to a NIfTI file.
#' @param mask Spatial mask image of matching dimensions (optional; defaults to all
#'   voxels if omitted). Continuous masks are thresholded at \code{> 0.5}.
#' @param Rd Roughness-volume product \code{V * dLh} (FSL \code{smoothest}: \code{VOLUME * DLH}).
#' @param V Number of voxels in mask (FSL \code{smoothest}: \code{VOLUME}).
#' @param resels Resel size in voxels per resel (FSL \code{smoothest}: \code{RESELS}).
#' @param residual Optional 4D residual image (or file path) for accurate smoothness
#'   estimation via \code{\link{smoothest}}.
#' @param dof Residual degrees of freedom (required when \code{residual} is supplied).
#' @param logpmin Minimum \code{-log(p)} threshold used in legacy discrete grid mode (default: \code{0}).
#' @param logpmax Maximum \code{-log(p)} threshold (optional; if supplied, activates legacy grid mode).
#' @param ZestThr Cluster-forming Z threshold below which \code{P(H > h | C)} is treated
#'   as unenhanced \code{P(H > h)} due to GRF percolation limits (default: \code{1.3}).
#' @param Nh Number of discrete thresholds for legacy grid mode (default: \code{NULL},
#'   which uses the fixed-\code{dh} continuous sub-threshold integration grid; set e.g.
#'   \code{Nh = 100} or \code{legacy = TRUE} to reproduce upstream grid spacing).
#' @param verbose Logical; print progress/diagnostic messages if \code{TRUE} (default).
#' @param sided Direction of enhancement: \code{"positive"} (default, upper tail \code{Z > 0}),
#'   \code{"two.sided"} (enhances both positive and negative tails symmetrically using a
#'   shared threshold grid and returns signed enhanced Z-scores and two-sided p-values),
#'   or \code{"negative"} (lower tail \code{Z < 0}).
#' @param dh Step size in \code{-log(p)} space above \code{ZestThr} for the default
#'   grid-invariant integration mode (default: \code{0.05}). Ignored when \code{Nh} or
#'   \code{legacy = TRUE} is specified.
#' @param connectivity Voxel neighborhood connectivity in 3D: \code{6} (faces touch, default;
#'   AFNI NN=1), \code{18} (faces or edges touch; AFNI NN=2), or \code{26} (faces, edges,
#'   or corners touch; AFNI NN=3).
#' @param fwe_p Numeric vector of Family-Wise Error Rate alpha levels for which
#'   Z-score thresholds are computed (default: \code{0.05}). When \code{sided = "two.sided"},
#'   thresholding at \code{abs(Z) >= fwe.Z} controls two-sided FWER at \code{fwe_p}.
#' @param legacy Logical; if \code{TRUE}, uses upstream \code{Nh}-point discrete summation
#'   from \code{logpmin} to \code{logpmax} instead of analytical sub-\code{ZestThr}
#'   integration with fixed \code{dh} (default: \code{FALSE}).
#' @param nonstationary Logical; if \code{TRUE}, adjusts cluster sizes for spatial
#'   smoothness non-stationarity using normalized local roughness weights
#'   \code{w(v) = dLh(v) / mean(dLh)} (Worsley et al., 1999; Hayasaka et al., 2004;
#'   Salimi-Khorshidi et al., 2011). Default is \code{FALSE} (stationary cluster voxel counts).
#' @param local_smoothness Optional precomputed local smoothness object for non-stationary
#'   adjustment: either a list returned by \code{smoothest(..., nonstationary = TRUE)},
#'   or a 2D/3D image / file path containing local \code{RPV} (e.g. SPM's
#'   \code{RPV.nii}), \code{dLh}, or positive voxel roughness weights. Supplying
#'   \code{local_smoothness} automatically enables \code{nonstationary = TRUE} unless
#'   explicitly set to \code{FALSE}, and when \code{RPV} or \code{dLh} units are supplied
#'   without \code{Rd}/\code{V}, global Riemannian \code{Rd} and \code{resels} are derived
#'   directly from \code{local_smoothness}.
#' @param local_smoothness_units Units of a 2D/3D image or file supplied via
#'   \code{local_smoothness}: \code{"auto"} (default; detects \code{smoothest()} attributes,
#'   filenames, or mean-1 weights and defaults other positive maps to SPM \code{"rpv"}),
#'   \code{"rpv"} (RESELs per voxel, e.g. SPM \code{RPV.nii}), \code{"dlh"}
#'   (\eqn{|\Lambda|^{1/2}}), or \code{"weights"} (dimensionless relative roughness weights).
#' @param smooth_fwhm Numeric scalar; spatial Gaussian smoothing FWHM in voxels applied
#'   within the mask to regularize voxelwise derivative statistics when estimating
#'   non-stationary smoothness (default: \code{3.0}; Salimi-Khorshidi et al., 2011).
#' @param fwe_method Method used to compute FWER-corrected Z-score thresholds
#'   (\code{fwer0.05.Z} and \code{fwe.Z}):
#'   \itemize{
#'     \item \code{"grf"} (default): Fast analytical Gaussian Random Field Euler
#'       characteristic threshold from \code{number_of_resels} via \code{\link{fwe.p2z}}.
#'       Conservative on pTFCE-enhanced maps because pTFCE smooths supra-threshold peaks.
#'     \item \code{"grf_sim"}: Parametric GRF Monte Carlo simulation (\code{n_perm} null
#'       replicates) synthesized on the exact \code{mask} geometry with stationary
#'       anisotropic or multi-scale non-stationary smoothness and enhanced by the exact
#'       same pTFCE operator.
#'     \item \code{"sign_flip"}: Nonparametric Rademacher (\eqn{\pm 1}) sign-flipping
#'       calibration (\code{n_perm} replicates) across 4D subject-level contrasts
#'       (\code{copes}) or 4D GLM residuals (\code{residual}), preserving empirical
#'       spatial covariance, non-stationarity, and heavy-tailed noise structure.
#'   }
#'   Supplying \code{copes} automatically selects \code{fwe_method = "sign_flip"} unless
#'   \code{fwe_method} is explicitly specified.
#' @param n_perm Integer number of Monte Carlo replicates used when
#'   \code{fwe_method = "grf_sim"} or \code{"sign_flip"} (default: \code{500L}).
#' @param copes Optional 4D image (or NIfTI file path) of independent subject-level
#'   contrast / COPE volumes for group-level Rademacher sign-flipping FWER calibration
#'   (and automatic 4D residual smoothness estimation when \code{residual} is omitted).
#'
#' @details
#' By default (\code{legacy = FALSE}, \code{Nh = NULL}), \code{ptfce()} evaluates the
#' unenhanced sub-\code{ZestThr} interval \code{[0, -log(Phi(-ZestThr))]} analytically
#' and integrates cluster-enhanced log-probabilities above \code{ZestThr} using a fixed
#' step \code{dh = 0.05}. This eliminates two artifacts of the original discrete grid:
#' (i) sensitivity of moderate clusters to distant high-Z peaks via \code{logpmax}, and
#' (ii) the collapse of small-magnitude voxels around \code{Z = 0}.
#'
#' When \code{sided = "two.sided"}, positive (\code{Z > 0}) and negative (\code{Z < 0})
#' excursion sets are enhanced independently on identical threshold grids:
#' \itemize{
#'   \item \code{$Z} contains signed enhanced Z-scores (\code{> 0} where input \code{Z > 0},
#'     \code{< 0} where input \code{Z < 0}, \code{0} outside mask). Thresholding at
#'     \code{abs(res$Z) >= res$fwer0.05.Z} controls two-sided FWER at \code{alpha = 0.05}.
#'   \item \code{$p} contains two-sided uncorrected enhanced p-values \code{2 * Phi(-abs(Z))}.
#'   \item \code{$logp} contains \code{-log(p)}.
#' }
#'
#' When \code{nonstationary = TRUE} (ideally paired with 4D \code{residual} and \code{dof},
#' 4D \code{copes}, or a precomputed \code{local_smoothness} map), each supra-threshold
#' cluster's extent is measured in local resel-scaled voxel units
#' \eqn{\tilde{c}(\mathcal{C}) = \sum_{v \in \mathcal{C}} w(v)} rather than raw voxel count.
#' Smooth brain regions (\eqn{w(v) < 1}) are appropriately tempered against noise-blob
#' inflation, while rougher regions (\eqn{w(v) > 1}, e.g. subcortex) gain sensitivity for
#' compact true signals.
#'
#' Because the pTFCE-enhanced field is smoother than the raw input Z-map, the classical
#' raw-field Euler characteristic formula (\code{fwe_method = "grf"}) is conservative
#' (typical empirical FWER \eqn{\approx 0.005\text{--}0.015} at nominal \eqn{\alpha = 0.05}).
#' Setting \code{fwe_method = "grf_sim"} or \code{fwe_method = "sign_flip"} calibrates the
#' exact distribution of \eqn{\max_{v \in \mathcal{M}} \tilde{Z}(v)} (or
#' \eqn{\max_{v \in \mathcal{M}} |\tilde{Z}(v)|} in \code{"two.sided"} mode) via a fast
#' lookup-table-accelerated C++ engine, recovering full nominal statistical power.
#'
#' @return An object of class \code{"ptfce"} (a list) containing:
#' \itemize{
#'   \item \strong{p}: Uncorrected pTFCE-enhanced p-values (same image class as \code{img}).
#'   \item \strong{logp}: Negative natural logarithm \code{-log(p)} of enhanced p-values.
#'   \item \strong{Z}: Enhanced Z-score image (signed when \code{sided = "two.sided"}).
#'   \item \strong{number_of_resels}: Estimated number of resels in the mask.
#'   \item \strong{fwer0.05.Z}: FWER Z-score threshold at \code{alpha = 0.05}.
#'   \item \strong{fwe.Z}: Named numeric vector of FWER Z-score thresholds corresponding
#'     to \code{fwe_p}.
#'   \item \strong{sided}: Enhancement sidedness used.
#'   \item \strong{nonstationary}: Logical indicating whether non-stationary smoothness
#'     adjustment was applied.
#'   \item \strong{fwe_method}: FWER calibration method used (\code{"grf"},
#'     \code{"grf_sim"}, or \code{"sign_flip"}).
#'   \item \strong{null_max_Z}: Numeric vector of length \code{n_perm} containing the
#'     simulated null maximum enhanced Z-scores when \code{fwe_method != "grf"}
#'     (\code{NULL} when \code{fwe_method = "grf"}).
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard bi-sided fMRI enhancement using 4D residuals for smoothness:
#' res <- ptfce("zstat1.nii.gz", mask = "mask.nii.gz",
#'              residual = "res4d.nii.gz", dof = 42, sided = "two.sided")
#'
#' # Non-stationary pTFCE with parametric GRF simulation FWER calibration:
#' res_sim <- ptfce("zstat1.nii.gz", mask = "mask.nii.gz",
#'                  residual = "res4d.nii.gz", dof = 42,
#'                  sided = "two.sided", nonstationary = TRUE,
#'                  fwe_method = "grf_sim", n_perm = 500)
#'
#' # Group-level non-stationary pTFCE with Rademacher sign-flipping calibration:
#' res_sf <- ptfce("zstat1.nii.gz", mask = "mask.nii.gz",
#'                 copes = "all_copes_4d.nii.gz",
#'                 sided = "two.sided", nonstationary = TRUE,
#'                 fwe_method = "sign_flip", n_perm = 500)
#' }
ptfce <- function(img,
                  mask = NULL,
                  Rd = NA,
                  V = NA,
                  resels = NA,
                  residual = NULL,
                  dof = NA,
                  logpmin = 0,
                  logpmax = NULL,
                  ZestThr = 1.3,
                  Nh = NULL,
                  verbose = TRUE,
                  sided = c("positive", "two.sided", "negative"),
                  dh = 0.05,
                  connectivity = 6L,
                  fwe_p = 0.05,
                  legacy = FALSE,
                  nonstationary = FALSE,
                  local_smoothness = NULL,
                  local_smoothness_units = c("auto", "rpv", "dlh", "weights"),
                  smooth_fwhm = 3.0,
                  fwe_method = c("grf", "grf_sim", "sign_flip"),
                  n_perm = 500L,
                  copes = NULL) {
  sided <- match.arg(sided)
  local_smoothness_units <- match.arg(local_smoothness_units)
  if (!is_missing_arg(copes) && missing(fwe_method)) {
    fwe_method <- "sign_flip"
  } else {
    fwe_method <- match.arg(fwe_method)
  }
  if (!is_missing_arg(local_smoothness) && missing(nonstationary)) {
    nonstationary <- TRUE
  }
  nonstationary <- isTRUE(nonstationary)
  if (!is.numeric(smooth_fwhm) || length(smooth_fwhm) != 1L || !is.finite(smooth_fwhm) || smooth_fwhm < 0) {
    stop("`smooth_fwhm` must be a non-negative finite number.", call. = FALSE)
  }
  if (!is.numeric(fwe_p) || length(fwe_p) < 1L || any(is.na(fwe_p)) || any(fwe_p <= 0 | fwe_p >= 1)) {
    stop("`fwe_p` must be a numeric vector of alpha levels strictly in (0, 1).", call. = FALSE)
  }
  if (fwe_method != "grf") {
    n_perm <- as.integer(n_perm)[1L]
    if (is.na(n_perm) || n_perm < 19L) {
      stop("`n_perm` must be an integer >= 19 when `fwe_method` is 'grf_sim' or 'sign_flip'.", call. = FALSE)
    }
  }

  connectivity <- as.integer(connectivity)
  if (length(connectivity) != 1L || !connectivity %in% c(6L, 18L, 26L)) {
    stop("`connectivity` must be 6, 18, or 26.", call. = FALSE)
  }
  if (!is.numeric(ZestThr) || length(ZestThr) != 1L || !is.finite(ZestThr) || ZestThr <= 0) {
    stop("`ZestThr` must be a positive finite scalar.", call. = FALSE)
  }

  img_meta <- read_volume_input(img, arg_name = "img", allow_4d = FALSE)
  img_arr <- img_meta$data
  d3 <- img_meta$dim3

  if (is_missing_arg(mask)) {
    mask_arr <- array(1.0, dim = d3)
  } else {
    mask_meta <- read_volume_input(mask, arg_name = "mask", allow_4d = FALSE)
    if (!identical(mask_meta$dim3, d3)) {
      stop(sprintf("`mask` dimensions (%s) do not match `img` dimensions (%s).",
                   paste(mask_meta$dim3, collapse = "x"),
                   paste(d3, collapse = "x")), call. = FALSE)
    }
    mask_arr <- mask_meta$data
  }

  inside <- !is.na(mask_arr) & is.finite(mask_arr) & (mask_arr > 0.5)
  if (!any(inside)) {
    stop("`mask` contains zero active voxels (> 0.5).", call. = FALSE)
  }

  # Zero out-of-mask image values first so out-of-brain NAs/Infs (e.g. SPM outputs)
  # neither trigger spurious warnings nor leak into connected components or logpmax.
  img_arr[!inside] <- 0.0

  in_mask_vals <- img_arr[inside]
  if (any(is.na(in_mask_vals))) {
    warning("NAs detected inside mask and replaced with zero!")
    img_arr[inside & is.na(img_arr)] <- 0.0
  }
  if (any(is.infinite(img_arr[inside]))) {
    warning("Infinite values detected inside mask and replaced with zero!")
    img_arr[inside & is.infinite(img_arr)] <- 0.0
  }

  normalize_local_weights <- function(w_raw) {
    w <- numeric(length(img_arr))
    pos <- inside & is.finite(w_raw) & (w_raw > 0)
    if (!any(pos)) {
      stop("`local_smoothness` contains no positive finite values inside `mask`.", call. = FALSE)
    }
    m_raw <- mean(w_raw[pos])
    w[pos] <- pmin(50.0, pmax(0.02, w_raw[pos] / m_raw))
    w[inside & !pos] <- 1.0
    w[inside] <- w[inside] / mean(w[inside])
    w
  }

  mask_clean <- array(as.numeric(inside), dim = d3)
  autosmooth <- FALSE
  voxel_weights <- NULL
  smooth_computed <- NULL
  copes_meta <- NULL
  res_meta <- NULL

  if (!is_missing_arg(copes)) {
    copes_meta <- read_volume_input(copes, arg_name = "copes", allow_4d = TRUE)
    if (!identical(copes_meta$dim3, d3)) {
      stop(sprintf("`copes` spatial dimensions (%s) do not match `img` dimensions (%s).",
                   paste(copes_meta$dim3, collapse = "x"), paste(d3, collapse = "x")), call. = FALSE)
    }
    if (copes_meta$nt < 2L) {
      stop("`copes` must be a 4D image with at least 2 subject/contrast volumes.", call. = FALSE)
    }
    if (fwe_method == "sign_flip" && copes_meta$nt < 6L) {
      stop("`fwe_method = \"sign_flip\"` requires at least 6 volumes in `copes`.", call. = FALSE)
    }
  } else if (fwe_method == "sign_flip") {
    if (is_missing_arg(residual)) {
      stop("`fwe_method = \"sign_flip\"` requires either a 4D `copes` image or a 4D `residual` image.", call. = FALSE)
    }
    res_meta <- read_volume_input(residual, arg_name = "residual", allow_4d = TRUE)
    if (!identical(res_meta$dim3, d3)) {
      stop(sprintf("`residual` spatial dimensions (%s) do not match `img` dimensions (%s).",
                   paste(res_meta$dim3, collapse = "x"), paste(d3, collapse = "x")), call. = FALSE)
    }
    if (res_meta$nt < 6L) {
      stop("`fwe_method = \"sign_flip\"` requires a 4D `residual` image with at least 6 volumes.", call. = FALSE)
    }
  }

  # Pre-parse external 2D/3D `local_smoothness` image/file (if supplied) so global V, Rd,
  # and resels can be derived directly from RPV/dLh maps without falling back to 3D Z-map smoothest()
  ext_loc_meta <- NULL
  loc_units_eff <- local_smoothness_units
  if (!is_missing_arg(local_smoothness) && !is.list(local_smoothness)) {
    ext_loc_meta <- read_volume_input(local_smoothness, arg_name = "local_smoothness", allow_4d = FALSE)
    if (!identical(ext_loc_meta$dim3, d3)) {
      stop(sprintf("`local_smoothness` dimensions (%s) do not match `img` dimensions (%s).",
                   paste(ext_loc_meta$dim3, collapse = "x"), paste(d3, collapse = "x")), call. = FALSE)
    }
    if (loc_units_eff == "auto") {
      tag_u <- attr(local_smoothness, "ptfce_smoothness_units", exact = TRUE)
      if (is.character(tag_u) && length(tag_u) == 1L && tag_u %in% c("rpv", "dlh", "weights")) {
        loc_units_eff <- tag_u
      } else if (is.character(local_smoothness)) {
        bn <- basename(local_smoothness[1L])
        if (grepl("dlh", bn, ignore.case = TRUE)) {
          loc_units_eff <- "dlh"
        } else if (grepl("weight", bn, ignore.case = TRUE)) {
          loc_units_eff <- "weights"
        } else {
          loc_units_eff <- "rpv"
        }
      } else {
        w_tmp <- as.numeric(ext_loc_meta$data)
        pos_tmp <- inside & is.finite(w_tmp) & (w_tmp > 0)
        if (any(pos_tmp) && abs(mean(w_tmp[pos_tmp]) - 1.0) <= 0.05) {
          loc_units_eff <- "weights"
        } else {
          loc_units_eff <- "rpv"
        }
      }
    }
  }

  # If a precomputed smoothest(..., nonstationary = TRUE) list or RPV/dLh map was supplied, reuse its parameters
  if (is.list(local_smoothness) && all(c("volume", "dLh", "resels") %in% names(local_smoothness))) {
    if (is_missing_arg(Rd) || is_missing_arg(V)) {
      autosmooth <- TRUE
      V <- as.numeric(local_smoothness$volume)[1L]
      if (nonstationary && !is.null(local_smoothness$dLh_nonstationary)) {
        Rd <- as.numeric(local_smoothness$dLh_nonstationary)[1L] * V
        resels <- as.numeric(local_smoothness$resels_nonstationary)[1L]
      } else {
        Rd <- as.numeric(local_smoothness$dLh)[1L] * V
        resels <- as.numeric(local_smoothness$resels)[1L]
      }
    }
  } else if (!is.null(ext_loc_meta) && (is_missing_arg(Rd) || is_missing_arg(V)) &&
             is_missing_arg(residual) && is.null(copes_meta) &&
             loc_units_eff %in% c("rpv", "dlh")) {
    w_vec_loc <- as.numeric(ext_loc_meta$data)
    pos_loc <- inside & is.finite(w_vec_loc) & (w_vec_loc > 0)
    if (!any(pos_loc)) {
      stop("`local_smoothness` contains no positive finite values inside `mask`.", call. = FALSE)
    }
    usez_loc <- d3[3L] > 1L
    rpv_scale_loc <- (4.0 * log(2.0))^(if (usez_loc) 1.5 else 1.0)
    V <- as.numeric(sum(inside))
    if (loc_units_eff == "rpv") {
      mean_rpv <- mean(w_vec_loc[pos_loc])
      dlh_ns_loc <- mean_rpv * rpv_scale_loc
      resels <- 1.0 / mean_rpv
    } else {
      dlh_ns_loc <- mean(w_vec_loc[pos_loc])
      resels <- rpv_scale_loc / dlh_ns_loc
    }
    Rd <- V * dlh_ns_loc
    autosmooth <- TRUE
    if (verbose) {
      cat(sprintf("* Derived Riemannian global smoothness from `local_smoothness` (%s)...\n", loc_units_eff))
    }
  }

  if (is_missing_arg(Rd) || is_missing_arg(V)) {
    autosmooth <- TRUE
    need_local_here <- nonstationary && is_missing_arg(local_smoothness)
    if (!is_missing_arg(residual)) {
      res_for_sm <- if (!is.null(res_meta)) res_meta$data else residual
      if (!is_missing_arg(dof)) {
        smooth_computed <- smoothest(
          res_for_sm, mask_clean, dof = dof, verbose = verbose,
          nonstationary = need_local_here, smooth_fwhm = smooth_fwhm
        )
      } else {
        warning("Residual supplied without `dof`: defaulting back to Z-score map smoothness estimation.")
        if (need_local_here) {
          warning("`nonstationary = TRUE` without valid 4D `residual` + `dof`: estimating local smoothness from 3D Z-map is noisy and signal-sensitive; 4D residuals are strongly recommended.")
        }
        smooth_computed <- smoothest(
          img_arr, mask_clean, verbose = verbose,
          nonstationary = need_local_here, smooth_fwhm = if (need_local_here) max(smooth_fwhm, 4.0) else smooth_fwhm
        )
      }
    } else if (!is.null(copes_meta)) {
      dof_copes <- if (!is_missing_arg(dof)) as.numeric(dof)[1L] else as.numeric(copes_meta$nt - 1L)
      smooth_computed <- smoothest(
        copes_meta$data, mask_clean, dof = dof_copes, verbose = verbose,
        nonstationary = need_local_here, smooth_fwhm = smooth_fwhm
      )
    } else {
      if (need_local_here) {
        warning("`nonstationary = TRUE` without 4D `residual` or `local_smoothness`: estimating local smoothness from 3D Z-map is noisy and signal-sensitive; 4D residuals are strongly recommended.")
      }
      smooth_computed <- smoothest(
        img_arr, mask_clean, verbose = verbose,
        nonstationary = need_local_here, smooth_fwhm = if (need_local_here) max(smooth_fwhm, 4.0) else smooth_fwhm
      )
    }
    V <- smooth_computed$volume
    if (nonstationary && !is.null(smooth_computed$dLh_nonstationary)) {
      Rd <- smooth_computed$dLh_nonstationary * V
      resels <- smooth_computed$resels_nonstationary
    } else {
      Rd <- smooth_computed$dLh * V
      resels <- smooth_computed$resels
    }
  } else {
    V <- as.numeric(V)[1L]
    Rd <- as.numeric(Rd)[1L]
    if (!is.finite(V) || V <= 0 || !is.finite(Rd) || Rd <= 0) {
      stop("`V` and `Rd` must be positive finite numbers.", call. = FALSE)
    }
  }

  if (nonstationary) {
    if (!is_missing_arg(local_smoothness)) {
      if (is.list(local_smoothness)) {
        w_src <- if (!is.null(local_smoothness$voxel_weights)) {
          local_smoothness$voxel_weights
        } else if (!is.null(local_smoothness$dLh_map)) {
          local_smoothness$dLh_map
        } else if (!is.null(local_smoothness$rpv_map)) {
          local_smoothness$rpv_map
        } else {
          stop("`local_smoothness` list must contain `voxel_weights`, `dLh_map`, or `rpv_map`.", call. = FALSE)
        }
        w_meta <- read_volume_input(w_src, arg_name = "local_smoothness", allow_4d = FALSE)
        if (!identical(w_meta$dim3, d3)) {
          stop(sprintf("`local_smoothness` dimensions (%s) do not match `img` dimensions (%s).",
                       paste(w_meta$dim3, collapse = "x"), paste(d3, collapse = "x")), call. = FALSE)
        }
      } else {
        w_meta <- ext_loc_meta
      }
      voxel_weights <- normalize_local_weights(as.numeric(w_meta$data))
    } else if (!is.null(smooth_computed) && !is.null(smooth_computed$voxel_weights)) {
      voxel_weights <- normalize_local_weights(as.numeric(smooth_computed$voxel_weights))
    } else {
      # User supplied explicit V and Rd, plus nonstationary = TRUE without local_smoothness
      if (!is_missing_arg(residual) && !is_missing_arg(dof)) {
        res_for_sm <- if (!is.null(res_meta)) res_meta$data else residual
        sm_loc <- smoothest(res_for_sm, mask_clean, dof = dof, verbose = verbose,
                            nonstationary = TRUE, smooth_fwhm = smooth_fwhm)
      } else if (!is.null(copes_meta)) {
        dof_copes <- if (!is_missing_arg(dof)) as.numeric(dof)[1L] else as.numeric(copes_meta$nt - 1L)
        sm_loc <- smoothest(copes_meta$data, mask_clean, dof = dof_copes, verbose = verbose,
                            nonstationary = TRUE, smooth_fwhm = smooth_fwhm)
      } else {
        warning("`nonstationary = TRUE` without 4D `residual` or `local_smoothness`: estimating local smoothness from 3D Z-map is noisy and signal-sensitive; 4D residuals are strongly recommended.")
        sm_loc <- smoothest(img_arr, mask_clean, verbose = verbose,
                            nonstationary = TRUE, smooth_fwhm = max(smooth_fwhm, 4.0))
      }
      voxel_weights <- normalize_local_weights(as.numeric(sm_loc$voxel_weights))
      if (is.null(smooth_computed)) smooth_computed <- sm_loc
    }
  }

  use_legacy <- isTRUE(legacy) || !is.null(Nh) || !is.null(logpmax)
  in_vals <- img_arr[inside]
  z_peak <- switch(
    sided,
    "positive" = max(in_vals),
    "negative" = max(-in_vals),
    "two.sided" = max(abs(in_vals))
  )

  if (use_legacy) {
    Nh_eff <- if (is.null(Nh)) 100L else as.integer(Nh)[1L]
    if (is.na(Nh_eff) || Nh_eff < 2L) {
      stop("`Nh` must be an integer >= 2.", call. = FALSE)
    }
    if (is.null(logpmax)) {
      logpmax_eff <- -stats::pnorm(z_peak, lower.tail = FALSE, log.p = TRUE)
    } else {
      logpmax_eff <- as.numeric(logpmax)[1L]
    }
    if (!is.finite(logpmax_eff) || logpmax_eff <= logpmin) {
      logpmax_eff <- logpmin + 1e-6
    }
    logp_thres <- seq(logpmin, logpmax_eff, length.out = Nh_eff)
    dh_legacy <- logp_thres[2L] - logp_thres[1L]
    threshs <- stats::qnorm(-logp_thres, lower.tail = FALSE, log.p = TRUE)
    if (logpmin == 0) {
      threshs[1L] <- -9999999999
    }
    weights <- rep(1.0, Nh_eff)
    analytical_subthr <- FALSE
  } else {
    if (!is.numeric(dh) || length(dh) != 1L || !is.finite(dh) || dh <= 0) {
      stop("`dh` must be a positive finite number.", call. = FALSE)
    }
    x_zest <- -stats::pnorm(ZestThr, lower.tail = FALSE, log.p = TRUE)
    x_max <- -stats::pnorm(z_peak, lower.tail = FALSE, log.p = TRUE)

    if (!is.finite(x_max) || x_max <= x_zest) {
      threshs <- numeric(0)
      weights <- numeric(0)
    } else {
      # Fixed step `dh` across [x_zest, 100] (up to Z ~ 14), with bounded adaptive
      # panels above 100 so extreme Z > 14 peaks never explode panel count.
      x_fine_cap <- min(x_max, 100.0)
      n_lo <- max(1L, as.integer(ceiling((x_fine_cap - x_zest) / dh)))
      m_grid <- x_zest + (seq_len(n_lo) - 0.5) * dh
      w_grid <- rep(dh, n_lo)
      top_lo <- x_zest + n_lo * dh

      if (x_max > top_lo) {
        n_hi <- 200L
        dh_hi <- (x_max - top_lo) / n_hi
        m_hi <- top_lo + (seq_len(n_hi) - 0.5) * dh_hi
        w_hi <- rep(dh_hi, n_hi)
        m_grid <- c(m_grid, m_hi)
        w_grid <- c(w_grid, w_hi)
      }
      threshs <- stats::qnorm(-m_grid, lower.tail = FALSE, log.p = TRUE)
      weights <- w_grid
    }
    dh_legacy <- dh
    analytical_subthr <- TRUE
  }

  if (verbose) {
    mode_lbl <- if (nonstationary) paste0(sided, ", non-stationary") else sided
    cat(sprintf("* Performing pTFCE (%s)...\n", mode_lbl))
  }

  n_total <- length(img_arr)
  p_arr <- rep(1.0, n_total)
  logp_arr <- numeric(n_total)
  Z_arr <- numeric(n_total)

  run_single_tail <- function(tail_img) {
    nlp <- ptfce_tail_cpp(
      img = tail_img,
      mask = inside,
      nx = d3[1L],
      ny = d3[2L],
      nz = d3[3L],
      V = V,
      Rd = Rd,
      threshs = threshs,
      weights = weights,
      dh_legacy = dh_legacy,
      ZestThr = ZestThr,
      connectivity = connectivity,
      analytical_subthr = analytical_subthr,
      voxel_weights = voxel_weights
    )
    nlp
  }

  if (sided %in% c("positive", "negative")) {
    work_img <- if (sided == "positive") img_arr else -img_arr
    nlp <- run_single_tail(work_img)

    if (use_legacy) {
      logp_arr[inside] <- nlp[inside]
      p_in <- exp(-nlp[inside])
      p_in[p_in == 0] <- .Machine$double.xmin
      p_in[p_in == 1] <- 1 - .Machine$double.neg.eps
      p_arr[inside] <- p_in

      z_in <- numeric(sum(inside))
      pos_nlp <- nlp[inside] > 0
      z_in[pos_nlp] <- stats::qnorm(-nlp[inside][pos_nlp], lower.tail = FALSE, log.p = TRUE)
      z_in[!pos_nlp] <- stats::qnorm(1 - .Machine$double.neg.eps, lower.tail = FALSE)
      Z_arr[inside] <- z_in
    } else {
      # Modern mode: exact analytical preservation for voxels <= ZestThr
      sub_mask <- inside & (work_img <= ZestThr)
      sup_mask <- inside & (work_img > ZestThr)

      if (any(sub_mask)) {
        w_sub <- work_img[sub_mask]
        Z_arr[sub_mask] <- w_sub
        logp_arr[sub_mask] <- -stats::pnorm(w_sub, lower.tail = FALSE, log.p = TRUE)
        p_arr[sub_mask] <- stats::pnorm(w_sub, lower.tail = FALSE)
      }
      if (any(sup_mask)) {
        nlp_sup <- nlp[sup_mask]
        logp_arr[sup_mask] <- nlp_sup
        p_sup <- exp(-nlp_sup)
        p_sup[p_sup == 0] <- .Machine$double.xmin
        p_arr[sup_mask] <- p_sup
        Z_arr[sup_mask] <- stats::qnorm(-nlp_sup, lower.tail = FALSE, log.p = TRUE)
      }
    }
  } else {
    # Native bi-sided ("two.sided") mode
    nlp_pos <- run_single_tail(img_arr)
    nlp_neg <- run_single_tail(-img_arr)

    pos_vox <- inside & (img_arr > 0)
    neg_vox <- inside & (img_arr < 0)

    if (use_legacy) {
      if (any(pos_vox)) {
        nlp <- nlp_pos[pos_vox]
        z_p <- ifelse(nlp > 0, stats::qnorm(-nlp, lower.tail = FALSE, log.p = TRUE), -Inf)
        z_p <- pmax(0.0, z_p)
        Z_arr[pos_vox] <- z_p
        lp2 <- pmax(0.0, nlp - log(2))
        logp_arr[pos_vox] <- lp2
        p2 <- exp(-lp2)
        p2[p2 == 0] <- .Machine$double.xmin
        p_arr[pos_vox] <- p2
      }
      if (any(neg_vox)) {
        nlp <- nlp_neg[neg_vox]
        z_n <- ifelse(nlp > 0, stats::qnorm(-nlp, lower.tail = FALSE, log.p = TRUE), -Inf)
        z_n <- pmax(0.0, z_n)
        Z_arr[neg_vox] <- -z_n
        lp2 <- pmax(0.0, nlp - log(2))
        logp_arr[neg_vox] <- lp2
        p2 <- exp(-lp2)
        p2[p2 == 0] <- .Machine$double.xmin
        p_arr[neg_vox] <- p2
      }
    } else {
      # Modern bi-sided mode: exact identity on [-ZestThr, ZestThr], smooth boost beyond
      sub_vox <- inside & (abs(img_arr) <= ZestThr) & (img_arr != 0)
      if (any(sub_vox)) {
        z_sub <- img_arr[sub_vox]
        Z_arr[sub_vox] <- z_sub
        nlp_1s <- -stats::pnorm(abs(z_sub), lower.tail = FALSE, log.p = TRUE)
        lp2 <- pmax(0.0, nlp_1s - log(2))
        logp_arr[sub_vox] <- lp2
        p_arr[sub_vox] <- exp(-lp2)
      }

      pos_sup <- inside & (img_arr > ZestThr)
      if (any(pos_sup)) {
        nlp <- nlp_pos[pos_sup]
        Z_arr[pos_sup] <- stats::qnorm(-nlp, lower.tail = FALSE, log.p = TRUE)
        lp2 <- pmax(0.0, nlp - log(2))
        logp_arr[pos_sup] <- lp2
        p2 <- exp(-lp2)
        p2[p2 == 0] <- .Machine$double.xmin
        p_arr[pos_sup] <- p2
      }

      neg_sup <- inside & (img_arr < -ZestThr)
      if (any(neg_sup)) {
        nlp <- nlp_neg[neg_sup]
        Z_arr[neg_sup] <- -stats::qnorm(-nlp, lower.tail = FALSE, log.p = TRUE)
        lp2 <- pmax(0.0, nlp - log(2))
        logp_arr[neg_sup] <- lp2
        p2 <- exp(-lp2)
        p2[p2 == 0] <- .Machine$double.xmin
        p_arr[neg_sup] <- p2
      }
    }
  }

  # Compute number_of_resels (when available) and FWER Z-score thresholds
  if (autosmooth || !is_missing_arg(resels)) {
    resels_val <- as.numeric(resels)[1L]
    number_of_resels <- V / resels_val
  } else {
    number_of_resels <- NA_real_
  }

  null_max_Z <- NULL

  if (fwe_method == "grf") {
    if (!is.na(number_of_resels)) {
      alpha_scale <- if (sided == "two.sided") 0.5 else 1.0
      dim_fwe <- if (d3[3L] > 1L) 3L else 2L
      fwer0.05.Z <- fwe.p2z(number_of_resels, 0.05 * alpha_scale, D = dim_fwe)
      fwe_Z <- stats::setNames(
        fwe.p2z(number_of_resels, fwe_p * alpha_scale, D = dim_fwe),
        format(fwe_p, trim = TRUE)
      )
    } else {
      warning("For GRF-based FWER correction, please specify `resels` or allow automatic smoothness estimation.")
      fwer0.05.Z <- NA_real_
      fwe_Z <- stats::setNames(rep(NA_real_, length(fwe_p)), format(fwe_p, trim = TRUE))
    }
  } else {
    sided_int <- switch(sided, "positive" = 0L, "negative" = 1L, "two.sided" = 2L)
    Nh_pass <- if (use_legacy) Nh_eff else 100L
    lpmax_pass <- if (use_legacy && !is.null(logpmax)) as.numeric(logpmax)[1L] else NA_real_
    dh_pass <- if (use_legacy) dh_legacy else as.numeric(dh)[1L]

    if (fwe_method == "grf_sim") {
      if (verbose) {
        cat(sprintf("* Calibrating FWER threshold via parametric GRF simulation (B = %d)...\n", n_perm))
      }
      sigmasq_base <- NULL
      if (!is.null(smooth_computed) && !is.null(smooth_computed$sigmasq)) {
        sigmasq_base <- as.numeric(smooth_computed$sigmasq)
      } else if (is.list(local_smoothness) && !is.null(local_smoothness$sigmasq)) {
        sigmasq_base <- as.numeric(local_smoothness$sigmasq)
      }
      if (is.null(sigmasq_base) || length(sigmasq_base) < 3L || !all(is.finite(sigmasq_base[1:2]))) {
        dlh_iso <- Rd / V
        if (d3[3L] > 1L) {
          s2_iso <- 0.5 * (dlh_iso^(-2.0 / 3.0))
          sigmasq_base <- c(s2_iso, s2_iso, s2_iso)
        } else {
          s2_iso <- 0.5 / dlh_iso
          sigmasq_base <- c(s2_iso, s2_iso, 0.0)
        }
      }

      null_max_Z <- calibrate_fwer_grf_cpp(
        mask = inside,
        nx = d3[1L],
        ny = d3[2L],
        nz = d3[3L],
        V = V,
        Rd = Rd,
        sigmasq_base = sigmasq_base,
        ZestThr = ZestThr,
        connectivity = connectivity,
        sided_mode = sided_int,
        use_legacy = use_legacy,
        Nh_eff = Nh_pass,
        logpmin = as.numeric(logpmin)[1L],
        logpmax_fixed = lpmax_pass,
        dh = dh_pass,
        n_perm = n_perm,
        voxel_weights = voxel_weights
      )
    } else {
      # fwe_method == "sign_flip"
      use_copes_sf <- !is.null(copes_meta)
      sf_meta <- if (use_copes_sf) copes_meta else res_meta
      if (verbose) {
        src_lbl <- if (use_copes_sf) "4D copes" else "4D residuals"
        cat(sprintf("* Calibrating FWER threshold via Rademacher sign-flipping (%s, B = %d)...\n",
                    src_lbl, n_perm))
      }

      dof_sf_hint <- if (!is_missing_arg(dof)) as.numeric(dof)[1L] else NA_real_
      sf_res <- calibrate_fwer_signflip_cpp(
        data4d = sf_meta$data,
        mask = inside,
        nx = d3[1L],
        ny = d3[2L],
        nz = d3[3L],
        nt = sf_meta$nt,
        center_residuals = !use_copes_sf,
        V = V,
        Rd = Rd,
        ZestThr = ZestThr,
        connectivity = connectivity,
        sided_mode = sided_int,
        use_legacy = use_legacy,
        Nh_eff = Nh_pass,
        logpmin = as.numeric(logpmin)[1L],
        logpmax_fixed = lpmax_pass,
        dh = dh_pass,
        n_perm = n_perm,
        voxel_weights = voxel_weights,
        dof_hint = dof_sf_hint
      )

      if (!use_copes_sf && is.finite(sf_res$mean_temp_ar1) && abs(sf_res$mean_temp_ar1) > 0.25) {
        warning(sprintf(
          paste0("Substantial lag-1 autocorrelation (r = %.2f) detected along the 4th dimension ",
                 "of `residual`. Rademacher sign-flipping assumes independent/exchangeable volumes ",
                 "(e.g., group-level subjects or prewhitened residuals); consider `fwe_method = \"grf_sim\"` ",
                 "for unwhitened time series."),
          sf_res$mean_temp_ar1
        ))
      }
      null_max_Z <- sf_res$null_max_Z
    }

    m_sorted <- sort(as.numeric(null_max_Z), decreasing = FALSE)
    B_eff <- length(m_sorted)
    if (B_eff < 1L || anyNA(m_sorted)) {
      stop("Empirical FWER calibration produced invalid null maximum statistics.", call. = FALSE)
    }

    # Exact finite-sample Monte Carlo quantile M_(ceiling((1 - alpha) * (B + 1)))
    # guarantees P_H0(max_Z > thr) <= alpha by exchangeability (Phipson & Smyth, 2010).
    emp_p2z <- function(alphas, warn_small = TRUE) {
      res_z <- numeric(length(alphas))
      for (idx_a in seq_along(alphas)) {
        a <- alphas[idx_a]
        if ((B_eff + 1.0) * a < 1.0 - 1e-12) {
          if (warn_small) {
            warning(sprintf(
              paste0("`n_perm = %d` is too small for `fwe_p = %g` ((n_perm + 1) * alpha < 1); ",
                     "using maximum simulated null statistic. Increase `n_perm >= %d`."),
              B_eff, a, as.integer(ceiling(1.0 / a - 1.0))
            ))
          }
          res_z[idx_a] <- m_sorted[B_eff]
        } else {
          j <- as.integer(ceiling((1.0 - a) * (B_eff + 1.0) - 1e-12))
          j <- min(B_eff, max(1L, j))
          res_z[idx_a] <- m_sorted[j]
        }
      }
      res_z
    }

    fwer0.05.Z <- emp_p2z(0.05, warn_small = FALSE)
    fwe_Z <- stats::setNames(emp_p2z(fwe_p, warn_small = TRUE), format(fwe_p, trim = TRUE))
  }

  out <- list(
    p = wrap_volume_output(p_arr, img_meta),
    logp = wrap_volume_output(logp_arr, img_meta),
    Z = wrap_volume_output(Z_arr, img_meta),
    number_of_resels = number_of_resels,
    fwer0.05.Z = fwer0.05.Z,
    fwe.Z = fwe_Z,
    sided = sided,
    nonstationary = nonstationary,
    fwe_method = fwe_method,
    null_max_Z = null_max_Z
  )
  class(out) <- "ptfce"
  out
}

#' @export
print.ptfce <- function(x, ...) {
  ns_lbl <- if (isTRUE(x$nonstationary)) ", non-stationary" else ""
  fwe_lbl <- if (!is.null(x$fwe_method) && x$fwe_method != "grf") {
    sprintf(" [%s, B=%d]", x$fwe_method, length(x$null_max_Z))
  } else {
    " [GRF]"
  }
  cat("pTFCE result (", x$sided, ns_lbl, ")\n", sep = "")
  cat("  Number of resels : ", format(x$number_of_resels, digits = 5), "\n", sep = "")
  cat("  FWER 0.05 Z-thr  : ", format(x$fwer0.05.Z, digits = 5), fwe_lbl, "\n", sep = "")
  invisible(x)
}

#' Convert a Student's t-statistic image or vector to Z-scores in log-space
#'
#' Accurately maps t-statistics with \code{dof} degrees of freedom to standard
#' normal Z-scores without overflow at extreme t-values (\code{|t| >> 15}).
#'
#' @param img A numeric vector, array, \code{RNifti} image, \code{oro.nifti} object,
#'   or NIfTI file path containing t-statistics.
#' @param dof Positive degrees of freedom.
#'
#' @return Z-score object of the same class and dimensions as \code{img}.
#' @export
t2z <- function(img, dof) {
  if (!is.numeric(dof) || length(dof) != 1L || !is.finite(dof) || dof <= 0) {
    stop("`dof` must be a positive finite number.", call. = FALSE)
  }
  if (is.numeric(img) && is.null(dim(img))) {
    sgn <- sign(img)
    lp <- stats::pt(abs(img), df = dof, lower.tail = FALSE, log.p = TRUE)
    z <- -stats::qnorm(lp, lower.tail = TRUE, log.p = TRUE) * sgn
    z[img == 0] <- 0
    return(z)
  }
  meta <- read_volume_input(img, arg_name = "img", allow_4d = TRUE)
  arr <- meta$data
  sgn <- sign(arr)
  lp <- stats::pt(abs(arr), df = dof, lower.tail = FALSE, log.p = TRUE)
  z_arr <- -stats::qnorm(lp, lower.tail = TRUE, log.p = TRUE) * sgn
  z_arr[arr == 0] <- 0
  if (meta$nt > 1L) {
    dim(z_arr) <- meta$orig_dim
    if (identical(meta$kind, "rnifti")) {
      return(RNifti::updateNifti(z_arr, template = meta$template))
    }
    return(z_arr)
  }
  wrap_volume_output(z_arr, meta)
}

# Helper functions

#' Aggregate logpvals by equidistant incremental logarithmic probability aggregation
#'
#' @param logpvals Vector of enhanced \code{-log(p)} values across thresholds.
#' @param d Step size \code{delta -log(p)}.
#'
#' @return Aggregated \code{-log(p)} value.
#' @export
aggregate_logpvals <- function(logpvals, d) {
  logpvals[is.infinite(logpvals) & logpvals > 0] <- 745
  s <- max(0, sum(logpvals))
  0.5 * (sqrt(d * (8 * s + d)) - d)
}

#' Expected cluster size at threshold h in image of volume V with roughness-volume Rd
#'
#' Ported from FSL Gaussian Random Field cluster theory.
#'
#' @param h Z-score cluster-forming threshold(s).
#' @param V Number of voxels in mask.
#' @param Rd Roughness-volume product (\code{dLh * V}).
#' @param D Spatial dimension (\code{3} for 3D volumes, default; \code{2} for 2D planar images).
#'
#' @return Expected value of cluster size in voxels.
#' @export
Es <- function(h, V, Rd, D = 3L) {
  es_cpp(as.numeric(h), as.numeric(V)[1L], as.numeric(Rd)[1L], as.integer(D)[1L])
}

dvox <- function(h) {
  stats::dnorm(h)
}

pvox <- function(h) {
  stats::pnorm(h, lower.tail = FALSE)
}

dclust <- function(h, V, Rd, c, ZestThr = 1.3, CONST = 1, D = 3L) {
  if (length(c) == 1L && is.na(c)) {
    return(dvox(h))
  }
  dim_eff <- if (as.integer(D)[1L] == 2L) 2L else 3L
  beta_d <- if (dim_eff == 2L) 1.0 else (2.0 / 3.0)
  gam_d <- if (dim_eff == 2L) 1.0 else gamma(2.5)
  h_num <- as.numeric(h)
  lam <- (Es(h_num, V, Rd, D = dim_eff) / gam_d)^(-beta_d)
  log_d <- log(lam) - lam * (c^beta_d) + log(CONST)
  out <- exp(log_d)
  out[h_num < ZestThr] <- 0
  out
}

dvox.clust <- function(h, V, Rd, c, ZestThr = 1.3, D = 3L) {
  if (length(c) == 1L && is.na(c)) {
    return(dvox(h))
  }
  dvox_clust_cpp(as.numeric(h), as.numeric(V)[1L], as.numeric(Rd)[1L],
                 as.numeric(c)[1L], ZestThr, as.integer(D)[1L])
}

pvox.clust <- function(V, Rd, c, actH, ZestThr = 1.3, D = 3L) {
  nlp <- neglog_pvox_clust_cpp(
    as.numeric(V)[1L],
    as.numeric(Rd)[1L],
    as.numeric(c),
    as.numeric(actH),
    as.numeric(ZestThr)[1L],
    as.integer(D)[1L]
  )
  p <- exp(-nlp)
  p[p == 0] <- .Machine$double.xmin
  p
}

#' Estimate global (and optionally local non-stationary) image smoothness
#'
#' Fast C++ port of FSL's \code{smoothest} tool (Jenkinson, FMRIB Tech Report TR00DF1)
#' supporting both 2D/3D Z-score maps and 4D GLM residual time-series images, with
#' optional regularized voxelwise local roughness estimation (\code{dLh_map}, \code{rpv_map},
#' and normalized \code{voxel_weights}) for non-stationary cluster inference
#' (Worsley et al., 1999; Salimi-Khorshidi et al., 2011).
#'
#' @param img A 2D/3D Z-score image or 4D residual image (\code{RNifti}, \code{oro.nifti},
#'   \code{array}, or file path).
#' @param mask Spatial mask image (\code{RNifti}, \code{oro.nifti}, \code{array}, or file path).
#' @param dof Residual degrees of freedom (required when \code{img} is a 4D residual image).
#' @param verbose Logical; print diagnostic messages if \code{TRUE} (default).
#' @param nonstationary Logical; if \code{TRUE}, also compute spatially regularized
#'   voxelwise \code{dLh_map}, \code{rpv_map} (RESELs per voxel), and mean-1 normalized
#'   \code{voxel_weights} inside the mask (default: \code{FALSE}).
#' @param smooth_fwhm Spatial Gaussian smoothing FWHM in voxels applied within the mask
#'   to regularize voxelwise derivative sufficient statistics when \code{nonstationary = TRUE}
#'   (default: \code{3.0}).
#'
#' @return A list containing:
#' \itemize{
#'   \item \strong{volume}: Mask volume in voxels.
#'   \item \strong{sigmasq}: Estimated Gaussian smoothness variance in x, y, z (voxels^2).
#'   \item \strong{FWHM}: Full-width at half-maximum in x, y, z (voxels).
#'   \item \strong{dLh}: Determinant of \eqn{\Lambda^{1/2}} (\code{voxels^-3}).
#'   \item \strong{resels}: Resel size in voxels per resel.
#'   \item \strong{dLh_map}, \strong{rpv_map}, \strong{voxel_weights}: Included when
#'     \code{nonstationary = TRUE}, preserving the spatial image class of \code{mask}.
#' }
#' @export
smoothest <- function(img, mask, dof = NA, verbose = TRUE,
                      nonstationary = FALSE, smooth_fwhm = 3.0) {
  if (verbose) cat("* Estimating smoothness based on the data...\n")
  nonstationary <- isTRUE(nonstationary)
  if (!is.numeric(smooth_fwhm) || length(smooth_fwhm) != 1L || !is.finite(smooth_fwhm) || smooth_fwhm < 0) {
    stop("`smooth_fwhm` must be a non-negative finite number.", call. = FALSE)
  }

  img_meta <- read_volume_input(img, arg_name = "img", allow_4d = TRUE)
  d3 <- img_meta$dim3
  nt <- img_meta$nt

  if (nt > 1L && is_missing_arg(dof)) {
    stop("`dof` must be provided when estimating smoothness from a 4D residual image.", call. = FALSE)
  }

  mask_meta <- read_volume_input(mask, arg_name = "mask", allow_4d = FALSE)
  if (!identical(mask_meta$dim3, d3)) {
    stop(sprintf("`mask` dimensions (%s) do not match `img` spatial dimensions (%s).",
                 paste(mask_meta$dim3, collapse = "x"),
                 paste(d3, collapse = "x")), call. = FALSE)
  }

  if (d3[3L] <= 1L) {
    warning("Using 2D image mode.")
  }

  dof_corr <- if (nt > 1L) interpolate(as.numeric(dof)[1L]) else 1.0

  cpp_res <- smoothest_cpp(
    img = img_meta$data,
    mask = mask_meta$data,
    nx = d3[1L],
    ny = d3[2L],
    nz = d3[3L],
    nt = nt,
    compute_local = nonstationary,
    smooth_fwhm = smooth_fwhm,
    dof_correction = dof_corr
  )

  usez <- cpp_res$usez
  SSminus <- cpp_res$SSminus
  S2 <- cpp_res$S2

  if (cpp_res$N <= 0) {
    stop("Mask contains too few contiguous interior voxels to estimate spatial smoothness.", call. = FALSE)
  }

  if (nt > 1L && verbose) {
    v <- as.numeric(dof)[1L]
    cat(sprintf("  Non-edge voxels = %.0f, (dof - 2)/(dof - 1) = %.6f\n",
                cpp_res$N, (v - 2) / (v - 1)))
  }

  axes <- c("X", "Y", "Z")
  n_axes <- if (usez) 3L else 2L
  for (ax in seq_len(n_axes)) {
    if (SSminus[ax] >= 0.99999999 * S2[ax]) {
      SSminus[ax] <- 0.99999 * S2[ax]
      warning(sprintf("Extreme smoothness detected in %s - possibly biased global estimate.", axes[ax]))
    }
  }

  sigmasq <- c(x = NA_real_, y = NA_real_, z = NA_real_)
  sigmasq[1L] <- -1.0 / (4.0 * log(abs(SSminus[1L] / S2[1L])))
  sigmasq[2L] <- -1.0 / (4.0 * log(abs(SSminus[2L] / S2[2L])))
  sigmasq[3L] <- if (usez) -1.0 / (4.0 * log(abs(SSminus[3L] / S2[3L]))) else 0.0

  if (usez) {
    dLh <- ((sigmasq[1L] * sigmasq[2L] * sigmasq[3L])^-0.5) * (8.0^-0.5)
  } else {
    dLh <- ((sigmasq[1L] * sigmasq[2L])^-0.5) * (4.0^-0.5)
  }

  if (nt > 1L) {
    dLh <- dLh * dof_corr
  }

  FWHM <- c(
    x = sqrt(8.0 * log(2.0) * sigmasq[1L]),
    y = sqrt(8.0 * log(2.0) * sigmasq[2L]),
    z = if (usez) sqrt(8.0 * log(2.0) * sigmasq[3L]) else 0.0
  )

  resels <- FWHM[1L] * FWHM[2L] * (if (usez) FWHM[3L] else 1.0)

  out <- list(
    volume = cpp_res$volume,
    sigmasq = sigmasq,
    FWHM = FWHM,
    dLh = unname(dLh),
    resels = unname(resels)
  )

  if (nonstationary) {
    rpv_scale <- (4.0 * log(2.0))^(if (usez) 1.5 else 1.0)
    dlh_ns <- sum(cpp_res$dLh_map) / cpp_res$volume
    out$dLh_nonstationary <- unname(dlh_ns)
    out$resels_nonstationary <- unname(rpv_scale / dlh_ns)
    dlh_img <- wrap_volume_output(cpp_res$dLh_map, mask_meta)
    rpv_img <- wrap_volume_output(cpp_res$dLh_map / rpv_scale, mask_meta)
    vw_img  <- wrap_volume_output(cpp_res$voxel_weights, mask_meta)
    attr(dlh_img, "ptfce_smoothness_units") <- "dlh"
    attr(rpv_img, "ptfce_smoothness_units") <- "rpv"
    attr(vw_img,  "ptfce_smoothness_units") <- "weights"
    out$dLh_map <- dlh_img
    out$rpv_map <- rpv_img
    out$voxel_weights <- vw_img
  }

  out
}

interpolate <- function(v) {
  knots <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
             25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 500)
  vals <- c(1.5423138, 1.3757105, 1.2842680, 1.2272151, 1.1885232, 1.1606988,
            1.1398000, 1.1235677, 1.1106196, 1.1000651, 1.0913060, 1.0839261,
            1.0776276, 1.0721920, 1.0674553, 1.0632924, 1.0483053, 1.0390117,
            1.0281339, 1.0219834, 1.0180339, 1.0152850, 1.0132621, 1.0117115,
            1.0104851, 1.0068808, 1.0051200, 1.0033865, 1.0020191)

  if (v < 5) {
    retval <- vals[1L]
  } else if (v > 500) {
    retval <- 1.0321 / v + 1.0
  } else {
    retval <- stats::approx(x = knots, y = vals, xout = v)$y
  }
  retval^0.5
}

#' Convert Z-score value(s) to GRF FWER-corrected p-value(s)
#'
#' Uses the 2D or 3D Gaussian Random Field Expected Euler Characteristic approximation.
#' Vectorized over \code{Z} and \code{resel_count}.
#'
#' @param resel_count Number of resels in the search volume (\code{volume / resels}).
#' @param Z Numeric vector of Z-score values.
#' @param D Spatial dimension (\code{3} for 3D volumes, default; \code{2} for 2D planar images).
#'
#' @return Numeric vector of FWER-corrected p-values in \code{[0, 1]}.
#' @export
fwe.z2p <- function(resel_count, Z, D = 3L) {
  dim_eff <- if (as.integer(D)[1L] == 2L) 2L else 3L
  df <- data.frame(r = as.numeric(resel_count), z = as.numeric(Z))
  out <- rep(1.0, nrow(df))
  na_idx <- is.na(df$r) | is.na(df$z)
  out[na_idx] <- NA_real_

  if (dim_eff == 2L) {
    ec_const_2d <- (4.0 * log(2.0)) / ((2.0 * pi)^1.5)
    ok <- !na_idx & (df$z >= 1.0)
    if (any(ok)) {
      z_ok <- df$z[ok]
      log_p <- log(df$r[ok]) + log(ec_const_2d) - 0.5 * z_ok * z_ok + log(z_ok)
      out[ok] <- pmin(1.0, exp(log_p))
    }
    return(out)
  }

  ec_const <- 0.11694 # (4 * ln(2))^1.5 / (2 * pi)^2 (FSL convention)
  ok <- !na_idx & (df$z >= 2.0)
  if (any(ok)) {
    z_ok <- df$z[ok]
    log_p <- log(df$r[ok]) + log(ec_const) - 0.5 * z_ok * z_ok + log(z_ok * z_ok - 1.0)
    out[ok] <- pmin(1.0, exp(log_p))
  }
  out
}

#' Compute GRF FWER-corrected Z-score threshold(s) from p-value(s)
#'
#' Uses Newton-Raphson in log-space to invert the 2D or 3D Expected Euler Characteristic
#' FWER approximation accurately across arbitrary \code{FWEP} and \code{resel_count}.
#' Vectorized over \code{FWEP} and \code{resel_count}.
#'
#' @param resel_count Number of resels in the search volume (\code{volume / resels}).
#' @param FWEP Numeric vector of target FWER-corrected p-values (default: \code{0.05}).
#' @param D Spatial dimension (\code{3} for 3D volumes, default; \code{2} for 2D planar images).
#'
#' @return Numeric vector of Z-score thresholds.
#' @export
fwe.p2z <- function(resel_count, FWEP = 0.05, D = 3L) {
  dim_eff <- if (as.integer(D)[1L] == 2L) 2L else 3L
  ec_const <- if (dim_eff == 2L) (4.0 * log(2.0)) / ((2.0 * pi)^1.5) else 0.11694
  df <- data.frame(r = as.numeric(resel_count), fwep = as.numeric(FWEP))
  out <- numeric(nrow(df))

  for (i in seq_len(nrow(df))) {
    r <- df$r[i]
    fwep <- df$fwep[i]
    if (is.na(r) || is.na(fwep) || r <= 0 || fwep <= 0) {
      out[i] <- NA_real_
      next
    }
    log_target <- log(fwep) - log(r) - log(ec_const)
    if (log_target > log(0.5)) {
      out[i] <- 0.0
      next
    }

    if (dim_eff == 2L) {
      # 2D: g(z) = log(z) - 0.5 * z^2 - log_target; maximum at z = 1 is -0.5
      if (log_target >= -0.5) {
        out[i] <- 1.0
        next
      }
      z <- max(1.1, sqrt(-2.0 * log_target + log(max(1.0, -2.0 * log_target))))
      for (it in 1:30) {
        g <- log(z) - 0.5 * z * z - log_target
        gp <- (1.0 / z) - z
        step <- g / gp
        z_new <- max(1.0000001, z - step)
        if (abs(z_new - z) <= 1e-13 * z_new) {
          z <- z_new
          break
        }
        z <- z_new
      }
      out[i] <- z
    } else {
      # At z = 2, log(z^2 - 1) - 0.5*z^2 = log(3) - 2 = -0.9013877; FSL returns 0 for p > 0.5
      if (log_target >= log(3.0) - 2.0) {
        out[i] <- 2.0
        next
      }

      # Solve g(z) = log(z^2 - 1) - 0.5 * z^2 - log_target = 0 for z > 2
      # g'(z) = 2*z / (z^2 - 1) - z = -z * (z^2 - 3) / (z^2 - 1)
      # Note g(z) has maximum at z = sqrt(3) < 2 and is strictly decreasing & concave for z > 2!
      z <- max(2.1, sqrt(-2.0 * log_target + 2.0 * log(max(2.0, -2.0 * log_target))))
      for (it in 1:30) {
        z2 <- z * z
        g <- log(z2 - 1.0) - 0.5 * z2 - log_target
        gp <- (2.0 * z) / (z2 - 1.0) - z
        step <- g / gp
        z_new <- max(2.0000001, z - step)
        if (abs(z_new - z) <= 1e-13 * z_new) {
          z <- z_new
          break
        }
        z <- z_new
      }
      out[i] <- z
    }
  }
  out
}
