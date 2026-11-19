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
#' @param fwe_p Numeric vector of Family-Wise Error Rate alpha levels for which GRF
#'   Z-score thresholds are computed (default: \code{0.05}). In \code{sided = "two.sided"}
#'   mode, each tail is controlled at \code{fwe_p / 2} so \code{abs(Z) >= fwe.Z} controls
#'   two-sided FWER at \code{fwe_p}.
#' @param legacy Logical; if \code{TRUE}, uses upstream \code{Nh}-point discrete summation
#'   from \code{logpmin} to \code{logpmax} instead of analytical sub-\code{ZestThr}
#'   integration with fixed \code{dh} (default: \code{FALSE}).
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
#' @return An object of class \code{"ptfce"} (a list) containing:
#' \itemize{
#'   \item \strong{p}: Uncorrected pTFCE-enhanced p-values (same image class as \code{img}).
#'   \item \strong{logp}: Negative natural logarithm \code{-log(p)} of enhanced p-values.
#'   \item \strong{Z}: Enhanced Z-score image (signed when \code{sided = "two.sided"}).
#'   \item \strong{number_of_resels}: Estimated number of resels in the mask.
#'   \item \strong{fwer0.05.Z}: GRF FWER Z-score threshold at \code{alpha = 0.05}
#'     (adjusted to \code{0.025} per tail when \code{sided = "two.sided"}).
#'   \item \strong{fwe.Z}: Named numeric vector of GRF FWER Z-score thresholds corresponding
#'     to \code{fwe_p}.
#'   \item \strong{sided}: Enhancement sidedness used.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard bi-sided fMRI enhancement using 4D residuals for smoothness:
#' res <- ptfce("zstat1.nii.gz", mask = "mask.nii.gz",
#'              residual = "res4d.nii.gz", dof = 42, sided = "two.sided")
#' sig_map <- res$Z
#' sig_map[abs(sig_map) < res$fwer0.05.Z] <- 0
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
                  legacy = FALSE) {
  sided <- match.arg(sided)
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

  autosmooth <- FALSE
  if (is_missing_arg(Rd) || is_missing_arg(V)) {
    autosmooth <- TRUE
    mask_clean <- array(as.numeric(inside), dim = d3)
    if (!is_missing_arg(residual)) {
      if (!is_missing_arg(dof)) {
        smooth <- smoothest(residual, mask_clean, dof = dof, verbose = verbose)
      } else {
        warning("Residual supplied without `dof`: defaulting back to Z-score map smoothness estimation.")
        smooth <- smoothest(img_arr, mask_clean, verbose = verbose)
      }
    } else {
      smooth <- smoothest(img_arr, mask_clean, verbose = verbose)
    }
    V <- smooth$volume
    Rd <- smooth$dLh * V
    resels <- smooth$resels
  } else {
    V <- as.numeric(V)[1L]
    Rd <- as.numeric(Rd)[1L]
    if (!is.finite(V) || V <= 0 || !is.finite(Rd) || Rd <= 0) {
      stop("`V` and `Rd` must be positive finite numbers.", call. = FALSE)
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
    cat(sprintf("* Performing pTFCE (%s)...\n", sided))
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
      analytical_subthr = analytical_subthr
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

  # Compute GRF FWER thresholds
  if (autosmooth || !is_missing_arg(resels)) {
    resels_val <- as.numeric(resels)[1L]
    number_of_resels <- V / resels_val
    alpha_scale <- if (sided == "two.sided") 0.5 else 1.0
    fwer0.05.Z <- fwe.p2z(number_of_resels, 0.05 * alpha_scale)
    fwe_Z <- stats::setNames(
      fwe.p2z(number_of_resels, fwe_p * alpha_scale),
      format(fwe_p, trim = TRUE)
    )
  } else {
    warning("For GRF-based FWER correction, please specify `resels` or allow automatic smoothness estimation.")
    number_of_resels <- NA_real_
    fwer0.05.Z <- NA_real_
    fwe_Z <- stats::setNames(rep(NA_real_, length(fwe_p)), format(fwe_p, trim = TRUE))
  }

  out <- list(
    p = wrap_volume_output(p_arr, img_meta),
    logp = wrap_volume_output(logp_arr, img_meta),
    Z = wrap_volume_output(Z_arr, img_meta),
    number_of_resels = number_of_resels,
    fwer0.05.Z = fwer0.05.Z,
    fwe.Z = fwe_Z,
    sided = sided
  )
  class(out) <- "ptfce"
  out
}

#' @export
print.ptfce <- function(x, ...) {
  cat("pTFCE result (", x$sided, ")\n", sep = "")
  cat("  Number of resels : ", format(x$number_of_resels, digits = 5), "\n", sep = "")
  cat("  FWER 0.05 Z-thr  : ", format(x$fwer0.05.Z, digits = 5), "\n", sep = "")
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
#'
#' @return Expected value of cluster size in voxels.
#' @export
Es <- function(h, V, Rd) {
  es_cpp(as.numeric(h), as.numeric(V)[1L], as.numeric(Rd)[1L])
}

dvox <- function(h) {
  stats::dnorm(h)
}

pvox <- function(h) {
  stats::pnorm(h, lower.tail = FALSE)
}

dclust <- function(h, V, Rd, c, ZestThr = 1.3, CONST = 1) {
  if (length(c) == 1L && is.na(c)) {
    return(dvox(h))
  }
  h_num <- as.numeric(h)
  lam <- (Es(h_num, V, Rd) / gamma(2.5))^(-2 / 3)
  log_d <- log(lam) - lam * (c^(2 / 3)) + log(CONST)
  out <- exp(log_d)
  out[h_num < ZestThr] <- 0
  out
}

dvox.clust <- function(h, V, Rd, c, ZestThr = 1.3) {
  if (length(c) == 1L && is.na(c)) {
    return(dvox(h))
  }
  dvox_clust_cpp(as.numeric(h), as.numeric(V)[1L], as.numeric(Rd)[1L], as.numeric(c)[1L], ZestThr)
}

pvox.clust <- function(V, Rd, c, actH, ZestThr = 1.3) {
  nlp <- neglog_pvox_clust_cpp(
    as.numeric(V)[1L],
    as.numeric(Rd)[1L],
    as.numeric(c),
    as.numeric(actH),
    as.numeric(ZestThr)[1L]
  )
  p <- exp(-nlp)
  p[p == 0] <- .Machine$double.xmin
  p
}

#' Estimate global image smoothness
#'
#' Fast C++ port of FSL's \code{smoothest} tool (Jenkinson, FMRIB Tech Report TR00DF1)
#' supporting both 2D/3D Z-score maps and 4D GLM residual time-series images.
#'
#' @param img A 2D/3D Z-score image or 4D residual image (\code{RNifti}, \code{oro.nifti},
#'   \code{array}, or file path).
#' @param mask Spatial mask image (\code{RNifti}, \code{oro.nifti}, \code{array}, or file path).
#' @param dof Residual degrees of freedom (required when \code{img} is a 4D residual image).
#' @param verbose Logical; print diagnostic messages if \code{TRUE} (default).
#'
#' @return A list containing:
#' \itemize{
#'   \item \strong{volume}: Mask volume in voxels.
#'   \item \strong{sigmasq}: Estimated Gaussian smoothness variance in x, y, z (voxels^2).
#'   \item \strong{FWHM}: Full-width at half-maximum in x, y, z (voxels).
#'   \item \strong{dLh}: Determinant of \eqn{\Lambda^{1/2}} (\code{voxels^-3}).
#'   \item \strong{resels}: Resel size in voxels per resel.
#' }
#' @export
smoothest <- function(img, mask, dof = NA, verbose = TRUE) {
  if (verbose) cat("* Estimating smoothness based on the data...\n")

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

  cpp_res <- smoothest_cpp(
    img = img_meta$data,
    mask = mask_meta$data,
    nx = d3[1L],
    ny = d3[2L],
    nz = d3[3L],
    nt = nt
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
    dLh <- dLh * interpolate(as.numeric(dof)[1L])
  }

  FWHM <- c(
    x = sqrt(8.0 * log(2.0) * sigmasq[1L]),
    y = sqrt(8.0 * log(2.0) * sigmasq[2L]),
    z = if (usez) sqrt(8.0 * log(2.0) * sigmasq[3L]) else 0.0
  )

  resels <- FWHM[1L] * FWHM[2L] * (if (usez) FWHM[3L] else 1.0)

  list(
    volume = cpp_res$volume,
    sigmasq = sigmasq,
    FWHM = FWHM,
    dLh = unname(dLh),
    resels = unname(resels)
  )
}

interpolate <- function(v) {
  knots <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
             25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 500)
  vals <- c(1.5423138, 1.3757105, 1.2842680, 1.2272151, 1.1885232, 1.1606988,
            1.1398000, 1.1235677, 1.1106196, 1.1000651, 1.0913060, 1.0839261,
            1.0776276, 1.0721920, 1.0674553, 1.0632924, 1.0483053, 1.0390117,
            1.0281339, 1.0219834, 1.0180339, 1.0152850, 1.0132621, 1.0117115,
            1.0104851, 1.0068808, 1.0051200, 1.0033865, 1.0020191)

  if (v < 5) return(1.1)
  if (v > 500) {
    retval <- 1.0321 / v + 1.0
  } else {
    retval <- stats::approx(x = knots, y = vals, xout = v)$y
  }
  retval^0.5
}

#' Convert Z-score value(s) to GRF FWER-corrected p-value(s)
#'
#' Uses the 3D Gaussian Random Field Expected Euler Characteristic approximation.
#' Vectorized over \code{Z} and \code{resel_count}.
#'
#' @param resel_count Number of resels in the search volume (\code{volume / resels}).
#' @param Z Numeric vector of Z-score values.
#'
#' @return Numeric vector of FWER-corrected p-values in \code{[0, 1]}.
#' @export
fwe.z2p <- function(resel_count, Z) {
  ec_const <- 0.11694 # (4 * ln(2))^1.5 / (2 * pi)^2 (FSL convention)
  df <- data.frame(r = as.numeric(resel_count), z = as.numeric(Z))
  out <- rep(1.0, nrow(df))
  na_idx <- is.na(df$r) | is.na(df$z)
  out[na_idx] <- NA_real_

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
#' Uses Newton-Raphson in log-space to invert the 3D Expected Euler Characteristic
#' FWER approximation accurately across arbitrary \code{FWEP} and \code{resel_count}.
#' Vectorized over \code{FWEP} and \code{resel_count}.
#'
#' @param resel_count Number of resels in the search volume (\code{volume / resels}).
#' @param FWEP Numeric vector of target FWER-corrected p-values (default: \code{0.05}).
#'
#' @return Numeric vector of Z-score thresholds.
#' @export
fwe.p2z <- function(resel_count, FWEP = 0.05) {
  ec_const <- 0.11694 # (4 * ln(2))^1.5 / (2 * pi)^2 (FSL convention)
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
    # At z = 2, log(z^2 - 1) - 0.5*z^2 = log(3) - 2 = -0.9013877; FSL returns 0 for p > 0.5
    if (log_target > log(0.5)) {
      out[i] <- 0.0
      next
    }
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
  out
}
