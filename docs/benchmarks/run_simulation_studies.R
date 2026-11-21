# Reproducible Monte Carlo Simulation Studies for HallquistLab/pTFCE (v0.4.0 - v0.4.1)
#
# Evaluates:
#   Study 1: Spatial smoothness non-stationarity (regional false-positive disparity
#            across rough subcortex vs. smooth DMN cortex, and subcortical power)
#   Study 2: FWER threshold calibration ("grf" vs. "grf_sim" vs. "sign_flip" on
#            4D GLM residuals and 4D subject COPEs) under:
#              (2A) Non-stationary Gaussian group fMRI global null
#              (2B) Bi-sided group activation & deactivation power + smooth-DMN leakage
#              (2C) Heavy-tailed / spatial-outlier non-Gaussian group fMRI null
#              (2D) Whole-brain real NIfTI execution speed (301,461 mask voxels)
#
# Usage:
#   Rscript docs/benchmarks/run_simulation_studies.R

library(pTFCE)
library(RNifti)

cat("pTFCE Simulation Suite — Package Version:", as.character(packageVersion("pTFCE")), "\n\n")

make_psf_1d <- function(fwhm) {
  sigma <- fwhm / sqrt(8 * log(2))
  rad <- ceiling(3.5 * sigma)
  k <- exp(-0.5 * ((-rad:rad)^2) / (sigma^2))
  list(k_noise = k / sqrt(sum(k^2)), k_sig = k / sum(k), rad = rad)
}

conv3d <- function(pad, k, d) {
  rad <- (length(k) - 1L) %/% 2L
  pd <- dim(pad)
  off <- (pd[1L] - d[1L]) %/% 2L - rad
  sub <- pad[(off + 1L):(off + d[1L] + 2L * rad),
             (off + 1L):(off + d[2L] + 2L * rad),
             (off + 1L):(off + d[3L] + 2L * rad)]
  sp <- dim(sub)
  y1 <- array(0, dim = c(d[1L], sp[2L], sp[3L]))
  for (i in seq_along(k)) y1 <- y1 + k[i] * sub[i:(i + d[1L] - 1L), , ]
  y2 <- array(0, dim = c(d[1L], d[2L], sp[3L]))
  for (i in seq_along(k)) y2 <- y2 + k[i] * y1[, i:(i + d[2L] - 1L), ]
  y3 <- array(0, dim = d)
  for (i in seq_along(k)) y3 <- y3 + k[i] * y2[, , i:(i + d[3L] - 1L)]
  y3
}

# 3-zone anatomical smoothness layout:
#   Zone 1 (Subcortex / brainstem):  FWHM = 3.0 voxels (1 resel =  27 voxels)
#   Zone 2 (Lateral association):    FWHM = 4.8 voxels (1 resel = 111 voxels)
#   Zone 3 (Medial DMN / precuneus): FWHM = 7.2 voxels (1 resel = 373 voxels)
d <- c(36L, 26L, 18L)
mask <- array(1, dim = d)
psf1 <- make_psf_1d(3.0)
psf2 <- make_psf_1d(4.8)
psf3 <- make_psf_1d(7.2)
max_rad <- psf3$rad

x_coord <- 1:d[1L]
s_blend <- function(x, c0) 1 / (1 + exp(-(x - c0) / 1.3))
a1 <- 1 - s_blend(x_coord, 12.5)
a3 <- s_blend(x_coord, 24.5)
a2 <- pmax(0, 1 - a1 - a3)

set.seed(101)
cal_ssq <- numeric(d[1L])
for (b in 1:150) {
  pw <- array(rnorm(prod(d + 2L * max_rad)), dim = d + 2L * max_rad)
  comb <- sweep(conv3d(pw, psf1$k_noise, d), 1, a1, `*`) +
          sweep(conv3d(pw, psf2$k_noise, d), 1, a2, `*`) +
          sweep(conv3d(pw, psf3$k_noise, d), 1, a3, `*`)
  cal_ssq <- cal_ssq + apply(comb^2, 1, mean)
}
x_sd <- sqrt(cal_ssq / 150)

sim_subject_field <- function(df_innov = Inf) {
  n_tot <- prod(d + 2L * max_rad)
  w <- if (is.finite(df_innov)) rt(n_tot, df = df_innov) * sqrt((df_innov - 2) / df_innov) else rnorm(n_tot)
  pw <- array(w, dim = d + 2L * max_rad)
  comb <- sweep(conv3d(pw, psf1$k_noise, d), 1, a1, `*`) +
          sweep(conv3d(pw, psf2$k_noise, d), 1, a2, `*`) +
          sweep(conv3d(pw, psf3$k_noise, d), 1, a3, `*`)
  sweep(comb, 1, x_sd, `/`)
}

z1_mask <- array(FALSE, dim = d); z1_mask[1:11, , ]  <- TRUE
z2_mask <- array(FALSE, dim = d); z2_mask[14:23, , ] <- TRUE
z3_mask <- array(FALSE, dim = d); z3_mask[26:36, , ] <- TRUE

group_stats_from_copes <- function(copes4d) {
  ns <- dim(copes4d)[4L]
  mat <- matrix(copes4d, nrow = prod(d), ncol = ns)
  m <- rowMeans(mat)
  res_mat <- mat - m
  sd_v <- sqrt(rowSums(res_mat^2) / (ns - 1L))
  t_arr <- array(sqrt(ns) * m / pmax(1e-12, sd_v), dim = d)
  list(z = t2z(t_arr, dof = ns - 1L), res4d = array(res_mat, dim = c(d, ns)), dof = ns - 1L)
}

# --- STUDY 1 & 2A: Global Null Calibration & Regional False-Positive Parity ---
set.seed(2026041)
R_null <- 200L; N_subj <- 25L; B_perm <- 199L
methods <- c("grf_stat", "grf_sim_stat", "grf_ns", "grf_sim_ns", "sf_res_ns", "sf_cope_ns")
fwer_hit <- matrix(FALSE, R_null, length(methods), dimnames = list(NULL, methods))
thr_mat  <- matrix(0.0,   R_null, length(methods), dimnames = list(NULL, methods))
p001_st  <- matrix(0.0,   R_null, 3L)
p001_ns  <- matrix(0.0,   R_null, 3L)

for (r in seq_len(R_null)) {
  copes <- array(0, dim = c(d, N_subj))
  for (s in seq_len(N_subj)) copes[, , , s] <- sim_subject_field()
  gs <- group_stats_from_copes(copes)
  sm <- smoothest(gs$res4d, mask, dof = gs$dof, nonstationary = TRUE, smooth_fwhm = 3.0, verbose = FALSE)

  o_st_grf <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = FALSE, fwe_method = "grf", verbose = FALSE)
  o_st_sim <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = FALSE, fwe_method = "grf_sim", n_perm = B_perm, verbose = FALSE)
  o_ns_grf <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE,  fwe_method = "grf", verbose = FALSE)
  o_ns_sim <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE,  fwe_method = "grf_sim", n_perm = B_perm, verbose = FALSE)
  o_ns_sfr <- ptfce(gs$z, mask, residual = gs$res4d, dof = gs$dof, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE, fwe_method = "sign_flip", n_perm = B_perm, verbose = FALSE)
  o_ns_sfc <- ptfce(gs$z, mask, copes = copes, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE, fwe_method = "sign_flip", n_perm = B_perm, verbose = FALSE)

  p001_st[r, ] <- c(mean(o_st_grf$p[z1_mask] < 0.001), mean(o_st_grf$p[z2_mask] < 0.001), mean(o_st_grf$p[z3_mask] < 0.001))
  p001_ns[r, ] <- c(mean(o_ns_grf$p[z1_mask] < 0.001), mean(o_ns_grf$p[z2_mask] < 0.001), mean(o_ns_grf$p[z3_mask] < 0.001))

  objs <- list(o_st_grf, o_st_sim, o_ns_grf, o_ns_sim, o_ns_sfr, o_ns_sfc)
  for (m_i in seq_along(objs)) {
    th <- objs[[m_i]]$fwer0.05.Z
    thr_mat[r, m_i]  <- th
    fwer_hit[r, m_i] <- any(abs(objs[[m_i]]$Z) >= th)
  }
}

cat("=== STUDY 1: Regional Uncorrected p < 0.001 Rate (x10^3) [Rough / Mod / Smooth] ===\n")
cat(sprintf("Stationary pTFCE:     %.3f / %.3f / %.3f (Smooth/Rough ratio = %.1fx)\n",
            1e3 * mean(p001_st[, 1]), 1e3 * mean(p001_st[, 2]), 1e3 * mean(p001_st[, 3]),
            mean(p001_st[, 3]) / max(1e-8, mean(p001_st[, 1]))))
cat(sprintf("Non-Stationary pTFCE: %.3f / %.3f / %.3f (Smooth/Rough ratio = %.1fx)\n\n",
            1e3 * mean(p001_ns[, 1]), 1e3 * mean(p001_ns[, 2]), 1e3 * mean(p001_ns[, 3]),
            mean(p001_ns[, 3]) / max(1e-8, mean(p001_ns[, 1]))))

cat("=== STUDY 2A: Whole-Brain FWER Under Non-Stationary Gaussian Null (Nominal = 0.050) ===\n")
print(data.frame(Method = methods, Mean_Thr_Z = round(colMeans(thr_mat), 3),
                 WholeBrain_FWER = round(colMeans(fwer_hit), 3)), row.names = FALSE)

# --- STUDY 2B: Bi-Sided Power & Smooth-DMN Leakage ---
pad_zero <- function() array(0, dim = d + 2L * max_rad); off <- max_rad
p_sub <- pad_zero(); p_sub[(off + 5):(off + 7), (off + 12):(off + 14), (off + 8):(off + 10)] <- 1
s_sub <- conv3d(p_sub, psf1$k_sig, d); s_sub <- (4.35 / sqrt(N_subj)) * (s_sub / max(s_sub))
p_lat <- pad_zero(); p_lat[(off + 16):(off + 20), (off + 11):(off + 15), (off + 7):(off + 11)] <- 1
s_lat <- conv3d(p_lat, psf2$k_sig, d); s_lat <- (3.75 / sqrt(N_subj)) * (s_lat / max(s_lat))
p_dmn <- pad_zero(); p_dmn[(off + 28):(off + 34), (off + 10):(off + 16), (off + 6):(off + 12)] <- 1
s_dmn <- conv3d(p_dmn, psf3$k_sig, d); s_dmn <- (-3.75 / sqrt(N_subj)) * (s_dmn / max(s_dmn))

sig_subj <- s_sub + s_lat + s_dmn
m_sub <- s_sub >= 0.5 * max(s_sub); m_lat <- s_lat >= 0.5 * max(s_lat); m_dmn <- s_dmn <= 0.5 * min(s_dmn)
m_null_dmn <- z3_mask & (abs(sig_subj * sqrt(N_subj)) < 0.15)

R_sig <- 150L
sig_methods <- c("RawZ_GRF", "Stat_GRF", "Stat_GRFsim", "NS_GRF", "NS_GRFsim", "NS_SignFlipRes", "NS_SignFlipCope")
M_sig <- length(sig_methods)
pow_sub <- tpr_sub <- pow_lat <- tpr_lat <- pow_dmn <- tpr_dmn <- dmn_fp <-
  matrix(0.0, R_sig, M_sig, dimnames = list(NULL, sig_methods))

set.seed(2026044)
for (r in seq_len(R_sig)) {
  copes <- array(0, dim = c(d, N_subj))
  for (s in seq_len(N_subj)) copes[, , , s] <- sim_subject_field() + sig_subj
  gs <- group_stats_from_copes(copes)
  sm <- smoothest(gs$res4d, mask, dof = gs$dof, nonstationary = TRUE, smooth_fwhm = 3.0, verbose = FALSE)

  o_st_grf <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = FALSE, fwe_method = "grf", verbose = FALSE)
  o_st_sim <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = FALSE, fwe_method = "grf_sim", n_perm = B_perm, verbose = FALSE)
  o_ns_grf <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE,  fwe_method = "grf", verbose = FALSE)
  o_ns_sim <- ptfce(gs$z, mask, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE,  fwe_method = "grf_sim", n_perm = B_perm, verbose = FALSE)
  o_ns_sfr <- ptfce(gs$z, mask, residual = gs$res4d, dof = gs$dof, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE, fwe_method = "sign_flip", n_perm = B_perm, verbose = FALSE)
  o_ns_sfc <- ptfce(gs$z, mask, copes = copes, local_smoothness = sm, sided = "two.sided", nonstationary = TRUE, fwe_method = "sign_flip", n_perm = B_perm, verbose = FALSE)

  maps <- list(gs$z, o_st_grf$Z, o_st_sim$Z, o_ns_grf$Z, o_ns_sim$Z, o_ns_sfr$Z, o_ns_sfc$Z)
  thrs <- c(o_st_grf$fwer0.05.Z, o_st_grf$fwer0.05.Z, o_st_sim$fwer0.05.Z,
            o_ns_grf$fwer0.05.Z, o_ns_sim$fwer0.05.Z, o_ns_sfr$fwer0.05.Z, o_ns_sfc$fwer0.05.Z)
  for (m_i in seq_len(M_sig)) {
    zm <- maps[[m_i]]; th <- thrs[m_i]
    pow_sub[r, m_i] <- any(zm[m_sub] >= th);  tpr_sub[r, m_i] <- mean(zm[m_sub] >= th)
    pow_lat[r, m_i] <- any(zm[m_lat] >= th);  tpr_lat[r, m_i] <- mean(zm[m_lat] >= th)
    pow_dmn[r, m_i] <- any(zm[m_dmn] <= -th); tpr_dmn[r, m_i] <- mean(zm[m_dmn] <= -th)
    dmn_fp[r, m_i]  <- any(abs(zm[m_null_dmn]) >= th)
  }
}

cat("\n=== STUDY 2B: Bi-Sided Power & Smooth-DMN Leakage (FWER = 0.05) ===\n")
print(data.frame(
  Method = sig_methods,
  Sub_Pow = round(100 * colMeans(pow_sub), 1), Sub_TPR = round(100 * colMeans(tpr_sub), 1),
  Lat_Pow = round(100 * colMeans(pow_lat), 1), Lat_TPR = round(100 * colMeans(tpr_lat), 1),
  DMN_Pow = round(100 * colMeans(pow_dmn), 1), DMN_TPR = round(100 * colMeans(tpr_dmn), 1),
  NullDMN_FP = round(100 * colMeans(dmn_fp), 1)
), row.names = FALSE)
