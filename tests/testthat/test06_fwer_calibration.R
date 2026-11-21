# Tests for empirical FWER threshold calibration ("grf_sim" and "sign_flip") in v0.4.1

smooth_white_3d_cal <- function(x, sigma) {
  if (sigma <= 0) return(x)
  rad <- ceiling(3.5 * sigma)
  k <- exp(-0.5 * ((-rad:rad)^2) / (sigma^2))
  k <- k / sqrt(sum(k^2))
  d <- dim(x)
  pd <- d + 2L * rad
  pad <- array(stats::rnorm(prod(pd)), dim = pd)
  pad[(rad + 1L):(rad + d[1L]),
      (rad + 1L):(rad + d[2L]),
      (rad + 1L):(rad + d[3L])] <- x
  y1 <- array(0, dim = c(d[1L], pd[2L], pd[3L]))
  for (i in seq_along(k)) y1 <- y1 + k[i] * pad[i:(i + d[1L] - 1L), , ]
  y2 <- array(0, dim = c(d[1L], d[2L], pd[3L]))
  for (i in seq_along(k)) y2 <- y2 + k[i] * y1[, i:(i + d[2L] - 1L), ]
  y3 <- array(0, dim = d)
  for (i in seq_along(k)) y3 <- y3 + k[i] * y2[, , i:(i + d[3L] - 1L)]
  y3
}

test_that("Default fwe_method = 'grf' is backward-compatible and grf_sim is reproducible", {
  z_path <- system.file("extdata", "test02A_in.nii.gz", package = "pTFCE")
  m_path <- system.file("extdata", "test02_mask.nii.gz", package = "pTFCE")
  z_nii <- RNifti::readNifti(z_path)
  m_nii <- RNifti::readNifti(m_path)

  res_grf <- ptfce(z_nii, m_nii, sided = "two.sided", fwe_p = c(0.05, 0.01), verbose = FALSE)
  expect_identical(res_grf$fwe_method, "grf")
  expect_null(res_grf$null_max_Z)
  expect_equal(res_grf$fwer0.05.Z, fwe.p2z(res_grf$number_of_resels, 0.025), tolerance = 1e-12)

  set.seed(101)
  res_sim1 <- ptfce(z_nii, m_nii, sided = "two.sided", fwe_p = c(0.05, 0.01),
                    fwe_method = "grf_sim", n_perm = 99L, verbose = FALSE)
  set.seed(101)
  res_sim2 <- ptfce(z_nii, m_nii, sided = "two.sided", fwe_p = c(0.05, 0.01),
                    fwe_method = "grf_sim", n_perm = 99L, verbose = FALSE)

  expect_identical(res_sim1$fwe_method, "grf_sim")
  expect_length(res_sim1$null_max_Z, 99L)
  expect_equal(res_sim1$null_max_Z, res_sim2$null_max_Z, tolerance = 1e-14)
  expect_equal(res_sim1$fwer0.05.Z, res_sim2$fwer0.05.Z, tolerance = 1e-14)
  expect_equal(as.array(res_sim1$Z), as.array(res_grf$Z), tolerance = 1e-14)

  # Empirical pTFCE threshold relieves raw-field GRF conservativeness while staying valid
  expect_lt(res_sim1$fwer0.05.Z, res_grf$fwer0.05.Z)
  expect_gt(res_sim1$fwer0.05.Z, 2.8)
  expect_gt(res_sim1$fwe.Z[["0.01"]], res_sim1$fwe.Z[["0.05"]])
})

test_that("grf_sim and sign_flip (4D copes & 4D residuals with covariates) agree under Gaussian data", {
  set.seed(2027)
  d <- c(24L, 20L, 14L)
  n_subj <- 24L
  mask <- array(0, dim = d)
  mask[3:22, 3:18, 3:12] <- 1

  # Simulate non-stationary 4D subject COPEs: left rough (sigma = 0.8), right smooth (sigma = 1.9)
  copes4d <- array(0, dim = c(d, n_subj))
  for (s in seq_len(n_subj)) {
    w <- array(stats::rnorm(prod(d)), dim = d)
    v_r <- smooth_white_3d_cal(w, sigma = 0.80)
    v_s <- smooth_white_3d_cal(w, sigma = 1.90)
    vol <- v_r
    vol[13:24, , ] <- v_s[13:24, , ]
    copes4d[, , , s] <- vol
  }
  # Inject true focal activation in rough hemisphere
  copes4d[6:8, 9:11, 6:8, ] <- copes4d[6:8, 9:11, 6:8, ] + 0.85

  # Group one-sample t -> Z map and GLM residuals with 2 nuisance covariates (dof = n_subj - 3)
  X <- cbind(1, scale(seq_len(n_subj)), scale(sin(seq_len(n_subj))))
  P_perp <- diag(n_subj) - X %*% solve(crossprod(X), t(X))
  dof_glm <- n_subj - ncol(X)

  mat_y <- matrix(copes4d, nrow = prod(d), ncol = n_subj)
  mat_res <- mat_y %*% P_perp
  res4d <- array(mat_res, dim = c(d, n_subj))

  m_bar <- rowMeans(mat_y)
  sd_y <- sqrt(rowSums((mat_y - m_bar)^2) / (n_subj - 1L))
  t_obs <- array(sqrt(n_subj) * m_bar / pmax(1e-12, sd_y), dim = d)
  z_obs <- t2z(t_obs, dof = n_subj - 1L)

  # 1. Analytical GRF (non-stationary)
  r_grf <- ptfce(z_obs, mask, residual = res4d, dof = dof_glm,
                 sided = "two.sided", nonstationary = TRUE,
                 fwe_method = "grf", verbose = FALSE)

  # 2. Parametric GRF simulation (non-stationary)
  set.seed(51)
  r_sim <- ptfce(z_obs, mask, residual = res4d, dof = dof_glm,
                 sided = "two.sided", nonstationary = TRUE,
                 fwe_method = "grf_sim", n_perm = 199L, verbose = FALSE)

  # 3. Rademacher sign-flipping on 4D GLM residuals (with 3 projected-out columns)
  set.seed(52)
  r_sf_res <- ptfce(z_obs, mask, residual = res4d, dof = dof_glm,
                    sided = "two.sided", nonstationary = TRUE,
                    fwe_method = "sign_flip", n_perm = 199L, verbose = FALSE)

  # 4. Rademacher sign-flipping on 4D COPEs (auto-selects sign_flip and estimates 4D smoothness)
  set.seed(53)
  r_sf_copes <- ptfce(z_obs, mask, copes = copes4d,
                      sided = "two.sided", nonstationary = TRUE,
                      n_perm = 199L, verbose = FALSE)

  expect_identical(r_sf_copes$fwe_method, "sign_flip")
  expect_true(all(is.finite(r_sim$null_max_Z)))
  expect_true(all(is.finite(r_sf_res$null_max_Z)))
  expect_true(all(is.finite(r_sf_copes$null_max_Z)))

  # All three calibrated methods relieve raw-field GRF conservativeness and agree closely
  expect_lt(r_sim$fwer0.05.Z, r_grf$fwer0.05.Z - 0.15)
  expect_lt(r_sf_res$fwer0.05.Z, r_grf$fwer0.05.Z - 0.15)
  expect_lt(r_sf_copes$fwer0.05.Z, r_grf$fwer0.05.Z - 0.15)
  expect_lt(abs(r_sim$fwer0.05.Z - r_sf_res$fwer0.05.Z), 0.35)
  expect_lt(abs(r_sim$fwer0.05.Z - r_sf_copes$fwer0.05.Z), 0.35)

  # Reproducibility of sign_flip with set.seed()
  set.seed(52)
  r_sf_res2 <- ptfce(z_obs, mask, residual = res4d, dof = dof_glm,
                     sided = "two.sided", nonstationary = TRUE,
                     fwe_method = "sign_flip", n_perm = 199L, verbose = FALSE)
  expect_equal(r_sf_res$null_max_Z, r_sf_res2$null_max_Z, tolerance = 1e-14)
})

test_that("Positive/negative tails, manual (V, Rd), legacy mode, and input guards work cleanly", {
  set.seed(99)
  d <- c(16L, 16L, 12L)
  mask <- array(1, dim = d)
  z_img <- smooth_white_3d_cal(array(stats::rnorm(prod(d)), dim = d), sigma = 1.2)

  # Manual (V, Rd) without resels still produces valid calibrated fwer0.05.Z via grf_sim
  r_man <- ptfce(z_img, mask, V = prod(d), Rd = 150, sided = "positive",
                 fwe_method = "grf_sim", n_perm = 49L, verbose = FALSE)
  expect_true(is.na(r_man$number_of_resels))
  expect_true(is.finite(r_man$fwer0.05.Z) && r_man$fwer0.05.Z > 2.5)

  # Negative tail + legacy mode
  r_neg_leg <- ptfce(z_img, mask, sided = "negative", legacy = TRUE, Nh = 40L,
                     fwe_method = "grf_sim", n_perm = 29L, verbose = FALSE)
  expect_true(is.finite(r_neg_leg$fwer0.05.Z))

  # Error guards for sign_flip without 4D inputs or too few volumes
  expect_error(
    ptfce(z_img, mask, fwe_method = "sign_flip", verbose = FALSE),
    "requires either a 4D `copes` image or a 4D `residual` image"
  )
  tiny_4d <- array(stats::rnorm(prod(d) * 4L), dim = c(d, 4L))
  expect_error(
    ptfce(z_img, mask, copes = tiny_4d, fwe_method = "sign_flip", verbose = FALSE),
    "at least 6 volumes"
  )

  # Warning when 4D residuals have strong temporal AR(1) autocorrelation
  ar_res <- array(0, dim = c(d, 20L))
  ar_res[, , , 1L] <- stats::rnorm(prod(d))
  for (t in 2:20) {
    ar_res[, , , t] <- 0.75 * ar_res[, , , t - 1L] + sqrt(1 - 0.75^2) * stats::rnorm(prod(d))
  }
  expect_warning(
    ptfce(z_img, mask, residual = ar_res, dof = 19L,
          fwe_method = "sign_flip", n_perm = 39L, verbose = FALSE),
    "Substantial lag-1 autocorrelation"
  )

  # Warning when n_perm is too small for requested extreme fwe_p
  expect_warning(
    ptfce(z_img, mask, fwe_method = "grf_sim", n_perm = 39L, fwe_p = 0.001, verbose = FALSE),
    "is too small for `fwe_p"
  )
})
