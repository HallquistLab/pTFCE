# Tests for optional non-stationary smoothness adjustment in smoothest() and ptfce()

# Helper: separable 3D Gaussian filter preserving unit variance
smooth_white_3d <- function(x, sigma) {
  if (sigma <= 0) return(x)
  rad <- ceiling(3.5 * sigma)
  k <- exp(-0.5 * ((-rad:rad)^2) / (sigma^2))
  k <- k / sqrt(sum(k^2))
  d <- dim(x)

  # Pad with fresh N(0,1) draws so boundaries remain stationary unit-variance
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

test_that("Uniform local_smoothness weights match stationary ptfce() to high precision", {
  z_path <- system.file("extdata", "test02A_in.nii.gz", package = "pTFCE")
  m_path <- system.file("extdata", "test02_mask.nii.gz", package = "pTFCE")

  z_nii <- RNifti::readNifti(z_path)
  m_nii <- RNifti::readNifti(m_path)
  sm <- smoothest(z_nii, m_nii, verbose = FALSE)
  V <- sm$volume
  Rd <- sm$dLh * V

  uniform_w <- array(as.numeric(as.array(m_nii) > 0.5), dim = dim(m_nii))

  res_stat <- ptfce(z_nii, m_nii, V = V, Rd = Rd, resels = sm$resels,
                    sided = "two.sided", nonstationary = FALSE, verbose = FALSE)
  res_unif <- ptfce(z_nii, m_nii, V = V, Rd = Rd, resels = sm$resels,
                    sided = "two.sided", nonstationary = TRUE,
                    local_smoothness = uniform_w, verbose = FALSE)

  expect_false(res_stat$nonstationary)
  expect_true(res_unif$nonstationary)
  expect_equal(as.array(res_unif$Z), as.array(res_stat$Z), tolerance = 1e-5)
  expect_equal(as.array(res_unif$p), as.array(res_stat$p), tolerance = 1e-5)
})

test_that("Non-stationary 4D residual smoothness boosts rough-region clusters and tempers smooth-region clusters", {
  set.seed(2026)
  d <- c(32L, 28L, 18L)
  nt <- 36L
  dof <- 35L
  mask <- array(1, dim = d)

  # Construct 4D residuals: Left half (x = 1..16) rough (sigma = 0.75, FWHM ~ 1.77),
  # Right half (x = 17..32) smooth (sigma = 2.0, FWHM ~ 4.71)
  res4d <- array(0, dim = c(d, nt))
  for (t in seq_len(nt)) {
    w_raw <- array(stats::rnorm(prod(d)), dim = d)
    rough_vol <- smooth_white_3d(w_raw, sigma = 0.75)
    smooth_vol <- smooth_white_3d(w_raw, sigma = 2.00)
    vol <- rough_vol
    vol[17:32, , ] <- smooth_vol[17:32, , ]
    res4d[, , , t] <- vol
  }

  sm_stat <- smoothest(res4d, mask, dof = dof, nonstationary = FALSE, verbose = FALSE)
  sm_ns   <- smoothest(res4d, mask, dof = dof, nonstationary = TRUE, smooth_fwhm = 3.0, verbose = FALSE)

  # Global smoothness fields must be 100% unchanged when nonstationary = TRUE
  expect_equal(sm_ns$volume, sm_stat$volume)
  expect_equal(sm_ns$sigmasq, sm_stat$sigmasq, tolerance = 1e-12)
  expect_equal(sm_ns$FWHM, sm_stat$FWHM, tolerance = 1e-12)
  expect_equal(sm_ns$dLh, sm_stat$dLh, tolerance = 1e-12)
  expect_equal(sm_ns$resels, sm_stat$resels, tolerance = 1e-12)

  # Local voxel_weights must have mean 1.0 and clearly distinguish rough vs smooth halves
  vw <- sm_ns$voxel_weights
  expect_equal(mean(vw), 1.0, tolerance = 1e-10)
  w_rough  <- mean(vw[4:13, 6:22, 5:14])
  w_smooth <- mean(vw[20:29, 6:22, 5:14])
  expect_gt(w_rough, 1.3)
  expect_lt(w_smooth, 0.5)
  expect_gt(w_rough / w_smooth, 3.0)

  # Create test Z-map with identical positive & negative 3x3x3 clusters in rough vs smooth halves
  z_img <- array(0, dim = d)
  # Positive clusters (Z = 3.4, peak 3.7)
  z_img[7:9, 9:11, 8:10]   <- 3.4; z_img[8, 10, 9]  <- 3.7 # rough half
  z_img[23:25, 9:11, 8:10] <- 3.4; z_img[24, 10, 9] <- 3.7 # smooth half
  # Negative clusters (Z = -3.4, peak -3.7)
  z_img[7:9, 18:20, 8:10]   <- -3.4; z_img[8, 19, 9]  <- -3.7 # rough half
  z_img[23:25, 18:20, 8:10] <- -3.4; z_img[24, 19, 9] <- -3.7 # smooth half

  out_stat <- ptfce(z_img, mask, residual = res4d, dof = dof,
                    sided = "two.sided", nonstationary = FALSE, verbose = FALSE)
  out_ns   <- ptfce(z_img, mask, residual = res4d, dof = dof,
                    sided = "two.sided", nonstationary = TRUE, verbose = FALSE)

  # Stationary mode treats both hemispheres identically
  expect_equal(out_stat$Z[8, 10, 9], out_stat$Z[24, 10, 9], tolerance = 1e-12)
  expect_equal(out_stat$Z[8, 19, 9], out_stat$Z[24, 19, 9], tolerance = 1e-12)

  # Non-stationary mode boosts the rough-hemisphere cluster and tempers the smooth-hemisphere cluster
  expect_gt(out_ns$Z[8, 10, 9], out_stat$Z[8, 10, 9])
  expect_lt(out_ns$Z[24, 10, 9], out_stat$Z[24, 10, 9])
  expect_gt(out_ns$Z[8, 10, 9] - out_ns$Z[24, 10, 9], 0.4)

  # Negative clusters exhibit exact mirror behavior and contrast anti-symmetry
  expect_lt(out_ns$Z[8, 19, 9], out_stat$Z[8, 19, 9])   # more negative in rough half
  expect_gt(out_ns$Z[24, 19, 9], out_stat$Z[24, 19, 9]) # tempered toward 0 in smooth half

  out_ns_flip <- ptfce(-z_img, mask, residual = res4d, dof = dof,
                       sided = "two.sided", nonstationary = TRUE, verbose = FALSE)
  expect_equal(out_ns$Z, -out_ns_flip$Z, tolerance = 1e-12)

  # Precomputed local_smoothness list or RPV image gives identical results and auto-enables nonstationary
  out_pre_list <- ptfce(z_img, mask, local_smoothness = sm_ns,
                        sided = "two.sided", verbose = FALSE)
  expect_true(out_pre_list$nonstationary)
  expect_equal(out_pre_list$Z, out_ns$Z, tolerance = 1e-12)

  out_pre_rpv <- ptfce(z_img, mask, V = sm_ns$volume,
                       Rd = sm_ns$dLh_nonstationary * sm_ns$volume,
                       resels = sm_ns$resels_nonstationary,
                       local_smoothness = sm_ns$rpv_map,
                       sided = "two.sided", verbose = FALSE)
  expect_equal(out_pre_rpv$Z, out_ns$Z, tolerance = 1e-6)

  # Explicitly setting nonstationary = FALSE overrides local_smoothness for enhancement
  out_override <- ptfce(z_img, mask, local_smoothness = sm_ns,
                        sided = "two.sided", nonstationary = FALSE, verbose = FALSE)
  expect_false(out_override$nonstationary)
  expect_equal(out_override$Z, out_stat$Z, tolerance = 1e-12)
})

test_that("Non-stationary mode works cleanly on whole-brain RNifti 4D residuals and warns on 3D fallback", {
  res4d_path <- system.file("extdata", "test03A_res4d.nii.gz", package = "pTFCE")
  mask_path  <- system.file("extdata", "test03A_mask.nii.gz", package = "pTFCE")
  z_path     <- system.file("extdata", "test03B_zstat1.nii.gz", package = "pTFCE")
  dof_val    <- as.numeric(readLines(system.file("extdata", "test03A_dof.txt", package = "pTFCE"), warn = FALSE)[1L])

  r4d <- RNifti::readNifti(res4d_path)
  m3a <- RNifti::readNifti(mask_path)
  z3b <- RNifti::readNifti(z_path)

  res_ns <- ptfce(z3b, m3a, residual = r4d, dof = dof_val,
                  sided = "two.sided", nonstationary = TRUE, verbose = FALSE)
  expect_true(res_ns$nonstationary)
  expect_s3_class(res_ns$Z, "niftiImage")
  expect_true(all(is.finite(as.array(res_ns$Z))))
  expect_true(all(as.array(res_ns$Z)[as.array(m3a) <= 0.5] == 0))

  # Requesting nonstationary = TRUE without 4D residuals or local_smoothness warns clearly
  d_small <- c(14L, 14L, 14L)
  z_small <- array(stats::rnorm(prod(d_small)), dim = d_small)
  m_small <- array(1, dim = d_small)
  expect_warning(
    ptfce(z_small, m_small, nonstationary = TRUE, verbose = FALSE),
    "4D residuals are strongly recommended"
  )
})
