# Comprehensive tests for bi-sided pTFCE, numerical stability, masking, smoothness, and I/O

test_that("Bi-sided pTFCE is strictly anti-symmetric, preserves sub-ZestThr values, and has no dead zone", {
  z_path <- system.file("extdata", "test03B_zstat1.nii.gz", package = "pTFCE")
  m_path <- system.file("extdata", "test03B_mask.nii.gz", package = "pTFCE")

  z_nii <- RNifti::readNifti(z_path)
  m_nii <- RNifti::readNifti(m_path)
  inside <- as.array(m_nii) > 0.5
  z_arr <- as.array(z_nii)

  # Precompute smoothness once for speed
  sm <- smoothest(z_nii, m_nii, verbose = FALSE)
  V <- sm$volume
  Rd <- sm$dLh * V

  res_bi <- ptfce(z_nii, m_nii, Rd = Rd, V = V, resels = sm$resels,
                  sided = "two.sided", fwe_p = c(0.05, 0.01), verbose = FALSE)
  res_bi_flip <- ptfce(-z_arr, as.array(m_nii), Rd = Rd, V = V, resels = sm$resels,
                       sided = "two.sided", fwe_p = c(0.05, 0.01), verbose = FALSE)

  z_out <- as.array(res_bi$Z)
  p_out <- as.array(res_bi$p)
  lp_out <- as.array(res_bi$logp)

  # 1. Exact sign anti-symmetry when input contrast is negated
  expect_equal(z_out, -res_bi_flip$Z, tolerance = 1e-12)
  expect_equal(p_out, res_bi_flip$p, tolerance = 1e-12)

  # 2. Out-of-mask voxels are strictly zero in Z and logp, 1 in p (never -8.21)
  expect_true(all(z_out[!inside] == 0))
  expect_true(all(lp_out[!inside] == 0))
  expect_true(all(p_out[!inside] == 1))

  # 3. No 20% zero dead-zone: all non-zero sub-ZestThr voxels retain exact input Z
  sub_vox <- inside & (abs(z_arr) <= 1.3) & (z_arr != 0)
  expect_equal(z_out[sub_vox], z_arr[sub_vox], tolerance = 1e-12)

  # 4. Signs of all in-mask voxels match input Z signs (no sign inversions)
  expect_equal(sign(z_out[inside]), sign(z_arr[inside]))

  # 5. Consistency between signed Z, two-sided p, logp, and two-sided FWER threshold
  pos_check <- inside & (abs(z_out) > 0.1) & (abs(z_out) < 8)
  expect_equal(2 * stats::pnorm(abs(z_out[pos_check]), lower.tail = FALSE),
               p_out[pos_check], tolerance = 1e-10)
  expect_equal(-log(p_out[pos_check]), lp_out[pos_check], tolerance = 1e-10)
  expect_equal(res_bi$fwer0.05.Z, fwe.p2z(res_bi$number_of_resels, 0.025), tolerance = 1e-12)
  expect_named(res_bi$fwe.Z, c("0.05", "0.01"))
})

test_that("Default fixed-dh grid is invariant to distant extreme Z peaks (including Z > 37.5)", {
  d <- c(24L, 24L, 24L)
  mask <- array(1, dim = d)
  img1 <- array(0, dim = d)

  # Moderate cluster in corner A: 5x5x5 block of Z = 3.6
  img1[3:7, 3:7, 3:7] <- 3.6
  img1[5, 5, 5] <- 3.9

  # Identical cluster in corner A, plus an extreme Z = 55 peak in distant corner B
  # (Z = 55 previously crashed pnorm(max(img), lower.tail = FALSE) with Inf, and
  #  even with log.p = TRUE destroyed dh resolution under fixed Nh = 100)
  img2 <- img1
  img2[20:22, 20:22, 20:22] <- 25.0
  img2[21, 21, 21] <- 55.0

  res1 <- ptfce(img1, mask, V = prod(d), Rd = 1500, resels = 50,
                sided = "two.sided", verbose = FALSE)
  res2 <- ptfce(img2, mask, V = prod(d), Rd = 1500, resels = 50,
                sided = "two.sided", verbose = FALSE)

  # Cluster A must receive the exact same enhanced Z regardless of the distant Z = 55 peak!
  expect_equal(res1$Z[3:7, 3:7, 3:7], res2$Z[3:7, 3:7, 3:7], tolerance = 1e-12)
  expect_true(is.finite(res2$Z[21, 21, 21]) && res2$Z[21, 21, 21] > 50)
})

test_that("All-positive and all-negative maps work cleanly in two.sided mode", {
  d <- c(15L, 15L, 15L)
  mask <- array(0, dim = d)
  mask[3:13, 3:13, 3:13] <- 1

  # Strictly positive inside mask (e.g. group intercept contrast)
  img_pos <- array(0, dim = d)
  img_pos[mask == 1] <- 2.5
  img_pos[7:9, 7:9, 7:9] <- 4.2

  res_pos <- ptfce(img_pos, mask, V = sum(mask), Rd = 400, resels = 20,
                   sided = "two.sided", verbose = FALSE)
  expect_true(all(is.finite(res_pos$Z)))
  expect_true(all(res_pos$Z[mask == 1] > 0))
  expect_true(all(res_pos$Z[mask == 0] == 0))

  # Strictly negative inside mask
  res_neg <- ptfce(-img_pos, mask, V = sum(mask), Rd = 400, resels = 20,
                   sided = "two.sided", verbose = FALSE)
  expect_equal(res_neg$Z, -res_pos$Z, tolerance = 1e-12)
})

test_that("Log-space posterior quadrature prevents NaN and p > 1 at low smoothness & large cluster size", {
  # Previously: pvox.clust(300000, 100000, 20000, 1.31) returned 1.01e+19 -> -log(p) = -43.76 -> NaN!
  p_val <- pTFCE:::pvox.clust(300000, 100000, 20000, 1.31)
  nlp_val <- pTFCE:::neglog_pvox_clust_cpp(300000, 100000, 20000, 1.31, 1.3)

  expect_true(is.finite(nlp_val) && nlp_val > 0)
  expect_equal(nlp_val, 2.189116, tolerance = 1e-5)
  expect_true(p_val > 0 && p_val < 1)

  # High threshold where dclust * 1e40 previously underflowed to 0 (capped at 745)
  nlp_hi <- pTFCE:::neglog_pvox_clust_cpp(300000, 100000, 20000, 5.50, 1.3)
  expect_equal(nlp_hi, 1992.122, tolerance = 1e-4)

  # Full volume stress test with low smoothness and huge cluster border near Z = 1.31
  d <- c(30L, 30L, 30L)
  mask <- array(1, dim = d)
  img <- array(1.32, dim = d)
  img[1, , ] <- 0.5
  img[15, 15, 15] <- 5.0
  res <- ptfce(img, mask, V = 300000, Rd = 100000, resels = 3,
               sided = "two.sided", verbose = FALSE)
  expect_false(any(is.na(res$Z)))
  expect_true(all(is.finite(res$Z)))
})

test_that("Out-of-mask NAs/Infs and out-of-mask clusters are strictly isolated", {
  d <- c(16L, 16L, 16L)
  mask <- array(0, dim = d)
  mask[4:12, 4:12, 4:12] <- 1

  img <- array(0, dim = d)
  # Put NAs, Infs, and a huge supra-threshold bridge outside the ROI mask
  img[1, 1, 1] <- NA_real_
  img[1, 1, 2] <- Inf
  img[3, 4:12, 4:12] <- 10.0 # adjacent to ROI boundary, outside mask!
  img[4, 8, 8] <- 2.5        # single voxel inside ROI touching out-of-mask 10.0 slab

  # Reference where exterior is cleanly zero
  img_clean <- img
  img_clean[mask == 0] <- 0

  expect_silent(res_dirty <- ptfce(img, mask, V = sum(mask), Rd = 100, resels = 10,
                                   sided = "two.sided", verbose = FALSE))
  res_clean <- ptfce(img_clean, mask, V = sum(mask), Rd = 100, resels = 10,
                     sided = "two.sided", verbose = FALSE)

  expect_equal(res_dirty$Z, res_clean$Z, tolerance = 1e-12)
})

test_that("4D residual smoothness estimation and 3D trailing-1 dim work in R >= 4.2", {
  res4d_path <- system.file("extdata", "test03A_res4d.nii.gz", package = "pTFCE")
  mask_path <- system.file("extdata", "test03A_mask.nii.gz", package = "pTFCE")
  dof_val <- as.numeric(readLines(system.file("extdata", "test03A_dof.txt", package = "pTFCE"), warn = FALSE)[1L])

  r4d <- RNifti::readNifti(res4d_path)
  m3a <- RNifti::readNifti(mask_path)

  sm <- smoothest(r4d, m3a, dof = dof_val, verbose = FALSE)
  expect_true(is.finite(sm$dLh) && sm$dLh > 0)
  expect_true(is.finite(sm$resels) && sm$resels > 0)
  # test03A_res4d has zero-variance voxels inside m3a which smoothest excludes
  expect_true(sm$volume > 0 && sm$volume <= sum(as.array(m3a) > 0.5))

  # Passing residual = 4D image directly into ptfce() must not fail on R >= 4.2
  z_synth <- RNifti::updateNifti(as.array(r4d)[, , , 1L], template = m3a)
  res <- ptfce(z_synth, m3a, residual = r4d, dof = dof_val,
               sided = "two.sided", verbose = FALSE)
  expect_s3_class(res$Z, "niftiImage")
  expect_equal(res$number_of_resels, sm$volume / sm$resels, tolerance = 1e-10)

  # 3D array with trailing 4th dim of length 1 (Issue #19)
  arr_4d1 <- array(rnorm(12 * 12 * 12), dim = c(12L, 12L, 12L, 1L))
  m_3d <- array(1, dim = c(12L, 12L, 12L))
  sm_4d1 <- smoothest(arr_4d1, m_3d, verbose = FALSE)
  sm_3d <- smoothest(array(arr_4d1, dim = c(12L, 12L, 12L)), m_3d, verbose = FALSE)
  expect_equal(sm_4d1, sm_3d, tolerance = 1e-12)
})

test_that("Input classes (file path, RNifti, oro.nifti, array), connectivity, and t2z work", {
  z_path <- system.file("extdata", "test01A_in.nii.gz", package = "pTFCE")
  m_path <- system.file("extdata", "test01A_mask.nii.gz", package = "pTFCE")

  # File path -> returns niftiImage
  r_path <- ptfce(z_path, m_path, V = 8000, Rd = 20, resels = 100, verbose = FALSE)
  expect_s3_class(r_path$Z, "niftiImage")

  # Plain 3D array -> returns plain array
  nii_in <- RNifti::readNifti(z_path)
  nii_mk <- RNifti::readNifti(m_path)
  arr <- array(as.numeric(nii_in), dim = dim(nii_in))
  m_arr <- array(as.numeric(nii_mk), dim = dim(nii_mk))
  r_arr <- ptfce(arr, m_arr, V = 8000, Rd = 20, resels = 100, verbose = FALSE)
  expect_true(is.array(r_arr$Z) && !inherits(r_arr$Z, "niftiImage"))
  expect_equal(as.numeric(r_path$Z), as.numeric(r_arr$Z), tolerance = 1e-12)

  # oro.nifti -> returns oro.nifti S4 object when installed
  if (requireNamespace("oro.nifti", quietly = TRUE)) {
    z_oro <- oro.nifti::readNIfTI(z_path, reorient = FALSE)
    m_oro <- oro.nifti::readNIfTI(m_path, reorient = FALSE)
    r_oro <- ptfce(z_oro, m_oro, V = 8000, Rd = 20, resels = 100, verbose = FALSE)
    expect_true(methods::is(r_oro$Z, "nifti"))
  }

  # Diagonal connectivity: two corner-touching voxels merge under 26-connectivity but not 6-connectivity
  diag_img <- array(0, dim = c(10L, 10L, 10L))
  diag_mask <- array(1, dim = c(10L, 10L, 10L))
  diag_img[5, 5, 5] <- 3.5
  diag_img[6, 6, 6] <- 3.5
  r6 <- ptfce(diag_img, diag_mask, V = 1000, Rd = 200, resels = 20, connectivity = 6, verbose = FALSE)
  r26 <- ptfce(diag_img, diag_mask, V = 1000, Rd = 200, resels = 20, connectivity = 26, verbose = FALSE)
  expect_gt(r26$Z[5, 5, 5], r6$Z[5, 5, 5])

  # t2z log-space conversion handles extreme t-values without Inf
  z_vec <- t2z(c(-30, -2, 0, 2, 30), dof = 25)
  expect_true(all(is.finite(z_vec)))
  expect_equal(z_vec[1], -z_vec[5], tolerance = 1e-12)
  expect_equal(z_vec[3], 0)
  expect_equal(pnorm(z_vec[4], lower.tail = FALSE), pt(2, df = 25, lower.tail = FALSE), tolerance = 1e-12)
})
