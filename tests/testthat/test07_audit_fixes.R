# Tests for v0.4.2 audit fixes:
#   1. External RPV.nii / dLh map global Riemannian Rd & resels extraction without V/Rd (Recipe C)
#   2. Exact 2D (D = 2) GRF cluster-extent, Euler characteristic FWER, and 2D pTFCE calibration
#   3. Small-ROI (n_act < nt) and dof-mismatch guards in sign_flip
#   4. Single-map 3D local roughness winsorization & monotone low-DOF interpolate()

test_that("External RPV.nii and dLh maps derive global Rd and resels automatically without V/Rd", {
  set.seed(402)
  d <- c(24L, 20L, 14L)
  nt <- 24L
  mask <- array(1, dim = d)
  res4d <- array(stats::rnorm(prod(d) * nt), dim = c(d, nt))
  z_img <- array(stats::rnorm(prod(d)), dim = d)
  z_img[9:12, 8:11, 6:8] <- 3.8

  sm_ns <- smoothest(res4d, mask, dof = nt - 1L, nonstationary = TRUE, verbose = FALSE)
  r_list <- ptfce(z_img, mask, local_smoothness = sm_ns, sided = "two.sided", verbose = FALSE)

  # 1. Passing in-memory rpv_map or dLh_map directly WITHOUT V, Rd, or resels
  r_rpv_mem <- ptfce(z_img, mask, local_smoothness = sm_ns$rpv_map,
                     sided = "two.sided", verbose = FALSE)
  r_dlh_mem <- ptfce(z_img, mask, local_smoothness = sm_ns$dLh_map,
                     sided = "two.sided", verbose = FALSE)

  expect_equal(r_rpv_mem$number_of_resels, r_list$number_of_resels, tolerance = 1e-6)
  expect_equal(r_rpv_mem$fwer0.05.Z, r_list$fwer0.05.Z, tolerance = 1e-6)
  expect_equal(as.array(r_rpv_mem$Z), as.array(r_list$Z), tolerance = 1e-6)
  expect_equal(as.array(r_dlh_mem$Z), as.array(r_list$Z), tolerance = 1e-6)

  # 2. Passing an external NIfTI file path "RPV.nii.gz" WITHOUT V, Rd, or resels (Recipe C)
  tmp_rpv <- file.path(tempdir(), "RPV.nii.gz")
  on.exit(unlink(tmp_rpv), add = TRUE)
  RNifti::writeNifti(as.array(sm_ns$rpv_map), tmp_rpv)

  set.seed(77)
  r_list_sim <- ptfce(z_img, mask, local_smoothness = sm_ns, sided = "two.sided",
                      fwe_method = "grf_sim", n_perm = 49L, verbose = FALSE)
  set.seed(77)
  r_file_rpv <- ptfce(z_img, mask, local_smoothness = tmp_rpv, sided = "two.sided",
                      fwe_method = "grf_sim", n_perm = 49L, verbose = FALSE)

  expect_equal(r_file_rpv$number_of_resels, r_list$number_of_resels, tolerance = 1e-5)
  expect_equal(as.array(r_file_rpv$Z), as.array(r_list$Z), tolerance = 1e-5)
  expect_gt(r_file_rpv$fwer0.05.Z, 2.5)
})

test_that("2D (D = 2) GRF Es, fwe.p2z/fwe.z2p, and 2D ptfce work accurately", {
  # 2D Euler characteristic round-trip inversion
  alphas <- c(0.05, 0.025, 0.01, 1e-4)
  z2d <- fwe.p2z(200, alphas, D = 2L)
  expect_equal(fwe.z2p(200, z2d, D = 2L), alphas, tolerance = 1e-12)
  # 2D threshold at 200 resels differs from 3D threshold at 200 resels
  expect_false(isTRUE(all.equal(z2d[1L], fwe.p2z(200, 0.05, D = 3L))))

  # 2D expected cluster size matches exact analytical formula V * Phi(-h) / (Rd * (2*pi)^(-1.5) * h * exp(-0.5*h^2))
  h <- 2.5
  V <- 4096
  Rd <- 120
  es_2d_exact <- (V * stats::pnorm(h, lower.tail = FALSE)) / (Rd * ((2 * pi)^(-1.5)) * h * exp(-0.5 * h^2))
  expect_equal(Es(h, V, Rd, D = 2L), es_2d_exact, tolerance = 1e-12)

  # Full 2D image pipeline (stationary + grf_sim)
  set.seed(202)
  d2 <- c(40L, 40L)
  m2 <- matrix(1, nrow = 40L, ncol = 40L)
  z2 <- matrix(stats::rnorm(1600L), nrow = 40L, ncol = 40L)
  z2[18:22, 18:22] <- 3.6

  expect_warning(
    res_2d <- ptfce(z2, m2, sided = "two.sided", fwe_method = "grf_sim",
                    n_perm = 49L, verbose = FALSE),
    "Using 2D image mode"
  )
  expect_equal(dim(res_2d$Z), d2)
  expect_true(all(is.finite(res_2d$Z)))
  expect_true(is.finite(res_2d$fwer0.05.Z) && res_2d$fwer0.05.Z > 2.5)
})

test_that("Small-ROI rank cap and dof mismatch emit informative warnings in sign_flip", {
  set.seed(303)
  d <- c(10L, 10L, 10L)
  mask_tiny <- array(0, dim = d)
  mask_tiny[4:6, 4:5, 4:5] <- 1 # only 12 active voxels!
  nt <- 20L
  res_tiny <- array(stats::rnorm(prod(d) * nt), dim = c(d, nt))
  z_tiny <- array(0, dim = d)
  z_tiny[mask_tiny == 1] <- stats::rnorm(sum(mask_tiny))

  expect_warning(
    ptfce(z_tiny, mask_tiny, residual = res_tiny, dof = 19L,
          sided = "two.sided", fwe_method = "sign_flip", n_perm = 29L, verbose = FALSE),
    "Active mask voxel count"
  )

  # Full mask with severely mismatched supplied dof
  mask_full <- array(1, dim = d)
  expect_warning(
    ptfce(array(stats::rnorm(prod(d)), dim = d), mask_full,
          residual = res_tiny, dof = 10L, # true centered rank is 19
          sided = "two.sided", fwe_method = "sign_flip", n_perm = 29L, verbose = FALSE),
    "Empirical residual-subspace rank"
  )
})

test_that("3D single-map local smoothness resists mean-shift self-suppression and interpolate() is monotone", {
  set.seed(404)
  d <- c(24L, 24L, 16L)
  mask <- array(1, dim = d)
  z_map <- array(stats::rnorm(prod(d)), dim = d)
  # Add a strong +15 mean activation shift on top of the unit-variance noise field
  z_map[9:15, 9:15, 6:10] <- z_map[9:15, 9:15, 6:10] + 15.0

  sm_3d <- smoothest(z_map, mask, nonstationary = TRUE, smooth_fwhm = 4.0, verbose = FALSE)
  # Interior cluster noise differences (Z_u - Z_v) equal (eps_u - eps_v), so weight stays ~ 1.0 (not 0.02!)
  expect_gt(sm_3d$voxel_weights[12, 12, 8], 0.5)

  # Low-DOF interpolate() is monotone non-increasing across dof = 3..10
  dofs <- 3:10
  facs <- vapply(dofs, pTFCE:::interpolate, numeric(1L))
  expect_true(all(diff(facs) <= 0))
})
