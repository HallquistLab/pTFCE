# Test cases for GRF FWER conversion helpers

test_that("FWER p to Z is accurate and vectorized", {
  # Exact roots of resel_count * 0.11694 * exp(-0.5*z^2) * (z^2 - 1) = FWEP
  expect_equal(fwe.p2z(100, 0.05), 4.046618, tolerance = 1e-6)
  expect_equal(fwe.p2z(10, 0.01), 3.842918, tolerance = 1e-6)
  expect_equal(fwe.p2z(100000, 0.00001), 7.037302, tolerance = 1e-6)

  # Vectorized evaluation and exact round-trip inversion with fwe.z2p
  alphas <- c(0.05, 0.025, 0.01, 0.001, 1e-6)
  zs <- fwe.p2z(650, alphas)
  expect_length(zs, length(alphas))
  expect_equal(fwe.z2p(650, zs), alphas, tolerance = 1e-12)

  # Extreme small FWEP beyond IEEE exp(-0.5*z^2) double underflow (z > 38)
  z_extreme <- fwe.p2z(1000, 1e-300)
  expect_true(is.finite(z_extreme) && z_extreme > 37)
})

test_that("FWER Z to p is accurate and vectorized", {
  expect_equal(fwe.z2p(100, 2.3), 1, tolerance = 1e-6)
  expect_equal(fwe.z2p(100, 5), 0.001045908, tolerance = 1e-6)
  expect_equal(fwe.z2p(10000, 6), 0.0006233478, tolerance = 1e-6)

  # Vectorized input including Z < 2
  vals <- fwe.z2p(100, c(1.0, 2.3, 5.0))
  expect_equal(vals, c(1.0, 1.0, 0.001045908), tolerance = 1e-6)
})
