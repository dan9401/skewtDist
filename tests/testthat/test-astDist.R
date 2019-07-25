test_that("dast", {
  mu = 0.12; sigma = 0.6; alpha = 0.6; nu1 = 3; nu2 = 5
  pars = c(mu, sigma, alpha, nu1, nu2)
  expect_equal(round(dast(0, 0.12, 0.6, 0.6, 3, 5),6), 1.459715)
  expect_equal(round(dast(0, pars = pars),6), 1.459715)
  expect_equal(dast(0, 0.12, 0.6, 0.6, 3, 5), dast(0, pars = pars))
})

test_that("past", {
  mu = 0.12; sigma = 0.6; alpha = 0.6; nu1 = 3; nu2 = 5
  pars = c(mu, sigma, alpha, nu1, nu2)
  expect_equal(round(dast(0, 0.12, 0.6, 0.6, 3, 5),6), 1.459715)
  expect_equal(round(dast(0, pars = pars),6), 1.459715)
  expect_equal(dast(0, 0.12, 0.6, 0.6, 3, 5), dast(0, pars = pars))
})

test_that("qast", {
  mu = 0.12; sigma = 0.6; alpha = 0.6; nu1 = 3; nu2 = 5
  pars = c(mu, sigma, alpha, nu1, nu2)
  expect_equal(round(dast(0, 0.12, 0.6, 0.6, 3, 5),6), 1.459715)
  expect_equal(round(dast(0, pars = pars),6), 1.459715)
  expect_equal(dast(0, 0.12, 0.6, 0.6, 3, 5), dast(0, pars = pars))
})

test_that("rast", {
  mu = 0.12; sigma = 0.6; alpha = 0.6; nu1 = 3; nu2 = 5
  pars = c(mu, sigma, alpha, nu1, nu2)
  expect_equal(round(dast(0, 0.12, 0.6, 0.6, 3, 5),6), 1.459715)
  expect_equal(round(dast(0, pars = pars),6), 1.459715)
  expect_equal(dast(0, 0.12, 0.6, 0.6, 3, 5), dast(0, pars = pars))
})

