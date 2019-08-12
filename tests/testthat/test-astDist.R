test_that("dast", {
  mu <- 0.12; sigma <- 0.6; alpha <- 0.6; nu1 <- 3; nu2 <- 5
  pars <- c(mu, sigma, alpha, nu1, nu2)
  x <- 0; answer <- 1.459715
  expect_equal(round(dast(x, mu, sigma, alpha, nu1, nu2), 6), answer)
  expect_equal(round(dast(x, pars = pars), 6), answer)
  expect_equal(dast(x, mu, sigma, alpha, nu1, nu2), dast(x, pars = pars))
})

test_that("past", {
  mu <- 0.12; sigma <- 0.6; alpha <- 0.6; nu1 <- 3; nu2 <- 5
  pars <- c(mu, sigma, alpha, nu1, nu2)
  q <- 0; answer <- 0.408609
  expect_equal(round(past(q, mu, sigma, alpha, nu1, nu2), 6), answer)
  expect_equal(round(past(q, pars = pars), 6), answer)
  expect_equal(past(q, mu, sigma, alpha, nu1, nu2), past(q, pars = pars))
})

test_that("past", {
  mu <- 0.12; sigma <- 0.6; alpha <- 0.6; nu1 <- 3; nu2 <- 5
  pars <- c(mu, sigma, alpha, nu1, nu2)
  p <- 0.408609; answer <- 0
  expect_equal(round(qast(p, mu, sigma, alpha, nu1, nu2), 6), answer)
  expect_equal(round(qast(p, pars = pars), 6), answer)
  expect_equal(qast(p, mu, sigma, alpha, nu1, nu2), qast(p, pars = pars))
})

test_that("rast", {
  mu = 0.12; sigma = 0.6; alpha = 0.6; nu1 = 3; nu2 = 5
  pars = c(mu, sigma, alpha, nu1, nu2)
  expect_equal(round(dast(0, 0.12, 0.6, 0.6, 3, 5),6), 1.459715)
  expect_equal(round(dast(0, pars = pars),6), 1.459715)
  expect_equal(dast(0, 0.12, 0.6, 0.6, 3, 5), dast(0, pars = pars))
})

