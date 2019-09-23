test_that("dast", {
  mu <- 0.12; s <- 0.6; alpha <- 0.6; nu1 <- 5; nu2 <- 6
  pars <- c(mu, s, alpha, nu1, nu2)
  x <- 0; answer <- 1.487862
  expect_equal(round(dast(x, mu, s, alpha, nu1, nu2), 6), answer)
  expect_equal(round(dast(x, pars = pars), 6), answer)
})

test_that("past", {
  mu <- 0.12; s <- 0.6; alpha <- 0.6; nu1 <- 5; nu2 <- 6
  pars <- c(mu, s, alpha, nu1, nu2)
  q <- 0; answer <- 0.4073696
  expect_equal(round(past(q, mu, s, alpha, nu1, nu2), 7), answer)
  expect_equal(round(past(q, pars = pars), 7), answer)
})

test_that("qast", {
  mu <- 0.12; s <- 0.6; alpha <- 0.6; nu1 <- 5; nu2 <- 6
  pars <- c(mu, s, alpha, nu1, nu2)
  p <- 0.4073696; answer <- 0
  expect_equal(round(qast(p, mu, s, alpha, nu1, nu2), 7), answer)
  expect_equal(round(qast(p, pars = pars), 7), answer)
})



# # note K(341) = 0.3986499, K(342) = 0, K(343) = NaN
# # gamma(171) = 7.257416e+306, gamma(172) = Inf
# test_that("K function of mu = ", {
#   expect_equal(K(1), 0.3183099)
#   expect_equal(K(341), 0.3986499)
#   # is.na(K(343)) out of bound
# })
#
