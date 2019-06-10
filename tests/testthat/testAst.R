context("ast")
library(sn)

# note K(341) = 0.3986499, K(342) = 0, K(343) = NaN
# gamma(171) = 7.257416e+306, gamma(172) = Inf
test_that("K function of mu = ", {
  expect_equal(K(1), 0.3183099)
  expect_equal(K(341), 0.3986499)
  # is.na(K(343)) out of bound
})

test_that("density of mu = , sigma = , alpha = , nu1 = , nu2 = ", {
  expect_equal(dast(0, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), 0)
  expect_equal(dast(0, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), 0)
  expect_equal(dast(0, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), 0)
})

test_that("probability of mu = , sigma = , alpha = , nu1 = , nu2 = ", {
  expect_equal(past(0, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), 0)
  expect_equal(past(1.5, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), 0.8)
  expect_equal(past(Inf, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), 1)
})

test_that("quantile of mu = , sigma = , alpha = , nu1 = , nu2 = ", {
  expect_equal(qast(0, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), -Inf)
  expect_equal(qast(0.8, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), 1.5)
  expect_equal(qast(1, mu = 1.5, sigma = 1.2, alpha = 0.8, nu1 = 3, nu2 = 4), Inf)
})

test_that("random generation of mu = , sigma = , alpha = , nu1 = , nu2 = ", {
  # tbd
})

