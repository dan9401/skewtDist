test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


pars = c(0.12, 0.6, 0.6, 3, 5)
data = rast(10^6, pars = pars)
infoMat_ast(pars = pars, method = "expected")
infoMat_ast(pars = pars, data = data, method = "observed")

nu12
mu2
# small change
alpha2
sigma2


# auto test for bad fit
# infoMat symbolic
# added
#

