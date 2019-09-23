test_that("astMLE",
          {
            pars <- c(0.12, 0.6, 0.6, 6, 5)
            data <- rast(1000, pars = pars)
            fit <- astMLE(data)
          })
