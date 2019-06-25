# # # # Clean Session
# # rm(list = ls())
# # library(devtools)
# # dev_mode(on = T)
# # # Session Info
# # sessionInfo()
# # # remove dev packages
# # pkgs = installed.packages()
# # dev_pkgs = pkgs[, 1][pkgs[,2] == '/Users/zhiyexia/R-dev']
# # remove.packages(dev_pkgs)
# #
# # ## oo notes
# # # jia class
# # a = c(1,2,3)
# # is.object(a)
# # # in pryr
# # # otype for objects
# # # ftype for functions
# # # some functions do dispatch in R using UseMethod()
# # # while some other dipatches in C, they're called internal generics
# # methods(mean)
# # methods(class = 'ts')
# # # (Apart from methods defined in the base package, most S3 methods will not be visible: use getS3method() to read their
# # # source code.)
# # test <- function(x) { structure(list(), class = 'test') }
# # runner = test()
# # class(runner)
# # t.test <- function(x) 'Class test'
# # t.test(runner)
# # y <- 1
# # g <- function(x) { y <- 2
# # UseMethod('g') }
# # g.numeric <- function(x) y
# # g(10)
# # h <- function(x) { x <- 10
# # UseMethod('h') }
# # h.character <- function(x) paste('char', x)
# # h.numeric <- function(x)
# # paste('num', x)
# # h('a')
# # h(10)
# # f <- function() 'this is a foo func'
# # f1 <- list(a = c(1,2), b = function(x) 'foo')
# # g <- function() 2
# # class(g) <- 'function'
# # class(f)
# # class(f1)
# # class(g)
# # length.function <- function(x) 'function'
# # length(f)
# # length(f1)
# # length(g)
# # library(stats4)
# # # From example
# # (mle) y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
# # nLL <- function(lambda) - sum(dpois(y, lambda, log = TRUE))
# # fit <- mle(nLL, start = list(lambda = 5), nobs = length(y))
# # # An S4 object
# # isS4(fit)
# # is(fit)
# # getClasses()
# # getGenerics()
#
# # # test file # optimizer # example test # f(x, y, z) = 3 * x^3 - 2x^2 + xz - 2y^2 # x in (1,6) # y in (-3, 1) # z in (-2,
# # 5) fn = function(pars) { x = pars[1] y = pars[2] z = pars[3] } # example 2 f <- function(x) 2*(x[1]-1)^2 + 5*(x[2]-3)^2 +
# # 10 r <- optim(c(1, 1), f) r$convergence == 0 r$par r$value # ast problem # sigma > 0 # alpha in (0,1) # nu1 and nu2 > 0
# # library(nloptr) ## Rosenbrock Banana function eval_f <- function(x) { return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
# # ) } ## Gradient of Rosenbrock Banana function eval_grad_f <- function(x) { return( c( -400 * x[1] * (x[2] - x[1] * x[1])
# # - 2 * (1 - x[1]), 200 * (x[2] - x[1] * x[1]) ) ) } x0 <- c( -1.2, 1 ) opts <- list('algorithm'='NLOPT_LD_LBFGS',
# # 'xtol_rel'=1.0e-8) # solve Rosenbrock Banana function res <- nloptr( x0=x0, eval_f=eval_f, eval_grad_f=eval_grad_f,
# # opts=opts) print( res ) ## Rosenbrock Banana function and gradient in one function eval_f_list <- function(x) {
# # common_term <- x[2] - x[1] * x[1] return( list( 'objective' = 100 * common_term^2 + (1 - x[1])^2, 'gradient' = c( -400 *
# # x[1] * common_term - 2 * (1 - x[1]), 200 * common_term) ) ) } res <- nloptr( x0=x0, eval_f=eval_f_list, opts=opts) print(
# # res ) K = function(nu) { gamma(0.5*(nu+1)) / (sqrt(pi*nu)*gamma(0.5*nu)) } llast = function(pars, y) { mu = pars[1] sigma
# # = pars[2] alpha = pars[3] nu1 = pars[4] nu2 = pars[5] T_ = length(y) y1 = y[y <= mu] y2 = y[y > mu] logl = T_ *
# # log(sigma) + 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) + 0.5 * (nu2 + 1) * sum(log(1
# # + ((y2 - mu)/(2 * (1-alpha) * sigma * K(nu2)))^2/nu2)) logl } y = rast(1000, 1.5, 1, 0.8, 3, 4) llast(c(1.5, 1, 0.8, 3,
# # 4), y) restest = nloptr(x0 = c(0,1,0.5,1,1), eval_f = llast, lb = c(-Inf, 0, 0, 0, 0), ub = c(Inf, Inf, 1, Inf, Inf),
# # opts = list('algorithm'='NLOPT_LN_COBYLA', 'maxeval' = 100000, 'xtol_rel'=1.0e-8), y = y) restest
#
# ############################################################
# # with the ln_cobyla algo
# # 1. the dof parameter at the side of skewness would be a problem
# # 2. sigma and 2 dof would be more inaccurate
# # 3. requires more time with dof smaller than 1
# # 4. the est. time seems to be less when dof is larger
# # 5. (1, 1.5) seems to be an exception for no.4
# # 6. (1, 2) too
# # 7.
# ############################################################
#
#
#
# infoMat_t <- function(sigma, nu) {
#   I11 <- 1/4*(trigamma((nu+1)/2) - trigamma(nu/2)) - 1/nu*(1/(nu+1) - 1/(2*(nu+3)) )
#   I12 <- 1/sigma*(1/(nu+3)-1/(nu+1))
#   I22 <- 2/sigma^2*nu/(nu+3)
#   I13 <- I23 <- 0
#   I33 <- 1/sigma^2*(1-2/(nu+3))
#   infoMat <- matrix(c(I11, I12, I13,
#                       I12, I22, I23,
#                       I13, I23, I33),
#                     nrow = 3, ncol = 3)
#   rownames(infoMat) = colnames(infoMat) = c("nu", "sigma", "mu")
#   infoMat
# }
#
# est_idx = c(1,3,5)
# ipars = data.frame(fixed_pars = c(1,2,3,4,5))
# rownames(ipars) = c("alpha", "mu", "nu1", "nu2", "sigma")
# solver_control <- list('algorithm' = 'NLOPT_LN_COBYLA', 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8)
#
# data = rast(1000, 0.12, 0.6, 0.7, 3, 5)
# #spec = astspec(data)
# #fixed_names<- rownames(ipars)[-est_idx]
# ipars = data.frame(fixed_pars = c(1,2,3,4,5))
# rownames(ipars) = c("alpha", "mu", "nu1", "nu2", "sigma")
# fixed_names <- c("alpha")
#
# data = rast(1000, 0.12, 0.6, 0.7, 3, 5)
# if ("nu1" %in% fixed_names) {
#   nu1vec <- ipars["nu1", "fixed_pars"]
# } else {
#   nu1vec <- seq(2, 20, by = 4)
# }
# if ("nu2" %in% fixed_names) {
#   nu2vec <- ipars["nu2", "fixed_pars"]
# } else {
#   nu2vec <- seq(2, 20, by = 4)
# }
# grid <- array(c(rep(nu1vec, length(nu2vec)), rep(nu2vec, length(nu1vec))), c(length(nu1vec), length(nu2vec), 2))
# if (dim(grid)[1] == 5 & dim(grid)[2] == 5) {
#   grid[,,2] = t(grid[,,2])
# }
#
# valueGrid <- apply(grid, 1:2, objective_value, data = data, solver = "nloptr", solver_control = solver_control)
# idx <- as.numeric(which(valueGrid == min(valueGrid), arr.ind = T))
# idx <- c(grid[idx[1], idx[2], 1], grid[idx[1], idx[2], 2])
#
# if (!("nu1" %in% fixed_names)) { nu1vec <- seq(idx[1]-2, idx[1]+2, by=1) }
# if (!("nu2" %in% fixed_names)) { nu2vec <- seq(idx[2]-2, idx[2]+2, by=1) }
# grid <- array(c(rep(nu1vec, length(nu2vec)), rep(nu2vec, length(nu1vec))), c(length(nu1vec), length(nu2vec), 2))
# if (dim(grid)[1] == 5 & dim(grid)[2] == 5) {
#   grid[,,2] = t(grid[,,2])
# }
# valueGrid <- apply(grid, 1:2, objective_value, data = data, solver = "nloptr", solver_control = solver_control)
# idx <- as.numeric(which(valueGrid == min(valueGrid), arr.ind = T))
# idx <- c(grid[idx[1], idx[2], 1], grid[idx[1], idx[2], 2])
#
#
#
#
# solver <- "nloptr"
# spec <- astspec(data, fixed_pars = c("nu1" = idx[1], "nu2" = idx[2]))
# fit <- astfit(spec, solver, solver_control)
# fit$fitted_pars
#
# start_pars <- c(fit$sol_res$solution, c(idx[1], idx[2]))
# names(start_pars) <- c("alpha", "mu", "sigma", "nu1", "nu2")
# spec <- astspec(data, start_pars = start_pars)
# fit <- astfit(spec, solver, solver_control)
# fit$fitted_pars
#
#
#
#
#
# objective_value <- function(nus, data, solver, solver_control) {
#   if (nus[1] == 0 || nus[2] == 0) {
#     obj = 10^5
#   } else {
#     spec <- astspec(data, fixed_pars = c("nu1" = nus[1], "nu2" = nus[2]))
#     obj <- astfit(spec, solver, solver_control)$sol_res$objective
#   }
#   obj
# }
#
#
#
# # if ("nu1" %in% fixed_names) {
# #   nu1vec <- ipars["nu1", "fixed_pars"]
# # } else {
# #   nu1vec <- seq(1, 20, by = 1)
# # }
# # if ("nu2" %in% fixed_names) {
# #   nu2vec <- ipars["nu2", "fixed_pars"]
# # } else {
# #   nu2vec <- seq(1, 20, by = 1)
# # }
# # grid <- array(c(rep(nu1vec, length(nu2vec)), rep(nu2vec, length(nu1vec))), c(length(nu1vec), length(nu2vec), 2))
# # if (dim(grid)[1] == dim(grid)[2] ) {
# #   grid[,,2] = t(grid[,,2])
# # }
# # valueGrid <- apply(grid, 1:2, objective_value, data = data, solver = "nloptr", solver_control = solver_control)
# # idx <- as.numeric(which(valueGrid == min(valueGrid), arr.ind = T))
# # idx <- c(grid[idx[1], idx[2], 1], grid[idx[1], idx[2], 2])
#
#
#
#
# nuTest = function(n, mu, sigma, alpha, nu1, nu2) {
#   tt = Sys.time()
#   data <- rast(n, mu, sigma, alpha, nu1, nu2)
#   grid <- matrix(nrow =5, ncol = 5)
#   solver_control <- list('algorithm' = 'NLOPT_LN_COBYLA', 'maxeval' = 1.0e4, 'xtol_rel' = 1.0e-8)
#
#   for (i in seq(2,20,by=4)) {
#     for (j in seq(2,20,by=4)) {
#       spec_t <- astspec(data, fixed_pars = c("nu1" = i, "nu2" = j))
#       grid[(i+2)/4, (j+2)/4] <- astfit(spec_t, 'nloptr', solver_control)$sol_res$objective
#     }
#   }
#   idx1 = as.numeric(which(grid == min(grid), arr.ind = T)) * 4 - 2
#   spec_1 <- astspec(data, fixed_pars = c("nu1" = idx1[1], "nu2" = idx1[2]))
#   fit_1 <- astfit(spec_1, 'nloptr', solver_control)
#   tt1 = Sys.time()
#   res1 = c(fit_1$fitted_pars[c(1,2)], idx1, fit_1$fitted_pars[3],
#            fit_1$sol_res$iterations, difftime(tt1, tt, units = "secs"))
#   names(res1)[c(3,4,6,7)] = c("nu1", "nu2", "iter", "time")
#
#   seq1 = seq(idx1[1]-2,idx1[1]+2,by=1)
#   seq2 = seq(idx1[2]-2,idx1[2]+2,by=1)
#   ii = 1
#   for (i in seq1) {
#     ij = 1
#     for (j in seq2) {
#       if (i ==0 || j == 0) {
#         grid[ii, ij] = 10^5
#       } else {
#         spec_t <- astspec(data, fixed_pars = c("nu1" = i, "nu2" = j))
#         grid[ii, ij] <- astfit(spec_t, 'nloptr', solver_control)$sol_res$objective
#       }
#       ij = ij + 1
#     }
#     ii = ii + 1
#   }
#   idx_t = as.numeric(which(grid == min(grid), arr.ind = T))
#   idx2 = c(seq1[idx_t[1]], seq2[idx_t[2]])
#   spec_2 <- astspec(data, fixed_pars = c("nu1" = idx2[1], "nu2" = idx2[2]))
#   fit_2 <- astfit(spec_2, 'nloptr', solver_control)
#   tt2 = Sys.time()
#   res2 = c(fit_2$fitted_pars[c(1,2)], idx2, fit_2$fitted_pars[3],
#            fit_2$sol_res$iterations, difftime(tt2, tt1, units = "secs"))
#   names(res2)[c(3,4,6,7)] = c("nu1", "nu2", "iter", "time")
#
#   spec_3 <- astspec(data, start_pars = c(fit_2$fitted_pars, c("nu1" = idx2[1], "nu2" = idx2[2])))
#   #solver_control_n <- list('algorithm' = 'NLOPT_LN_COBYLA', 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8)
#   fit_3 <- astfit(spec_3, 'nloptr', solver_control)
#   #fit_n$sol_res
#   #fp3 = fit_n$fitted_pars
#   #fi3 = fit_n$sol_res$iterations
#   tt3 = Sys.time()
#   #t3 = difftime(tt3, tt1, units = "secs")
#   res3 = c(fit_3$fitted_pars,
#            fit_3$sol_res$iterations, difftime(tt3, tt2, units = "secs"))
#   names(res3)[c(6,7)] = c("iter", "time")
#   print(fit_3$sol_res$message)
#   print(fit_3$sol_res$objective)
#   deviation <- c(alpha = (res3["alpha"] - alpha)/alpha, mu = (res3["mu"] - mu)/mu,
#                  nu1 = (res3["nu1"] - nu1)/nu1, nu2 = (res3["nu2"]- nu2)/nu2,
#                  sigma = (res3["sigma"] - sigma)/sigma, iter = NA, time = NA)
#
#   res = t(data.frame(grid1 = res1, grid2 = res2, result = res3, deviation = deviation))
#   return(fit_3)
# }
# #
# # mu = 0.12; sigma = 0.6; nu1 = 3; nu2 = 5
# # n = 1000
# # alpha = 0.6
# # nuTest(n, mu, sigma, alpha, nu1, nu2)
# #
# # fit_t = nuTest(1000, 0.12, 0.6, 0.5, 3, 5)
# #
# # avg = numeric(5)
# # names(avg) = c("alpha", "mu", "nu1", "nu2", "sigma")
# # for (i in 1:10) {
# #   avg = avg + abs(nuTest(1000, 0.12, 0.6, 0.5, 3, 5)["deviation", 1:5])
# # }
# # avg/10
# #
# #
# # # grid surface
# # data <- rast(n, mu, sigma, alpha, nu1, nu2)
# # grid <- matrix(nrow = 20, ncol = 20)
# # solver_control <- list('algorithm' = 'NLOPT_LN_COBYLA', 'maxeval' = 1.0e4, 'xtol_rel' = 1.0e-8)
# #
# # for (i in 1:20) {
# #   for (j in 1:20) {
# #     spec_t <- astspec(data, fixed_pars = c("nu1" = i, "nu2" = j))
# #     grid[i, j] <- astfit(spec_t, 'nloptr', solver_control)$sol_res$objective
# #   }
# # }
# # idx = as.numeric(which(grid == min(grid), arr.ind = T))
# # spec <- astspec(data, fixed_pars = c("nu1" = idx[1], "nu2" = idx[2]))
# # fit <- astfit(spec, 'nloptr', solver_control)
# #
# # library(plotly)
# # Sys.setenv("plotly_username"="dan94")
# # Sys.setenv("plotly_api_key"="Iwc6tVefUrrO1muMvIkB")
# # p <- plot_ly(z = ~grid) %>% add_surface()
# # chart_link = api_create(p, filename="grid")
# # chart_link
# # #
# # #
# #
# #
# # #
# # fit_t = nuTest(10000, 0.12, 0.6, 0.5, 3, 5)
# # infoMat_ast(fit_t)
# # infoMat_ast(fit_t, method = "observed")
#
#
# qqplot_asta <- function(fit, dist = "normal", ...) {
#   y <- fit$data
#   pars <- fit$fitted_pars
#   mu <- pars["mu"]
#   sigma <- pars["sigma"]
#   alpha <- pars["alpha"]
#   nu1 <- pars["nu1"]
#   nu2 <- pars["nu2"]
#
#   p <- past(y, mu, sigma, alpha, nu1, nu2)
#   if (dist == "normal") {
#     x <- qnorm(p, mean = mu, sd = sigma)
#     px <- qnorm(c(0.25, 0.75), mean = mu, sd = sigma)
#   } else if (dist == "t") {
#     df = 0.5 * (nu1 + nu2)
#     x <- qt(p, df = df)
#     px <- qast(c(0.25, 0.75), mu, sigma, 0.5, df, df)
#   } else if (dist == "ast") {
#     x <- qast(p, mu, sigma, alpha, nu1, nu2)
#     px <- qast(c(0.25, 0.75), mu, sigma, alpha, nu1, nu2)
#   }
#   qqplot(x, y, ...)
#   px = quantile(x, c(0.25, 0.75))
#   py <- quantile(y, c(0.25, 0.75))
#   slope <- diff(py)/diff(px)
#   int <- py[1L] - slope * px[1L]
#   abline(int, slope)
#   points(px, py, col = 2)
# }
#
#
#
