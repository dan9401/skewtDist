# # test file
#
# # optimizer
#
# # example test
# # f(x, y, z) = 3 * x^3 - 2x^2 + xz - 2y^2
# # x in (1,6)
# # y in (-3, 1)
# # z in (-2, 5)
#
# fn = function(pars) {
#   x = pars[1]
#   y = pars[2]
#   z = pars[3]
#
# }
#
# # example 2
# f <- function(x) 2*(x[1]-1)^2 + 5*(x[2]-3)^2 + 10
# r <- optim(c(1, 1), f)
# r$convergence == 0
# r$par
# r$value
#
# # ast problem
# # sigma > 0
# # alpha in (0,1)
# # nu1 and nu2 > 0
#
# library(nloptr)
#
# ## Rosenbrock Banana function
# eval_f <- function(x) {
#   return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
# }
# ## Gradient of Rosenbrock Banana function
# eval_grad_f <- function(x) {
#   return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#              200 * (x[2] - x[1] * x[1]) ) )
# }
#
#
# x0 <- c( -1.2, 1 )
#
#
# opts <- list("algorithm"="NLOPT_LD_LBFGS",
#              "xtol_rel"=1.0e-8)
#
#
# # solve Rosenbrock Banana function
# res <- nloptr( x0=x0,
#                eval_f=eval_f,
#                eval_grad_f=eval_grad_f,
#                opts=opts)
# print( res )
#
# ## Rosenbrock Banana function and gradient in one function
# eval_f_list <- function(x) {
#   common_term <- x[2] - x[1] * x[1]
#   return( list( "objective" = 100 * common_term^2 + (1 - x[1])^2,
#                 "gradient" = c( -400 * x[1] * common_term - 2 * (1 - x[1]),
#                                 200 * common_term) ) )
# }
# res <- nloptr( x0=x0,
#                eval_f=eval_f_list,
#                opts=opts)
# print( res )
#
# K = function(nu) {
#   gamma(0.5*(nu+1)) / (sqrt(pi*nu)*gamma(0.5*nu))
# }
#
# llast = function(pars, y) {
#   mu = pars[1]
#   sigma = pars[2]
#   alpha = pars[3]
#   nu1 = pars[4]
#   nu2 = pars[5]
#   T_ = length(y)
#   y1 = y[y <= mu]
#   y2 = y[y > mu]
#   logl = T_ * log(sigma) + 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) +
#     0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1-alpha) * sigma * K(nu2)))^2/nu2))
#   logl
# }
#
# y = rast(1000, 1.5, 1, 0.8, 3, 4)
# llast(c(1.5, 1, 0.8, 3, 4), y)
#
# restest = nloptr(x0 = c(0,1,0.5,1,1),
#                  eval_f = llast,
#                  lb = c(-Inf, 0, 0, 0, 0),
#                  ub = c(Inf, Inf, 1, Inf, Inf),
#                  opts = list("algorithm"="NLOPT_LN_COBYLA",
#                              "maxeval" = 100000,
#                              "xtol_rel"=1.0e-8),
#                  y = y)
#
# restest
#
#
#
