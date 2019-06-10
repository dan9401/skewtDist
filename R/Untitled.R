# # Clean Session
# rm(list = ls())
# library(devtools)
# dev_mode(on = T)
#
# # Session Info
# sessionInfo()
#
# # remove dev packages
# pkgs = installed.packages()
# dev_pkgs = pkgs[, 1][pkgs[,2] == "/Users/zhiyexia/R-dev"]
# remove.packages(dev_pkgs)
#
# # Session Info
# sessionInfo()
#
# # install.packages()
#
#
# # Session Info
# sessionInfo()
#
# # load packages
# # library()
#
# ## oo notes
# # jia class
# a = c(1,2,3)
# is.object(a)
#
# # in pryr
# # otype for objects
# # ftype for functions
#
# # some functions do dispatch in R using UseMethod()
# # while some other dipatches in C, they're called internal generics
#
# methods(mean)
# methods(class = "ts")
# # (Apart from methods defined in the base package, most S3 methods will not be visible: use getS3method() to read their source code.)
#
# test <- function(x) {
#   structure(list(), class = "test")
# }
# runner = test()
# class(runner)
#
# t.test <- function(x) "Class test"
#
# t.test(runner)
#
# y <- 1
# g <- function(x) {
#   y <- 2
#   UseMethod("g")
# }
# g.numeric <- function(x) y
# g(10)
#
#
# h <- function(x) {
#   x <- 10
#   UseMethod("h")
# }
# h.character <- function(x) paste("char", x)
# h.numeric <- function(x) paste("num", x)
#
# h("a")
# h(10)
#
#
# f <- function() "this is a foo func"
# f1 <- list(a = c(1,2), b = function(x) "foo")
# g <- function() 2
# class(g) <- "function"
#
# class(f)
# class(f1)
# class(g)
#
# length.function <- function(x) "function"
# length(f)
# length(f1)
# length(g)
#
# library(stats4)
#
# # From example(mle)
# y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
# nLL <- function(lambda) - sum(dpois(y, lambda, log = TRUE))
# fit <- mle(nLL, start = list(lambda = 5), nobs = length(y))
#
# # An S4 object
# isS4(fit)
#
# is(fit)
# getClasses()
# getGenerics()
