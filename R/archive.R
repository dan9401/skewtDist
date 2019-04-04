# ## functions under development
#
#
# dast = function(x, location, scale, skewness, dof1, dof2) {
#   K = function(dof) {
#     gamma(0.5*(dof+1)) / (sqrt(pi*dof)*gamma(0.5*dof))
#   }
#   alpha_star = skewness*K(dof1) / (skewness * K(dof1)+(1-skewness)*K(dof2))
#   x1 = x[x <= location]
#   x2 = x[x > location]
#   dens = numeric(length(x))
#   dens[x <= location] = K(dof1) * skewness/ alpha_star / scale * (1 + ((x1-location)/(2*alpha_star*scale))^2/dof1)^(-0.5*(dof1+1))
#   dens[x > location] = K(dof2) * (1-skewness) / (1- alpha_star)/ scale * (1 + ((x2-location)/(2*(1-alpha_star)*scale))^2/dof1)^(-0.5*(dof2+1))
#   dens
# }
#
# dast2 = function(x, location, scale, skewness, dof1, dof2) {
#   K = function(dof) {
#     gamma(0.5*(dof+1)) / (sqrt(pi*dof)*gamma(0.5*dof))
#   }
#   alpha_star = skewness*K(dof1) / (skewness * K(dof1)+(1-skewness)*K(dof2))
#   x1 = x[x <= location]
#   x2 = x[x > location]
#   dens = numeric(length(x))
#   dens[x <= location] = K(dof1) * skewness/ alpha_star / scale * (1 + ((x1-location)/(2*alpha_star*scale))^2/dof1)^(-0.5*(dof1+1))
#   dens[x > location] = K(dof2) * (1-skewness) / (1- alpha_star)/ scale * (1 + ((x2-location)/(2*(1-alpha_star)*scale))^2/dof1)^(-0.5*(dof2+1))
#   dens
# }
#
# past = function(x, location, scale, skewness, dof1, dof2) {
#   x1 = x[x <= location]
#   x2 = x[x > location]
#   probs = numeric(length(x))
#   if(length(x1) > 0) {
#     probs[x <= location] = sapply(x1, function(y) integrate(dast, -Inf, y, location, scale, skewness, dof1, dof2)$value)
#   }
#   if(length(x2) > 0) {
#     probs[x > location] = integrate(dast, -Inf, location, location, scale, skewness, dof1, dof2)$value + sapply(x2, function(y) integrate(dast, location, y, location, scale, skewness, dof1, dof2)$value)
#   }
#   probs
# }
#
# past2 = function(x, location, scale, skewness, dof1, dof2) {
#   x1 = x[x <= location]
#   x2 = x[x > location]
#   probs = numeric(length(x))
#   if(length(x1) > 0) {
#     probs[x <= location] = sapply(x1, function(y) integrate(dast2, -Inf, y, location, scale, skewness, dof1, dof2)$value)
#   }
#   if(length(x2) > 0) {
#     probs[x > location] = integrate(dast2, -Inf, location, location, scale, skewness, dof1, dof2)$value + sapply(x2, function(y) integrate(dast, location, y, location, scale, skewness, dof1, dof2)$value)
#   }
#   probs
# }
#
#
# location = 0; scale = 1; dof1 = dof2 = 1; skewness = 0.5
# dast(0, 0, 1, 0.5, 1, 1)
# dast2(0, 0, 1, 0.5, 1, 1)
# dt(0, 1)
# a = seq(-7, 7, 0.01)
# plot(a, dast(a, 0, 1, 0.5, 1, 1))
# plot(a, dast2(a, 0, 1, 0.5, 1, 1))
# integrate(dast, -Inf, Inf, 0, 1, 0.5, 1, 1)
# integrate(dast2, -Inf, Inf, 0, 1, 0.5, 1, 1)
# past(Inf,  0, 1, 0.5, 1, 1)
# past(0.1,    0.1, 1, 0.5, 1, 1)
# past(-Inf, 0, 1, 0.5, 1, 1)
#
#
# # know issues
# # first = 1, last = 1, a problem of integrate in R
# integrate(dast, -Inf, -Inf, 0, 1, 0.5, 1, 1)
# integrate(dast, -Inf, 0, 0, 1, 0.5, 1, 1)
# integrate(dast, -Inf, Inf, 0, 1, 0.5, 1, 1)
# integrate(dast, 0, 0, 0, 1, 0.5, 1, 1)
# integrate(dast, 0, Inf, 0, 1, 0.5, 1, 1)
# integrate(dast, Inf, Inf, 0, 1, 0.5, 1, 1)
# integrate(dast, Inf, -Inf, 0, 1, 0.5, 1, 1)
#
#
#
# dgat = function(x, location, scale, dof, skewness, tpa, sa) {
#   r = tpa; c = sa
#   g = (x-location)/scale + sqrt(1+((x-location)/scale)^2)
#   A = skewness * (1 + tpa^2) / (tpa * scale)
#   B = (((c*g)^(skewness*r))+((c*g)^(-skewness/r)))^(-dof/skewness)/beta(dof/skewness/(1+r^2), r^2*dof/skewness/(1+r^2))
#   C = (1+((x-location)/scale)^2)^(-0.5)
#   A * B * C
# }
#
# int_beta = function(q, a, b) {
#   q^(a-1)*(1-q)^(b-1)
# }
# i_Beta = function(c, a, b) {
#   integrate(int_beta, 0, c, a, b)$value
# }
# ri_Beta = function(c, a, b) {
#   i_Beta(c,a,b)/beta(a,b)
# }
# pbeta(0.5, 2, 3)
# ri_Beta(0.5, 2, 3)
#
#
# pgat = function(x, location, scale, dof, skewness, tpa, sa){
#   r = tpa; c = sa
#   q = 1 / (1 + c^(-skewness*(1 + r^2)/r)*(((x-location)/scale)+sqrt(1+(x-location)^2/scale^2))^(-skewness*(1+r^2)/r))
#   pbeta(q, dof/skewness/(1+r^2), r^2*dof/skewness/(1+r^2))
# }
#
# dgat(0, 0, 1, 1, 1, 1, 1)
# pgat(0, 0, 1, 1, 1, 1, 1)
#
# # random generation
# rgat = function(n, location, scale, dof, skewness, tpa, sa) {
#   r = tpa; c = sa
#   a = dof/skewness/(1+r^2)
#   b = dof*r^2/skewness/(1+r^2)
#   q = rbeta(n, a, b)
#   r = tpa; c = sa
#   del = r / skewness / (1 + r^2)
#   x = location + 0.5*scale*((q/(1-q))^del/c - c*(q/(1-q))^(-del))
#   x
# }
#
# res = rgat(10000, 0, 1, 1, 1, 1, 1)
# hist(res)
#
# a = seq(-10, 10, 0.01)
# plot(a, dgat(a, 0, 1, 1, 1, 1, 1))
#
#
#
