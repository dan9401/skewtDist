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
#
# pgat = function(x, location, scale, dof, skewness, tpa, sa){
#   r = tpa; c = sa
#   q = 1 / (1 + c^(-skewness*(1 + r^2)/r)*(((x-location)/scale)+sqrt(1+(x-location)^2/scale^2))^(-skewness*(1+r^2)/r))
#   pbeta(q, dof/skewness/(1+r^2), r^2*dof/skewness/(1+r^2))
# }
#
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
