#' use package robustbase for least trimmed squared regression (LTSR)
library(robustbase)
library(MASS)
mux <- 0; muy <- 0
sxx <- 3; syy <- 1
sxy <- syx <- 0.5
# sx <- 3; sy <- 1
# rho <- 0.5
# nstep <- 10
nsample <- 10000

#' MCMC method
# p_xy <- function(z, mux, muy, rho, sx, sy, direct = 1) {
#     rnorm(1, ifelse(direct, muy + rho * sy / sx * (z - mux), 
#                             mux + rho * sx / sy * (z - muy)),
#              ifelse(direct, sqrt(1 - rho ^ 2) * sy,
#                             sqrt(1 - rho ^ 2) * sx)
#     )
# }
# 
# x <- vector(length = nsample)
# y <- x
# 
# for (i in 1:nsample)
#     for (j in 1:nstep) {
#         x[i] = p_xy(y[i], mux, muy, rho, sx, sy, 0)
#         y[i] = p_xy(x[i], mux, muy, rho, sx, sy, 1)
#     }
# 
xy <- mvrnorm(nsample, mu = c(mux, muy), Sigma = matrix(c(sxx, sxy, syx, syy), ncol = 2), empirical = T)
x <- xy[, 1]
y <- xy[, 2]
plot(x, y, col = rgb(1, 0, 0, 0.1), pch = 20, cex = 2)

mvn_kdensity <- kde2d(x, y, h = 6)
contour(mvn_kdensity, add = T, lcol = rgb(0, 0, 0, 0.5), wd = 3)
abline(coef(lm(y ~ x)), col = rgb(0, 1, 0, 0.5), lwd = 3)
#abline(0, rho * sy / sx, col = rgb(1, 0, 0, 0.5), lwd = 3) #for MCMC method
abline(0, sxy / (sqrt(sxx) * sqrt(syy)), col = rgb(1, 0, 0, 0.5), lwd = 3)

lts <- ltsReg(y ~ x)
abline(coef(lts), col = rgb(0, 0, 1, 0.5), lwd = 3)
