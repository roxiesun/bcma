##########prior setting

Alpha0 <- matrix(rep(0, 2*(1 + z + s)), nrow = w, ncol = (1 + z + s), byrow = TRUE)
sig.alpha <- rep(0.001, n.Alpha/w)

rho.scale <- 3
rho0 <- rho.scale + w + 1
R0 <- rho.scale*diag(1, w)*0.5   # Note, this R0 is used in Inv-Wishart dist

Delta0 <- matrix(rep(0, 1 + z + s + w), nrow = q, ncol = (1 + z + s + w), byrow = TRUE)

sig.delta <- rep(0.001,n.se)         # prior precision for β

L.se0 <- cbind(Delta0)



alpha.psd <- 9      # α0_epsilon
beta.psd <- 4       # β0_epsilon


B0 <- matrix(c(
  0,
  0,
  0
), nrow = p, ncol = q, byrow = TRUE)
sig.b <- 0.001                     

L.me0 <- B0

alpha.psi <- 9      # α0_ζ
beta.psi <- 4       # β0_ζ

gam0 <- rep(0, s + z + w + q)
sig.gamma <- rep(0.001,n.gam)          



#### For the piecewise constant λ0
#### we choose G = 5, thus λ = (λ1, λ2，λ3，λ4，λ5), corresponding to time interval
####                          [d1 = 0,d2],(d2,d3],(d3,d4],(d4,d5],(d5,d6 = max(T)]; 
####                          with dj = (j - 1)/5th quantile of T, vector d = (d1, d2, ..., d6), 
####                          interval length corresponding to λj is d[j+1] - d[j] 


alpha.lam0 <- 1        # prior shape parameter for λj in piecewise constant baseline
beta.lam0 <- 0.01      # prior rate parameter for λj in piecewise constant baseline

