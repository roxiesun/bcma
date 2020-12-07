##########prior setting2

Alpha0 <- matrix(c(
  0, -1, 1, 0, 0, -0.5, -1,
  0, 1, -1, 0, 0, 0.5, 1
), byrow = TRUE, nrow = w, ncol = (1 + z + s))
sig.alpha <- 1/c(1000,2,2,1000,1000,1,1)

rho.scale <- 3
rho0 <- rho.scale + w + 1
R0 <- rho.scale*diag(1, w)*0.5   # Note, this R0 is used in Inv-Wishart dist


Delta0 <-matrix(c(
  0, 1, -1, 0, 0, 0.5, 1, -1, 1
), nrow = q, ncol = (1 + z + s + w), byrow = TRUE)


sig.delta <- 1/c(1000,2,2,1000,1000,1,1,1,1)         # prior precision for δ


L.se0 <- cbind(Delta0)



alpha.psd <- 9      # α0_eps
beta.psd <- 4       # β0_eps


B0 <- matrix(c(
  1,
  1,
  -1
), nrow = p, ncol = q, byrow = TRUE)
sig.b <- 1                     # precision = 1/var of normal prior for B

#L.me0 <- cbind(mu0, B0)
L.me0 <- B0

alpha.psi <- 9      # α0_ζ
beta.psi <- 4       # β0_ζ


gam0 <- c(1,-1,0,0,0.5,1,-1,1,1)
sig.gamma <- 1/c(2,2,1000,1000,1,1,1,1)          # prior precision for γ



#### For the piecewise constant λ0
#### we choose G = 5, thus λ = (λ1, λ2，λ3，λ4，λ5), corresponding to time interval
####                          [d1 = 0,d2],(d2,d3],(d3,d4],(d4,d5],(d5,d6 = max(T)]; 
####                          with dj = (j - 1)/5th quantile of T, vector d = (d1, d2, ..., d6), 
####                          interval length corresponding to λj is d[j+1] - d[j] 


alpha.lam0 <- 1     # prior shape parameter for λj in piecewise constant baseline
beta.lam0 <- 0.001      # prior rate parameter for λj in piecewise constant baseline
