### generate initial values for latent mediator M1
source('init1.R')
#source('init2.R')
#source('init3.R')

if (q > 0) {
  const <- rep(1,n)
  Dat <- cbind(const, Z, S, W)
  Coe <- cbind(Delta)
  M <- Dat %*% t(Coe) + rmvnorm(n, mean = rep(0, q), sigma = diag(psd, q))
  M <- M %*% t(solve((diag(1, q) - PI)))
    
}
  

