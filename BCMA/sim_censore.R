#############
#simulate data with censoring
source('def.con.R')
source('true.par.R')
library(mvtnorm)

for (rep in 1:rep_num) {
  set.seed(201 + rep*10)
  
  ################################
  ######  data generation   ######
  S <- t(rmultinom(n, 1, prob = c(0.5, 0.2, 0.3)))
  S <- S[, 2:3]
  Z <- rmvnorm(n, mean = rep(0,z), sigma = diag(1,z))
  M <- array(NA, dim = c(n, q))
  V <- array(NA, dim = c(n, p))
  W <- array(NA, dim = c(n, w))
  X <- array(NA, dim = c(n, 2))
  T <- array(NA, dim = c(n, 1))
  C <- runif(n, 0, c[2])
  
  
  Dat <- cbind(rep(1,n), Z, S)
  W <- Dat %*% t(Alpha) + rmvnorm(n, mean = rep(0,w), sigma = psw)
  Dat <- cbind(rep(1,n), Z, S, W)
  Coe <- cbind(Delta)
  M <- Dat %*% t(Coe) + rmvnorm(n, mean = rep(0, q), sigma = diag(psd, q))
  M.true <- M <- M %*% t(solve((diag(1, q) - PI)))
  const <- rep(1,n)
  Dat1 <- cbind(M)
  V <- Dat1 %*% t(L.me) + rmvnorm(n, mean = rep(0, p), sigma = diag(psi))
  
  # censoring time and survival time
  Dat2 <- cbind(Z, S, W, M)
  temp <- exp(Dat2 %*% gam)
  # with lambda_0(t) = 1
  T <- (-1/(1*temp))*log(runif(n))
  C <- runif(n, min = 0, max = c[1])
  # with lamda_0(t) = 2t + 1
  #T <- sqrt((-1/temp)*log(runif(n)) + 1/4) - 1/2
  #C <- runif(n, min = 0, max = quantile(T,0.985))
  
  X[, 1] <- (T <= C)*T + (T > C)*C
  X[, 2] <- ifelse(T <= C, 1, 0)
    
  
 
    
  write(S, File.S, ncolumns = s, append = TRUE)  
  write(Z, File.Z, ncolumns = z, append = TRUE)
  write(X, File.X, ncolumns = 2, append = TRUE)
  write(V, File.V, ncolumns = p, append = TRUE)
  write(W, File.W, ncolumns = w, append = TRUE)
  
  print(paste(rep,'th dataset of sample size ', n, ' is generated and stored'))
  
}

info <- paste('Done!! ', rep_num, 'datasets of sample size', n, 'are generated', sep = ' ') 
print(info)



  
  
  

