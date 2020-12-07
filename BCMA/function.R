##### some useful functions



##### calculate log-likelihood of the full conditional dist. of M|X, V, theta, & observed data
Loglike <- function(B, psi, lam0, gam, d, u, Delta, PI = matrix(0, ncol = q, nrow = q), psd, V, M, Z, S, X, W){
  
    short <- array()
    invP  <- solve(diag(psi))
    gam.z <- gam[1:z]
    gam.s <- gam[(z + 1):(z + s)]
    gam.w <- gam[(z + s + 1):(z + s + w)]
    gam.m <- gam[(z + s + w + 1):(z + s + w + q)]
    for (i in 1:n) {
      short[i] <- (V[i, ] - M[i, ] %*% t(B)) %*% invP %*% t(V[i,]  - M[i, ] %*% t(B))
    }
    log.like1 <- -0.5*(p*log(2*pi) + log(abs(det(diag(psi)))) + short)
    
    
    short1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      part1 <- u[,j]*X[,2]*(log(lam0[j]) + Z %*% gam.z + S %*% gam.s + W %*% gam.w + M %*% gam.m)
      
      temp <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          temp <- temp + as.numeric(lam0[g]*(d[g + 1] - d[g]))
        }
      }
      part2 <- -u[,j]*(lam0[j]*(X[,1] - d[j]) + temp)*exp(Z %*% gam.z + S %*% gam.s + W %*% gam.w + M %*% gam.m)
      short1[, j] <- part1 + part2
    }
    log.like2 <- apply(short1, 1, sum)
    
    
    const <- rep(1,n)
    Dat1 <- cbind(const, Z, S, W)
    Coe1 <- cbind(Delta)
    Mcen <- Dat1 %*% t(Coe1) 
    invPI0 <- solve(diag(1, q) - PI)
    short2 <- array()
    SIG <- invPI0 %*% diag(psd,q) %*% t(invPI0)
    ISIG <- solve(SIG)
    for (i in 1:n) {
      short2[i] <- t(M[i, ] - invPI0 %*% Mcen[i, ]) %*% ISIG %*% (M[i, ] - invPI0 %*% Mcen[i, ])
    }
    log.like3 <- -0.5*(q*log(2*pi) + log(abs(det(SIG))) + short2)
    
    log.like <- log.like1 + log.like2 + log.like3
    return(log.like)
}





### calculate the loglikelihood of full conditional dist. of gam|X, theta, & observed data

Loglike_GA <- function(lam0, gam, d, u, M, Z, S, X, W){
  short1 <- array(0, dim = c(n, G))
  gam.z <- gam[1:z]
  gam.s <- gam[(z + 1):(z + s)]
  gam.w <- gam[(z + s + 1):(z + s + w)]
  gam.m <- gam[(z + s + w + 1):(z + s + w + q)]
  
  sum1 <- 0
  for (j in 1:G) {
    part1 <- u[,j]*X[,2]*(log(lam0[j]) + Z %*% gam.z + S %*% gam.s + W %*% gam.w + M %*% gam.m)
    
    temp <- 0
    if (j > 1) {
      for (g in 1:(j - 1)) {
        temp <- temp + as.numeric(lam0[g]*(d[g + 1] - d[g]))
      }
    }
    part2 <- -u[,j]*(lam0[j]*(X[,1] - d[j]) + temp)*exp(Z %*% gam.z + S %*% gam.s + W %*% gam.w + M %*% gam.m)
    sum1 <- sum1 + sum(part1) + sum(part2) 
  }
  log.like1 <- sum1

  invP  <- diag(sig.gamma, n.gam)
  gam <- matrix(gam, ncol = 1)
  short2 <- t(gam - gam0) %*% invP %*% (gam - gam0)
  log.like2 <- -0.5*((n.gam)*log(2*pi) + log(abs(det(chol2inv(chol(invP))))) + short2)

  log.like <- log.like1 + log.like2
  return(log.like)
}

