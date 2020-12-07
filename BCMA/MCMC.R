##### for storing results
result <- list('B' = array(NA, dim = c(n.kept/n.thin, n.B)),
               'Alpha' = array(NA, dim = c(n.kept/n.thin, n.Alpha)),
               'L.se' = array(NA, dim = c(n.kept/n.thin, n.se)),
               'gam' = array(NA, dim = c(n.kept/n.thin, n.gam)),
               'psd' = array(NA, dim = c(n.kept/n.thin, q)),
               'psi' = array(NA, dim = c(n.kept/n.thin, n.psi)), 
               'psw' = array(NA, dim = c(n.kept/n.thin, w^2)),
               'lam0' = array(NA, dim = c(n.kept/n.thin, G)),
               'med' = array(NA, dim = c(n.kept/n.thin, n.med)))
n.accept.gam = 1
n.accept.M = rep(1,n)
for (it in 1:n.mcmc) {
  its <- it - n.burn
  
  # step1: update latent mediators M
  if (q > 0) {
    PI0 <- (diag(1, q) - PI)
    ISG <- crossprod(iv.sqrt.psi * B) + crossprod(iv.sqrt.psd * PI0)     # ∑_ω^-1
    SIG <- chol2inv(chol(ISG))                     
    SIG <- sig.mh*SIG
    cSIG <- chol(SIG)
  }
  #### Random Walk Metropolis, M_t + N[0, σ_mh^2*∑_ω_t]
  M.new <- M + t(crossprod(cSIG, matrix(rnorm(q*n), nrow = q)))
  ll1 <- Loglike(B, psi, lam0, gam, d, u, Delta,PI, psd, V, M, Z, S, X, W)
  ll2 <- Loglike(B, psi, lam0, gam, d, u, Delta, PI, psd, V, M.new, Z, S, X, W)
  ## acceptance ratio
  p.accept <- exp(ll2 - ll1)
  accept <- (runif(n) < p.accept)  
  M[accept, ] <- M.new[accept, ]
  n.accept.M <- n.accept.M + accept

  # step2: update loading matrix B and residual variance psi
  count.B <- 1
  L.me <- B
  omg.me <- M
  for (k in 1:p) {
    free <- Id.me[k, ]                  # which column on row K is free/fixed
    len <- n.me.row[k]
    Vcen <- t(L.me[k, !free, drop = F] %*% t(omg.me)[!free, , drop = F])
    Mk <- t(omg.me)[free, , drop = F]   # data of the latent variable that the kth indicator loads on
    Psiginv <- rep(sig.b, len)          # H_0yk


    Vk.star <- V[ , k] - Vcen
    alpha.psi.star <- alpha.psi + 0.5*n
    beta.psi.star <- beta.psi + 0.5*sum(Vk.star^2)
    if (len > 0) {
      A_vk <- chol2inv(chol(diag(Psiginv, len) + tcrossprod(Mk)))
      temp <- Psiginv*L.me0[k, free] + Mk %*% Vk.star
      a_vk <- A_vk %*% temp
      beta.psi.star <- beta.psi.star + 0.5*(sum(L.me0[k, free]*Psiginv*L.me0[k, free]) - sum(temp*a_vk))
    }

    iv.psi[k] <- rgamma(1, shape = alpha.psi.star, rate = beta.psi.star)
    psi[k] <- 1/iv.psi[k]                  
    iv.sqrt.psi[k] <- sqrt(iv.psi[k])

    if (len > 0) {
      L.me[k, free] <- rmvnorm(1, a_vk, sigma = psi[k]*A_vk)
      if (n.b.row[k] > 0) {
        B[k, ] <- L.me[k, ]
      }
      if (its > 0 && its %% n.thin == 0) {
        if (n.b.row[k] > 0) {
          result$B[its/n.thin, count.B:(count.B + n.b.row[k] - 1)] <- B[k, Id.B[k, ]]
        }
        
      }
      count.B <- count.B + n.b.row[k]
    }
    
  }

  # step3: update regression coefficients delta and mean_m μ and error variance psd
  ## Note!! Need revise if the number of latent mediator changed!
  count.se <- 1
  L.se <- cbind(Delta)
  omg.se <- cbind(rep(1,n), Z, S, W)
  for (k in 1:q) {
    free <- Id.se[k, ]  # which column on row K is free/fixed
    len <- n.se.row[k]
    Mcen <- t(L.se[k, !free, drop = F] %*% t(omg.se)[!free, , drop = F])
    Mk.star <- M[ , k] - as.vector(Mcen)
    alpha.psd.star <- alpha.psd + 0.5*n
    beta.psd.star <- beta.psd + 0.5*sum(Mk.star^2)

    if (len > 0) {
      Yk <- omg.se[, free, drop = F]  
      iH0dk <- diag(sig.delta, len)
      L.se0k <- L.se0[k, free]
      A_dk <- chol2inv(chol(iH0dk + crossprod(Yk)))
      temp <- iH0dk %*% L.se0k + t(Yk) %*% Mk.star
      a_dk <- A_dk %*% temp
      beta.psd.star <- beta.psd.star + 0.5*(crossprod(crossprod(iH0dk, L.se0k), L.se0k)
                                    - crossprod(a_dk, temp))
    }

    iv.psd[k] <- rgamma(1, shape = alpha.psd.star, rate = beta.psd.star)
    psd[k] <- 1/iv.psd[k]
    iv.sqrt.psd[k] <- sqrt(iv.psd[k])

    if (len > 0) {
      L.se[k, free] <- rmvnorm(1, a_dk, sigma = psd[k]*A_dk)
      if (its > 0 && its %% n.thin == 0) {
        result$L.se[its / n.thin, count.se:(count.se + len - 1)] <- L.se[k, free]
      }
      count.se <- count.se + len
    }
  }

  Delta <- L.se[, 1:(1 + z + s + w), drop = F]

  if (its > 0 && its %% n.thin == 0) {
    result$psd[its / n.thin, ] <- psd
    result$psi[its / n.thin, ] <- psi
  }

  # step4: update regression coefficients gam, M-H I guess
  #### Random Walk Metropolis, gam_t + N[0, σ_mh^2*I}

 gam.new <- gam + rmvnorm(1, rep(0, n.gam), sigma = diag(rep(0.0015, n.gam)))
 llg1 <- Loglike_GA(lam0, gam, d, u, M, Z, S, X, W)
 llg2 <- Loglike_GA(lam0, gam.new, d, u, M, Z, S, X, W)
 ## acceptance ratio
 p.accept <- exp(llg2 - llg1)
 accept <- (runif(1) < p.accept)
 if (accept) gam <- gam.new
 gam.z <- gam[1:z]
 gam.s <- gam[(z + 1):(z + s)]
 gam.w <- gam[(z + s + 1):(z + s + w)]
 gam.m <- gam[(z + s + w + 1):(z + s + w + q)]
 n.accept.gam <- n.accept.gam + accept
 
  if (its > 0 && its %% n.thin == 0) {
    result$gam[its / n.thin, ] <- gam
  }
  
  
  # step5: update piecewise constant hazard lam0, conjugate gamma prior?
   
  alpha.lam0.star <- rep(alpha.lam0, G) + as.vector(t(X[,2]) %*% u)

  temp2 <- array(0, dim = c(n, G))

  for (j in 1:G) {
    temp3 <- 0
        if (j < G) {
          for (g in (j + 1):G) {
            temp3 <- temp3 + u[, g]*(d[j + 1] - d[j])
          }
        }
    temp2[,j] <- exp(cbind(Z, S, W, M) %*% as.vector(gam))*(u[,j]*(X[,1] - d[j]) + temp3)
  }

  beta.lam0.star <- rep(beta.lam0, G) + colSums(temp2)

  for (j in 1:G) {
    lam0[j] <- rgamma(1, shape = alpha.lam0.star[j], rate = beta.lam0.star[j])
  }
  
   if (its > 0 && its %% n.thin == 0) {
     result$lam0[its / n.thin, ] <- lam0
   }
  
  temp <- c(gam.s, (gam.w + gam.m*Delta[1, 8:9]) %*% Alpha[,6:7],
            (gam.w[1] + gam.m*Delta[1, 8]) %*% Alpha[1,6:7],
            (gam.w[2] + gam.m*Delta[1, 9]) %*% Alpha[2,6:7],
            Delta[1,6:7]*as.vector(gam.m))
  if (its > 0 && its %% n.thin == 0) {
    result$med[its / n.thin, ] <- temp
  }
  
  
  #step6 update Alpha, regression coefficients of W on Z, S
  Dat6 <- cbind(rep(1,n), Z, S)
  Pswginv <- diag(sig.alpha, n.Alpha/w)  # this is precision marix, what we use
  temp1.1 <- chol2inv(chol(crossprod(Dat6) + Pswginv))
  temp1.2 <- crossprod(Dat6, W) + Pswginv %*% t(Alpha0)
  Bn <- temp1.1 %*% temp1.2
  temp1 <- crossprod(W - Dat6 %*% Bn)
  temp2 <- (t(Bn) - Alpha0) %*% Pswginv %*% (Bn - t(Alpha0))
  iv.psw <- rWishart(5, rho0 + n , chol2inv(chol(R0 + temp1 + temp2)))[, , 1]
  c.iv.psw <- chol(iv.psw)
  psw <- chol2inv(c.iv.psw)
  
  library(matrixsampling)
  Alpha <- t(matrix(rmatrixnormal(1, Bn, temp1.1, psw), ncol = w))
  if (its > 0 && its %% n.thin == 0) {
    result$psw[its / n.thin, ] <- as.vector(psw)
    result$Alpha[its / n.thin, ] <- as.vector(t(Alpha))
  }
  
  
  if (it %% 100 == 0) cat(it,',')
} # end of mcmc
