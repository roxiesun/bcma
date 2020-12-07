###init value 1 for parameters

B <- matrix(c(
  1,
  0,
  0
), nrow = p, ncol = q, byrow = TRUE)

L.me <- B

psi <- rep(1, p)  

Alpha <- matrix(rep(0, w*(1 + s + z)), nrow = w, byrow = TRUE)
psw <- array(0,dim = c(w,w))
diag(psw) <- 1


Delta <- matrix(rep(0, q*(1 + s + z + w)), nrow = q, ncol = (1 + z + s + w), byrow = TRUE)


L.se <- cbind(Delta)

psd <- rep(1, q)  

gam.z <- c(0, 0, 0, 0)
gam.s <- c(0, 0)
gam.w <- c(0, 0)
gam.m <- 0
gam <- c(gam.z, gam.s, gam.w, gam.m)


lam0 <- rep(1, G)

# to facilitate computation later 
iv.psi <- 1/psi
iv.sqrt.psi <- sqrt(iv.psi)

iv.psd <- 1/psd
iv.sqrt.psd <- sqrt(iv.psd)

iv.psw <- chol2inv(chol(psw))
c.iv.psw <- chol(iv.psw)   
