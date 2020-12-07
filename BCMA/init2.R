###init value 2 for parameters

B <- matrix(c(
  1,
  2,
  2
), nrow = p, ncol = q, byrow = TRUE)

L.me <- B

psi <- rep(2, p)  

Alpha <- matrix(rep(2, w*(1 + s + z)), nrow = w, byrow = TRUE)
psw <- array(0,dim = c(w,w))
diag(psw) <- 2


Delta <- matrix(rep(2, q*(1 + s + z + w)), nrow = q, ncol = (1 + z + s + w), byrow = TRUE)

L.se <- cbind(Delta)


psd <- rep(2, q)  

gam.z <- c(2, 2, 2, 2)
gam.s <- c(2, 2)
gam.w <- c(2, 2)
gam.m <- 2
gam <- c(gam.z, gam.s, gam.w, gam.m)


lam0 <- rep(2, G)

# to facilitate computation later 
iv.psi <- 1/psi
iv.sqrt.psi <- sqrt(iv.psi)

iv.psd <- 1/psd
iv.sqrt.psd <- sqrt(iv.psd)

iv.psw <- chol2inv(chol(psw))
c.iv.psw <- chol(iv.psw)   
