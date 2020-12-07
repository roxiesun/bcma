##########################
# true parameter values #

Alpha.true <- Alpha <- matrix(c(
  0.2, 0.5, 0, 0, 0, 0.6, 0.6,
  0, 0.6, 0, 0, 0, 0.5, 0.5
), byrow = TRUE, nrow = w, ncol = (1 + z + s))

psw.true <- psw <- matrix(c(
  0.2, 0.1,
  0.1, 0.2
), nrow = w)



Delta.true <- Delta <- matrix(c(
   0, 0.2, 0, 0.5, 0, 0.5, 0.5, 0.6, 0.6
), nrow = q, ncol = (1 + z + s + w), byrow = TRUE)

L.se.true <- L.se <- cbind(Delta.true)
psd.true <- psd <- rep(0.3, q)  

B.true <- B <- matrix(0, nrow = p, ncol = q) 
B.true[,1] <- B[,1] <- c(1, 0.9, 0.7)

L.me.true <- L.me <- B
psi.true <- psi <- rep(0.2, p)  

# for calculation

gam.z.true <- gam.z <- c(0.5, 0.5, 0, 0)
gam.s.true <- gam.s <- c(0.5, 1)
#gam.m.true <- gam.m <- c(0.5, 1)
gam.w <- c(0.6, 0.2)
gam.m <- 1
gam.true <- gam <- c(gam.z, gam.s, gam.w, gam.m)

lam0.true <- lam0 <- rep(1, G)

med.true <- med <- c(gam.s, (gam.w + as.vector(gam.m)*Delta[1, 8:9]) %*% Alpha[,6:7], 
                     (gam.w[1] + as.vector(gam.m)*Delta[1, 8]) %*% Alpha[1,6:7],
                     (gam.w[2] + as.vector(gam.m)*Delta[1, 9]) %*% Alpha[2,6:7],
                     Delta[1,6:7]*gam.m)
