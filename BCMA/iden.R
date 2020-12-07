########## Model Identification #########
### for the loading matrix of CFA model

Id.B <- matrix(c(
  0,
  1,
  1
), nrow = p, ncol = q, byrow = TRUE )

Id.B <- (Id.B > 0)
## free loadings in each row of the loading matrix B
n.b.row <- rowSums(Id.B)  
n.B <- sum(Id.B)

###for the measurment equation
Id.me <- Id.B
n.me.row <- rowSums(Id.me)  
n.me <- sum(Id.me)


Id.psi <- rep(1, p)
Id.psi <- as.logical(Id.psi)
n.psi <- sum(Id.psi)


Id.Delta <- matrix(rep(1, q*(1 + s + z + w)), nrow  = q)
Id.Delta <- (Id.Delta > 0)

## for the structual equation
Id.se <- cbind(Id.Delta)
n.se.row <- rowSums(Id.se)  
n.se <- sum(Id.se)

Id.psd <- rep(1, q)
Id.psd <- as.logical(Id.psd)
n.psd <- sum(Id.psd)

### for the PH model
Id.gam <- as.logical(rep(1, s + z + w + q))
n.gam <- sum(Id.gam)

## for the multivariate regression
Id.Alpha <- matrix(rep(1, w*(1 + z + s)), nrow = w, ncol = 1 + z + s)
Id.Alpha <- (Id.Alpha > 0)
n.Alpha <- sum(Id.Alpha)

## for the causal effects
n.med <- 10
