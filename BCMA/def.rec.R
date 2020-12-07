####### Store the mean and sd. of each parameter in all replications
mean.B <- sd.B <- array(0, dim = c(rep_num, n.B))
mean.Alpha <- sd.Alpha <- array(0, dim = c(rep_num, n.Alpha))
mean.psw <- sd.psw <- array(0, dim = c(rep_num, w^2))
mean.psi <- sd.psi <- array(0, dim = c(rep_num, n.psi))
mean.L.se <- sd.L.se <- array(0, dim = c(rep_num, n.se))
mean.psd <- sd.psd <- array(0, dim = c(rep_num, n.psd))
mean.gam <- sd.gam <- array(0, dim = c(rep_num, n.gam))
mean.lam0 <- sd.lam0 <- array(0, dim = c(rep_num, G))
mean.med <- sd.med <- array(0, dim = c(rep_num, n.med))


####### Store the quantile of each parameter in all replications
q.B <- array(0, dim = c(rep_num, n.B*2))
q.psi <- array(0, dim = c(rep_num, n.psi*2))
q.Alpha <- array(0, dim = c(rep_num, n.Alpha*2))
q.psw <- array(0, dim = c(rep_num, w^2*2))
q.L.se <- array(0, dim = c(rep_num, n.se*2))
q.psd <- array(0, dim = c(rep_num, n.psd*2))
q.gam <- array(0, dim = c(rep_num, n.gam*2))
q.lam0 <- array(0, dim = c(rep_num, G*2))
q.med <- array(0, dim = c(rep_num, n.med*2))

