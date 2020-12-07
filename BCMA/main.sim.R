library(mvtnorm)
source('def.con.R')        # define the constants
source('true.par.R')       # define true values of parameters
source('iden.R')           # define identification constraints
source('prior.R')          # set hyperparameters
source('function.R')       # some funtions that might be used later
source('def.rec.R')        # define some arrays for data storing

time.begin <- proc.time()

cat('Generating', rep_num, 'sets of data...  Please wait \n')
source('sim_censore.R')
for (crep in 1:rep_num) {
  source('readdata.R')    # read the crep_th dataset 
  source('initial.R')     # initialise parameters and the latent variable
  
  d <- c(quantile(X[ , 1], probs = seq(0, 1, 1/G)))
  d[1] <- 0
  u <- array(0,dim = c(n, G))
  for (j in 1:G) {
    al <- as.logical(d[j] < X[,1])*(X[,1] <= d[j + 1])
    u[, j] <- al
  }
  
  source('MCMC.R')
  cat('\r', 'the acceptance rate of gamma is', n.accept.gam/n.mcmc, '\n')
  cat('\r', crep, 'sets of data completed, ', rep_num - crep, 'sets remaining.')
  
  
  mean.B[crep, ] <- colMeans(result$B)
  mean.psi[crep, ] <- colMeans(result$psi)
  mean.L.se[crep, ] <- colMeans(result$L.se)
  mean.psd[crep, ] <- colMeans(result$psd)
  mean.gam[crep, ] <- colMeans(result$gam)
  mean.lam0[crep, ] <- colMeans(result$lam0)
  mean.med[crep, ] <- colMeans(result$med)
  mean.psw[crep, ] <- colMeans(result$psw)
  mean.Alpha[crep, ] <- colMeans(result$Alpha)
  
  
  
  sd.B[crep, ] <- apply(result$B, 2, sd)
  sd.psi[crep, ] <- apply(result$psi,2, sd)
  sd.L.se[crep, ] <- apply(result$L.se, 2, sd)
  sd.psd[crep, ] <- apply(result$psd, 2, sd)
  sd.gam[crep, ] <- apply(result$gam, 2, sd)
  sd.lam0[crep, ] <- apply(result$lam0, 2, sd)
  sd.med[crep, ] <- apply(result$med, 2, sd)
  sd.psw[crep, ] <- apply(result$psw, 2, sd)
  sd.Alpha[crep, ] <- apply(result$Alpha, 2, sd)
  
  q.B[crep, ] <- as.vector(apply(result$B, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.psi[crep, ] <- as.vector(apply(result$psi, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.L.se[crep, ] <- as.vector(apply(result$L.se, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.psd[crep, ] <- as.vector(apply(result$psd, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.gam[crep, ] <- as.vector(apply(result$gam, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.lam0[crep, ] <- as.vector(apply(result$lam0, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.med[crep, ] <- as.vector(apply(result$med, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.psw[crep, ] <- as.vector(apply(result$psw, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.Alpha[crep, ] <- as.vector(apply(result$Alpha, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  
}

source('true.par.R')

B.list <- list(Mean = colMeans(mean.B), Sd = colMeans(sd.B),
               Bias = colMeans(mean.B) - B[Id.B],
               Rms = sqrt(rowSums((t(mean.B) - B[Id.B])^2) / rep_num),
               p.CI95 = colMeans((q.B[ , (1:n.B)*2 - 1] - t(replicate(
                 rep_num, B[Id.B])) > 0) * (q.B[, (1:n.B)*2] - t(
                   replicate(rep_num, B[Id.B])) > 0) == 0))

psi.list <- list(Mean = colMeans(mean.psi), Sd = colMeans(sd.psi),
                 Bias = colMeans(mean.psi) - psi[Id.psi],
                 Rms = sqrt(rowSums((t(mean.psi) - psi[Id.psi])^2) / rep_num),
                 p.CI95 = colMeans((q.psi[ , (1:n.psi)*2 - 1] - t(replicate(
                   rep_num, psi[Id.psi])) > 0) * (q.psi[, (1:n.psi)*2] - t(
                     replicate(rep_num, psi[Id.psi])) > 0) == 0))

L.se.list <- list(Mean = colMeans(mean.L.se), Sd = colMeans(sd.L.se),
               Bias = colMeans(mean.L.se) - t(L.se)[t(Id.se)],
               Rms = sqrt(rowSums((t(mean.L.se) - t(L.se)[t(Id.se)])^2) / rep_num),
               p.CI95 = colMeans((q.L.se[ , (1:n.se)*2 - 1] - t(replicate(
                 rep_num, t(L.se)[t(Id.se)])) > 0) * (q.L.se[, (1:n.se)*2] - t(
                   replicate(rep_num, t(L.se)[t(Id.se)])) > 0) == 0))

psd.list <- list(Mean = colMeans(mean.psd), Sd = colMeans(sd.psd),
                 Bias = colMeans(mean.psd) - psd[Id.psd],
                 Rms = sqrt(rowSums((t(mean.psd) - psd[Id.psd])^2) / rep_num),
                 p.CI95 = rowMeans((q.psd[ , (1:n.psd)*2 - 1] - t(replicate(
                   rep_num, psd[Id.psd])) > 0) * (q.psd[, (1:n.psd)*2] - t(
                     replicate(rep_num, psd[Id.psd])) > 0) == 0))

gam.list <- list(Mean = colMeans(mean.gam), Sd = colMeans(sd.gam),
                 Bias = colMeans(mean.gam) - gam,
                 Rms = sqrt(rowSums((t(mean.gam) - gam)^2) / rep_num),
                 p.CI95 = colMeans((q.gam[ , (1:n.gam)*2 - 1] - t(replicate(
                   rep_num, gam)) > 0) * (q.gam[, (1:n.gam)*2] - t(
                     replicate(rep_num, gam)) > 0) == 0))

lam0.list <- list(Mean = colMeans(mean.lam0), Sd = colMeans(sd.lam0),
                  Bias = colMeans(mean.lam0) - lam0,
                  Rms = sqrt(rowSums((t(mean.lam0) - lam0)^2) / rep_num),
                  p.CI95 = colMeans((q.lam0[ , (1:G)*2 - 1] - t(replicate(
                    rep_num, lam0)) > 0) * (q.lam0[, (1:G)*2] - t(
                      replicate(rep_num, lam0)) > 0) == 0))

psw.list <- list(Mean = colMeans(mean.psw), Sd = colMeans(sd.psw),
                 Bias = colMeans(mean.psw) - as.vector(psw),
                 Rms = sqrt(rowSums((t(mean.psw) - as.vector(psw))^2) / rep_num),
                 p.CI95 = colMeans((q.psw[ , (1:w^2)*2 - 1] - t(replicate(
                   rep_num, as.vector(psw))) > 0) * (q.psw[, (1:w^2)*2] - t(
                     replicate(rep_num, as.vector(psw))) > 0) == 0))

Alpha.list <- list(Mean = colMeans(mean.Alpha), Sd = colMeans(sd.Alpha),
                  Bias = colMeans(mean.Alpha) - t(Alpha)[t(Id.Alpha)],
                  Rms = sqrt(rowSums((t(mean.Alpha) - t(Alpha)[t(Id.Alpha)])^2) / rep_num),
                  p.CI95 = colMeans((q.Alpha[ , (1:n.Alpha)*2 - 1] - t(replicate(
                    rep_num, t(Alpha)[t(Id.Alpha)])) > 0) * (q.Alpha[, (1:n.Alpha)*2] - t(
                      replicate(rep_num, t(Alpha)[t(Id.Alpha)])) > 0) == 0))


med.list <- list(Mean = colMeans(mean.med), Sd = colMeans(sd.med),
                 Bias = colMeans(mean.med) - med,
                 Rms = sqrt(rowSums((t(mean.med) - med)^2) / rep_num),
                 p.CI95 = colMeans((q.med[ , (1:n.med)*2 - 1] - t(replicate(
                   rep_num, med)) > 0) * (q.med[, (1:n.med)*2] - t(
                     replicate(rep_num, med)) > 0) == 0))

est.list <- list(B = B.list, psi = psi.list, L.se = L.se.list,
                 psd = psd.list, gam = gam.list, lam0 = lam0.list, psw = psw.list, Alpha = Alpha.list, 
                 med = med.list)

names(est.list)[3] <- "Beta"
for (i in 1:length(est.list)) {
  write.csv(est.list[[i]], paste(names(est.list)[i],".csv",sep = ""))
}

time.total <- proc.time() - time.begin
save.image(paste(crep, '.Rdata', sep = ''))
cat('Task complete in', time.total[3], 'seconds. Results saved in lists and csv files.\n')




