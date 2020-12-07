#####################
##defined constants##
rep_num <- 100
n <- 500
z <- 4
s <- 2
c <- c(3.5, 1, 1.3)
p <- 3
q <- 1
w <- 2
G <- 5 

n.mcmc = 15000
n.kept = 10000
n.burn = n.mcmc - n.kept
n.thin = 1
sig.mh = 3
sig.mh.g = 0.5

PI = 0   #for 1 latent mediator case

#File name
File.Z <- 'Z.txt'
File.S <- 'S.txt'
File.X <- 'X.txt'
File.V <- 'V.txt'
File.W <- 'W.txt'


####### total number of acceptance for each variable in MH algorithm

