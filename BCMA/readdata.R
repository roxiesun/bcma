#### read data from file

Z <- Zr <- array(scan(File.Z, skip = (crep - 1)*n, nlines = n, quiet = TRUE), dim = c(n, z))
S <- Sr <- array(scan(File.S, skip = (crep - 1)*n, nlines = n, quiet = TRUE), dim = c(n, s))
X <- Xr <- array(scan(File.X, skip = (crep - 1)*n, nlines = n, quiet = TRUE), dim = c(n, 2))
V <- Vr <- array(scan(File.V, skip = (crep - 1)*n, nlines = n, quiet = TRUE), dim = c(n, p))
W <- Wr <- array(scan(File.W, skip = (crep - 1)*n, nlines = n, quiet = TRUE), dim = c(n, w))
cat('the ', crep, 'th data read starts with:')
print(head(cbind(Zr, Sr, Xr, Vr, Wr)))
print(colMeans(head(cbind(Zr, Sr, Xr, Vr, Wr))))
cat('\r')



          