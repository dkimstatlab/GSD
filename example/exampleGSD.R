library(GSD)

## define vertex coordinate
x <- y <- seq(0, 1, length=30); xy <- expand.grid(x=x, y=y)

## weighted adjacency matrix by Gaussian kernel 
## for connecting vertices within a distance 0.04
A <- adjmatrix(xy, method = "dist", 0.04) 

## signal
# high-frequency component
signal1 <- rep(sin(12.5*pi*x - 1.25*pi), 30)
# low-frequency component
signal2 <- rep(sin(5*pi*x - 0.5*pi), 30)
# composite signal
signal0 <- signal1 + signal2
# noisy signal with SNR(signal-to-noise ratio)=5
signal <- signal0 + rnorm(900, 0, sqrt(var(signal0) / 5)) 

# graph with signal
gsig <- gsignal(vertex = cbind(xy, signal), edge = A, edgetype = "matrix")

# display a noisy graph signal 
gplot(gsig, size=3)
# display a composite graph signal 
gplot(gsig, signal0, size=3)
# display high-frequency component
gplot(gsig, signal1, size=3)
# display low-frequency component
gplot(gsig, signal2, size=3)

# graph empirical mode decomposition (GEMD) without boundary treatment
out1 <- sgemd(gsig, nimf=3, smoothing=FALSE, boundary=FALSE)

# display of the decomposed high-frequency component and 
# low-frequency component by GEMD
gplot(gsig, out1$imf[[2]], size=3, legend=FALSE) 
gplot(gsig, out1$imf[[3]], size=3, legend=FALSE) 

# statistical graph empirical mode decomposition (SGEMD) with boundary treatment
out2 <- sgemd(gsig, nimf=3, smoothing=TRUE, smlevels=1:3, boundary=TRUE)

# display of the decomposed high-frequency component and 
# low-frequency component by SGEMD
gplot(gsig, out2$imf[[2]], size=3, legend=FALSE) 
gplot(gsig, out2$imf[[3]], size=3, legend=FALSE) 

# display of the absolute values of the graph Fourier coefficients vs the eigenvalues 
gftplot(gsig)
outgft <- gftplot(gsig, K=5, plot=FALSE)
outgft$eigenvalues

# graph Fourier decomposition 
out3 <- gfdecomp(gsig, K=4)

# display of the decomposed high-frequency component and 
# low-frequency component by GFD
gplot(gsig, out3$fc[[3]]+out3$fc[[4]], size=3, legend=FALSE) 
gplot(gsig, out3$fc[[1]]+out3$fc[[2]], size=3, legend=FALSE) 

#Seoul subway ridership data
data(gsubway)

# standardizing the graph signal
V(gsubway)$z <- c(scale(V(gsubway)$z))

# GEMD without boundary treatment
out1 <- sgemd(gsubway, nimf=1, smoothing=FALSE, boundary=FALSE)

# SGEMD with boundary treatment
out2 <- sgemd(gsubway, nimf=1, smoothing=TRUE, boundary=TRUE, connweight="graph")

# display of a standardized signal
limits <- range(c(V(gsubway)$z, out1$imf[[1]], out2$imf[[1]]))
gplot(gsubway, size=3, limits=limits) 
# display of the first IMF by GEMD
gplot(gsubway, out1$imf[[1]], size=3, limits=limits) 
# display of the first IMF by SGEMD
gplot(gsubway, out2$imf[[1]], size=3, limits=limits) 
