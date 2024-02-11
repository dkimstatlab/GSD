\name{gsmoothing}
\alias{gsmoothing}
\title{Smoothing a Graph Signal}
\description{
This function denoises a graph signal.
}
\usage{
gsmoothing(ad_mat, signal)
}
\arguments{
\item{ad_mat}{an weighted adjacency matrix.}
\item{signal}{a graph signal.}
}
\details{
This function denoises a graph signal utilizing the graph Fourier transform and empirical Bayes thresholding.
}
\value{ 
a denoised signal.
}
\references{
Ortega, A., Frossard, P., Kovačević, J., Moura, J. M. F., and Vandergheynst, P. (2018).
Graph signal processing: overview, challenges, and applications. \emph{Proceedings of the IEEE 106}, 808--828.
\doi{https://doi.org/10.1109/JPROC.2018.2820126}

Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., and Vandergheynst, P. (2013). The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains. \emph{IEEE Signal Processing Magazine}, \bold{30(3)}, 83--98.
\doi{https://doi.org/10.1109/MSP.2012.2235192}

Johnstone, I. and Silverman, B.~W. (2004).
Needles and straw in haystacks: empirical Bayes estimates of possibly sparse sequences. \emph{The Annals of Statistics}, \bold{32}, 594--1649.
\doi{https://doi.org/10.1214/009053604000000030}
}
\seealso{
\code{\link{gextrema}}, \code{\link{ginterpolating}}, \code{\link{sgemd}}.
}
\examples{
#### example : composite of two components having different frequencies

## define vertex coordinate
x <- y <- seq(0, 1, length=30)
xy <- expand.grid(x=x, y=y)

## weighted adjacency matrix by Gaussian kernel 
## for connecting vertices within distance 0.04
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

# local extrema using graph object
extremaout <- gextrema(as_adjacency_matrix(gsig, attr="weight"), V(gsig)$z)
maxima <- extremaout$maxima_list; minima <- extremaout$minima_list

# Interpolation of upper, lower and mean envelope
uenvelope <- ginterpolating(as_adjacency_matrix(gsig, attr="weight"), 
  V(gsig)$z, maxima)
lenvelope <- ginterpolating(as_adjacency_matrix(gsig, attr="weight"), 
  V(gsig)$z, minima) 

# smoothing upper, lower and mean envelope
suenvelope <- gsmoothing(A, uenvelope)
slenvelope <- gsmoothing(A, lenvelope) 

smenvelope <- (suenvelope + slenvelope) / 2

# display a graph signal
gplot(gsig, size=3, legend=FALSE) 

# display mean envelope
gplot(gsig, smenvelope, size=3, legend=FALSE) 
}
\keyword{nonparametric}
