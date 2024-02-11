\name{ginterpolating}
\alias{ginterpolating}
\title{Interpolation of a Graph Signal}
\description{
This function interpolates a graph signal utilizing the Laplacian matrix.
}
\usage{
ginterpolating(ad_mat, signal, vertices)
}
\arguments{
\item{ad_mat}{an weighted adjacency matrix.}
\item{signal}{a graph signal.}
\item{vertices}{specifies vertices for the observed signal. A signal on \code{vertices} and Laplacian matrix
  is used for interpolating a signal outside \code{vertices}.}
}
\details{
This function interpolates a graph signal utilizing the Laplacian matrix.
}
\value{ 
a signal with interpolated points.
}
\references{
Ortega, A., Frossard, P., Kovačević, J., Moura, J. M. F., and Vandergheynst, P. (2018).
Graph signal processing: overview, challenges, and applications. \emph{Proceedings of the IEEE 106}, 808--828.
\doi{https://doi.org/10.1109/JPROC.2018.2820126}

Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., and Vandergheynst, P. (2013). The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains. \emph{IEEE Signal Processing Magazine}, \bold{30(3)}, 83--98.
\doi{https://doi.org/10.1109/MSP.2012.2235192}

Tremblay, N., Borgnat, P., and Flandrin, P. (2014). Graph empirical mode decomposition. \emph{22nd European Signal Processing Conference (EUSIPCO)}, 2350--2354.

Zeng, J., Cheung, G., and Ortega, A. (2017). Bipartite approximation for graph wavelet signal decomposition. \emph{IEEE Transactions on Signal Processing}, \bold{65(20)}, 5466--5480.
\doi{https://doi.org/10.1109/TSP.2017.2733489}
}
\seealso{
\code{\link{gextrema}}, \code{\link{gsmoothing}}, \code{\link{sgemd}}.
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
menvelope <- (uenvelope + lenvelope) / 2

# display a graph signal
gplot(gsig, size=3, legend=FALSE) 

# display mean envelope
gplot(gsig, menvelope, size=3, legend=FALSE) 
}
\keyword{nonparametric}
