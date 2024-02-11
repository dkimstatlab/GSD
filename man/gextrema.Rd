\name{gextrema}
\alias{gextrema}
\title{Finding Local Extrema of a Graph Signal}
\description{
This function finds the local extrema of a graph signal identifying the edge between neighboring vertices. 
}
\usage{
gextrema(ad_mat, signal)
}
\arguments{
\item{ad_mat}{an weighted adjacency matrix.}
\item{signal}{a graph signal.}
}
\details{
This function finds the local extrema of a graph signal identifying the edge between neighboring vertices. 
}
\value{ 
\item{maxima_list}{vertex index of local maxima of a \code{signal}.}
\item{minima_list}{vertex index of local minima of a \code{signal}.}
\item{n_extrema}{the number of local maxima and local minima of a \code{signal}.}
}
\references{
Tremblay, N., Borgnat, P., and Flandrin, P. (2014). Graph empirical mode decomposition. \emph{22nd European Signal Processing Conference (EUSIPCO)}, 2350--2354.
}
\seealso{
\code{\link{ginterpolating}}, \code{\link{gsmoothing}}, \code{\link{sgemd}}.
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

# local extrema
gextrema(A, signal)

# local extrema using graph object
gextrema(as_adjacency_matrix(gsig, attr="weight"), V(gsig)$z)
}
\keyword{nonparametric}
