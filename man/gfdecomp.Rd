\name{gfdecomp}
\alias{gfdecomp}
\title{Graph Fourier Decomposition}
\description{
This function performs the graph Fourier decomposition.
}
\usage{
gfdecomp(graph, K)
}
\arguments{
\item{graph}{an \pkg{igraph} graph object with vertex attributes of coordinates \code{x}, \code{y}, a signal \code{z}, and edge attribute of \code{weight}.}
\item{K}{specifies the number of frequency components.}
}
\details{
This function performs the graph Fourier decomposition.
}
\value{ 
\item{fc}{list of frequency components according to the frequencies with \code{fc[[1]]} the lowest-frequency component.}
\item{residue}{residue signal after extracting frequency components.}
}
\references{
Ortega, A., Frossard, P., Kovačević, J., Moura, J. M. F., and Vandergheynst, P. (2018).
Graph signal processing: overview, challenges, and applications. \emph{Proceedings of the IEEE 106}, 808--828.
\doi{https://doi.org/10.1109/JPROC.2018.2820126}

Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., and Vandergheynst, P. (2013). The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains. \emph{IEEE Signal Processing Magazine}, \bold{30(3)}, 83--98.
\doi{https://doi.org/10.1109/MSP.2012.2235192}
}
\seealso{
\code{\link{sgemd}}.
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
\donttest{
# display of absolute values of the graph Fourier coefficients vs the eigenvalues 
gftplot(gsig)
gftplot(gsig, K=5, size=3)
outgft <- gftplot(gsig, K=5, plot=FALSE)
outgft$eigenvalues
}
# graph Fourier decomposition 
out <- gfdecomp(gsig, K=4)
names(out)

# display of a signal, the extracted low- and high-frequency components by GFD
gplot(gsig, size=3) 
gplot(gsig, out$fc[[1]]+out$fc[[2]], size=3) 
gplot(gsig, out$fc[[3]]+out$fc[[4]], size=3) 
}
\keyword{nonparametric}
