\name{gftplot}
\alias{gftplot}
\title{Plot of the absolute values of the graph Fourier coefficients vs the eigenvalues}
\description{
This function displays the absolute values of the graph Fourier coefficients vs the eigenvalues.
}
\usage{
gftplot(graph, signal = NULL, K = NULL, size = 1, plot=TRUE) 
}
\arguments{
\item{graph}{an \pkg{igraph} graph object with vertex attributes of coordinates \code{x}, \code{y}, a signal \code{z}, and edge attribute of \code{weight}.}
\item{signal}{specifies a signal for the graph Fourier transform.
   When \code{signal=NULL}, a signal \code{z} of object \code{graph} is used.}
\item{K}{specifies the number of frequency components.}   
\item{size}{specifies point size.}
\item{plot}{specifies whether plot is displayed.}
}
\details{
This function displays the absolute values of the graph Fourier coefficients vs the eigenvalues for \code{signal}.
The red color denotes the nonnegative graph Fourier coefficients, and the blue color indicates the negative graph Fourier coefficients.
}
\value{ 
If \code{plot=TRUE}, plot of the absolute values of the graph Fourier coefficients vs the eigenvalues for \code{signal} over a graph using package \pkg{ggplot2}.
If \code{plot=FALSE}, a list with components:
\item{absgFCoeffs}{the absolute values of the graph Fourier coefficients in decreasing order.}
\item{eigenvalues}{the eigenvalues according to the absolute values of the graph Fourier coefficients.}
}
\seealso{
\code{\link{gsignal}}, \code{\link{gfdecomp}}.
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

# display a signal over graph
gplot(gsig, size=3) 

# display of absolute values of the graph Fourier coefficients vs the eigenvalues 
# for signal
gftplot(gsig)
\donttest{
gftplot(gsig, K=5, size=3)
out <- gftplot(gsig, K=5, plot=FALSE)
names(out)
}
## signal3
# high-frequency component
signal11 <- c(outer(sin(6*pi*x - 0.5*pi), sin(6*pi*y - 0.5*pi)))

# low-frequency component
signal22 <- c(outer(sin(2*pi*x - 0.5*pi), sin(2*pi*y - 0.5*pi)))

# composite signal
signal00 <- signal11 + signal22

# noisy signal
signal3 <- signal00 + rnorm(900, 0, sqrt(var(signal00) / 5)) 

# display signal3 over graph
gplot(gsig, signal=signal3, size=3) 
\donttest{
# display of absolute values of the graph Fourier coefficients vs the eigenvalues 
# for signal3
gftplot(gsig, signal=signal3)
gftplot(gsig, signal=signal3, K=10, size=2)
}
}
\keyword{nonparametric}
