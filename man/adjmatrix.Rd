\name{adjmatrix}
\alias{adjmatrix}
\title{Weighted Adjacency Matrix}
\description{
This function produces weighted adjacency matrix by Gaussian kernel.
}
\usage{
adjmatrix(xy, method = c("dist", "neighbor"), alpha) 
}
\arguments{
\item{xy}{matrix or data.frame containing vertex coordinate \code{x}, \code{y}.}
\item{method}{When \code{method="dist"}, edge weights are calculated by Gaussian kernel for connecting vertices within distance \code{alpha}.
   When \code{method="neighbor"}, edge weights are calcaulated by Gaussian kernel for connecting \code{alpha} neighboring vertices.}
\item{alpha}{specifies distance between vertices when \code{method="dist"}, and the number of neighboring vertices when \code{method="neighbor"}.}
}   
\details{
This function produces a sparse weighted adjacency matrix by Gaussian kernel based on the distance between vertices.
}
\value{
a sparse weighted adjacency matrix
}
\references{
Zeng, J., Cheung, G., and Ortega, A. (2017). Bipartite approximation for graph wavelet signal decomposition. \emph{IEEE Transactions on Signal Processing}, \bold{65(20)}, 5466--5480.
\doi{https://doi.org/10.1109/TSP.2017.2733489}
}

\seealso{
\code{\link{gsignal}}, \code{\link{gplot}}.
}
\examples{
## define vertex coordinate
x <- y <- seq(0, 1, length=30)
xy <- expand.grid(x=x, y=y)

## weighted adjacency matrix by Gaussian kernel 
## for connecting vertices within distance 0.04
A1 <- adjmatrix(xy, method = "dist", 0.04) 

## weighted adjacency matrix by Gaussian kernel 
## for connecting seven neighboring vertices
A2 <- adjmatrix(xy, method="neighbor", 7)
}
\keyword{nonparametric}
