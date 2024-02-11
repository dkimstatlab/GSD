\name{gsignal}
\alias{gsignal}
\title{Graph Object with a Signal}
\description{
This function constructs an \pkg{igraph} graph object with several vertex and edge attributes.
}
\usage{
gsignal(vertex, edge, edgetype = c("matrix", "list"))
}
\arguments{
\item{vertex}{matrix or data.frame containing vertex information.
   The first two columns are vertex coordinate \code{x}, \code{y}, and the third column is a signal \code{z} on each vertex.}
\item{edge}{When \code{edgetype="matrix"}, a square weighted adjacency matrix. 
   This can be a sparse matrix created with \pkg{Matrix} package.
   When \code{edgetype="list"}, matrix or data.frame of edge list with three columns.
   The first two columns are edge lists and the third column is an edge weight.}
\item{edgetype}{edges and weights information are provided by \code{"matrix"} or \code{"list"}.}
}
\details{
This function constructs an \pkg{igraph} graph object with \code{vertex} and \code{edge} attributes.
}
\value{ 
an \pkg{igraph} graph object. The vertex attributes are vertex coordinate \code{x}, \code{y}, and a signal \code{z} on each vertex.
The edge attribute is \code{weight}.
These vertex and edge attributes can be identified by \code{names(vertex_attr())} and \code{names(edge_attr())},
and are accessible by \code{V()}, \code{E()}, \code{as_edgelist()}, \code{as_adjacency_matrix()} or other \pkg{igraph} functions.
}
\seealso{
\code{\link{adjmatrix}}, \code{\link{gplot}}.
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

# vertex and edge attribute
names(vertex_attr(gsig)); names(edge_attr(gsig))

# edge list
# as_edgelist(gsig, name=FALSE)

# weighted adjacency matrix
# as_adjacency_matrix(gsig, attr="weight")

# display a noisy graph signal 
gplot(gsig, size=3)

# display a composite graph signal 
gplot(gsig, signal0, size=3)

# display high-frequency component
gplot(gsig, signal1, size=3)

# display low-frequency component
gplot(gsig, signal2, size=3)
}
\keyword{nonparametric}
