\name{gplot}
\alias{gplot}
\title{Plot of a Graph Signal}
\description{
This function displays a signal on a graph using a color palette.
}
\usage{
gplot(graph, signal = NULL, size = 1, limits = range(V(graph)$z), 
    gpalette = NULL, legend = TRUE)
}
\arguments{
\item{graph}{an \pkg{igraph} graph object with vertex attributes of coordinates \code{x}, \code{y}, a signal \code{z}, and edge attribute of \code{weight}.}
\item{signal}{specifies a signal to be displayed over object \code{graph} using a color palette.
   When \code{signal=NULL}, a signal \code{z} of object \code{graph} is used.}
\item{size}{specifies point size of vertex.}
\item{limits}{specifies color palette limits.}
\item{gpalette}{specifies a character vector of color palette. When \code{gpalette=NULL}, \code{c('#00008D', '#002AFF', '#00D4FF', '#2AFFD4', 
  '#FFFF00', '#FF8D00', '#FF0000')} is used as default palette.}
\item{legend}{if \code{legend=FALSE}, the legend is not included.}
}
\details{
This function displays a signal on a graph using a color palette.
}
\value{ 
plot of a signal over a graph using package \pkg{ggplot2}.
}
\seealso{
\code{\link{gsignal}}, \code{\link{adjmatrix}}.
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

# display a signal over graph with legend
gplot(gsig, size=3, legend=TRUE) 

# display a signal over graph without legend
gplotout <- gplot(gsig, size=3, legend=FALSE) 
gplotout

# adding labels using ggplot2 package
gplotout + 
  theme(axis.title=element_text(), 
    plot.title=element_text(hjust = 0.5, vjust = 0)) + 
    labs(x="x", y="y", title="visualization of a composite signal")

# deleting axis title, text and ticks using ggplot2 package
gplotout + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

# display high-frequency component
gplot(gsig, signal1, size=3, legend=FALSE) 

# display low-frequency component
gplot(gsig, signal2, size=3, legend=FALSE) 
}
\keyword{nonparametric}
