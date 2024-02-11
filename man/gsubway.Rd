\name{gsubway}
\alias{gsubway}
\docType{data}
\title{Seoul Subway Ridership Data}
\description{
This \pkg{igraph} graph object represents the number of subway passengers in Seoul, Korea for January 2021.
Consider the subway stations as vertices and railroads between the stations as the edges. 
The number of passengers for each station is regarded as the graph signal.
Seoul subway ridership data can be obtained at \code{https://www.seoulmetro.co.kr}.
}
\usage{
data(gsubway)
}
\format{
\pkg{igraph} graph object. The vertex attributes are longitude \code{x}, latitude \code{y}, and the number of passengers \code{z} for each station.
}
\examples{
data(gsubway)

# attributes
names(vertex_attr(gsubway)); names(edge_attr(gsubway)); names(graph_attr(gsubway))

# standardizing the graph signal
V(gsubway)$z <- c(scale(V(gsubway)$z))
 
# statistical graph empirical mode decomposition (SGEMD) with boundary treatment
out <- sgemd(gsubway, nimf=1, smoothing=TRUE, boundary=TRUE, connweight="graph")

# display of a signal, denoised signal by SGEMD
gplot(gsubway, size=3) 
gplot(gsubway, out$residue, size=3) 
}
\keyword{datasets}

