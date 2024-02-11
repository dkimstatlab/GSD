\name{sgemd}
\alias{sgemd}
\title{Statistical Graph Empirical Mode Decomposition}
\description{
This function performs statistical graph empirical mode decomposition.
}
\usage{
sgemd(graph, nimf, smoothing = FALSE, smlevels = c(1),  
    boundary = FALSE, reflperc = 0.3, reflaver = FALSE, 
    connperc = 0.05, connweight = "boundary", 
    tol = 0.1^3, max.sift = 50, verbose = FALSE)
}
\arguments{
\item{graph}{an \pkg{igraph} graph object with vertex attributes of coordinates \code{x}, \code{y}, a signal \code{z}, and edge attribute of \code{weight}.}
\item{nimf}{specifies the maximum number of intrinsic mode functions (IMF).}
\item{smoothing}{specifies whether intrinsic mode functions are constructed by interpolating or smoothing envelopes. 
   When \code{smoothing =TRUE}, denoise envelopes utilizing the graph Fourier transform and empirical Bayes thresholding.}
\item{smlevels}{specifies which level of the IMF is obtained by smoothing other than interpolation.}
\item{boundary}{When \code{boundary=TRUE}, a given graph is reflected for boundary treatment.} 
\item{reflperc}{expand a graph by adding specified percentage of a graph at the boundary when \code{boundary=TRUE}.}
\item{reflaver}{specifies the method assigning signal to reflected vertices. 
   When \code{reflaver=TRUE}, the signal on reflected vertices is produced by averaging signals on neighboring vertices on a given graph.
   Otherwise, signal on reflected vertices is the same to the signal on a given graph.}
\item{connperc}{specifies percentage of a graph for connecting a given graph and reflected graph when \code{boundary=TRUE}.} 
\item{connweight}{specifies the method assigning the edge weights for connecting a given graph and reflected graph 
   when \code{boundary=TRUE}.
   The edge weights are calculated by Gaussian kernel considering the relative distance between vertices.
   When \code{connweight="graph"}, the relative distance is calculated based on the maximum distance of all the neighboring edges 
   of a given graph.
   When \code{connweight="boundary"}, the relative distance is calculated based on the maximum distance of the connected vertices 
   between a given graph and reflected graph.}
\item{tol}{tolerance for stopping rule of sifting.}
\item{max.sift}{the maximum number of sifting.} 
\item{verbose}{specifies whether sifting steps are displayed.}
}
\details{
This function performs statistical graph empirical mode decomposition utilizing extrema detection of a graph signal.
}
\value{ 
\item{imf}{list of IMF's according to the frequencies with \code{imf[[1]]} the highest-frequency IMF.}
\item{residue}{residue signal after extracting IMF's.}
\item{nimf}{the number of IMF's.} 
\item{n_extrema}{Each row specifies the number of local maxima and local minima of the remaining signal after extracting the i-th IMF. 
  The first row represents the number of local maxima and local minima of a given signal.}
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. C., Shih, H. H., Zheng, Q., Yen, N.- C., Tung, C. C., and Liu, H. H. (1998). The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis. \emph{Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences}, \bold{454(1971)}, 903–-995.
\doi{https://doi.org/10.1098/rspa.1998.0193}

Johnstone, I. and Silverman, B.~W. (2004).
Needles and straw in haystacks: empirical Bayes estimates of possibly sparse sequences. \emph{The Annals of Statistics}, \bold{32}, 594--1649.
\doi{https://doi.org/10.1214/009053604000000030}

Kim, D., Kim, K. O., and Oh, H.-S. (2012a). Extending the scope of empirical mode decomposition by smoothing. \emph{EURASIP Journal on Advances in Signal Processing}, \bold{2012}, 1--17.
\doi{https://doi.org/10.1186/1687-6180-2012-168}

Kim, D., Park, M., and Oh, H.-S. (2012b). Bidimensional statistical empirical mode decomposition. \emph{IEEE Signal Processing Letters}, \bold{19(4)}, 191--194.
\doi{https://doi.org/10.1109/LSP.2012.2186566}

Ortega, A., Frossard, P., Kovačević, J., Moura, J. M. F., and Vandergheynst, P. (2018).
Graph signal processing: overview, challenges, and applications. \emph{Proceedings of the IEEE 106}, 808--828.
\doi{https://doi.org/10.1109/JPROC.2018.2820126}

Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., and Vandergheynst, P. (2013). The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains. \emph{IEEE Signal Processing Magazine}, \bold{30(3)}, 83--98.
\doi{https://doi.org/10.1109/MSP.2012.2235192}

Tremblay, N., Borgnat, P., and Flandrin, P. (2014). Graph empirical mode decomposition. \emph{22nd European Signal Processing Conference (EUSIPCO)}, 2350--2354

Zeng, J., Cheung, G., and Ortega, A. (2017). Bipartite approximation for graph wavelet signal decomposition. \emph{IEEE Transactions on Signal Processing}, \bold{65(20)}, 5466--5480.
\doi{https://doi.org/10.1109/TSP.2017.2733489}
}
\seealso{
\code{\link{gextrema}}, \code{\link{gsmoothing}}, \code{\link{ginterpolating}}.
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

# graph empirical mode decomposition (GEMD) without boundary treatment
out1 <- sgemd(gsig, nimf=3, smoothing=FALSE, boundary=FALSE)

# denoised signal by GEMD 
dsignal1 <- out1$imf[[2]] + out1$imf[[3]] + out1$residue

\donttest{
# statistical graph empirical mode decomposition (SGEMD) with boundary treatment
out2 <- sgemd(gsig, nimf=3, smoothing=TRUE, boundary=TRUE)
names(out2)

# denoised signal by SGEMD 
dsignal2 <- out2$imf[[2]] + out2$imf[[3]] + out2$residue

# display of a signal, denoised signal, imf2, imf3 and residue by SGEMD
gplot(gsig, size=3) 
gplot(gsig, dsignal2, size=3) 
gplot(gsig, out2$imf[[2]], size=3) 
gplot(gsig, out2$imf[[3]], size=3) 
gplot(gsig, out2$residue, size=3) 
}
}
\keyword{nonparametric}
