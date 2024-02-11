# version 1.0.0 (Github)
gsignal <- function(vertex, edge, edgetype=c("matrix", "list")) {
  
  if (edgetype == "matrix") g <- graph_from_adjacency_matrix(edge, mode="undirected", weighted=TRUE)
  if (edgetype == "list") {
    g <- graph_from_edgelist(as.matrix(edge[, 1:2]), directed = FALSE)
    edge_attr(g, "weight") <- edge[, 3] 
  }
  vertex_attr(g, "x") <- vertex[, 1]
  vertex_attr(g, "y") <- vertex[, 2]
  vertex_attr(g, "z") <- vertex[, 3] 
  
  g
} 

adjmatrix <- function(xy, method=c("dist", "neighbor"), alpha) {
  
  # alpha : distance when method == "dist" or number of neighbor vertices when method == "neighbor"
  if (method == "dist" && alpha < 0) stop("When method == \"dist\", \"alpha\" must be positive.") 
  if (method == "neighbor" && alpha < 0) stop("When method == \"neighbor\", \"alpha\" must be positive interger.") 
  if (method == "neighbor") alpha <- as.integer(alpha) 
  
  n <- nrow(xy)
  dst_mat <- as.matrix(dist(xy, diag=TRUE, upper=TRUE))
  if (method == "dist") {
    zeta <- max(dst_mat[dst_mat <= alpha])
    ad_mat <- exp(-dst_mat^2/(2*zeta^2)) * (dst_mat > 0 & dst_mat <= alpha)  
  } else if (method == "neighbor") {
    ad_mat <- matrix(0, n, n)
    idx_connected_vertices <- t(apply(dst_mat, 1, order))[, 1:alpha+1]
    
    for (i in 1:n) 
      ad_mat[i, idx_connected_vertices[i,]] <- ad_mat[idx_connected_vertices[i,], i] <- dst_mat[i, idx_connected_vertices[i,]]
    
    maxad <- max(ad_mat)
    ad_mat <- exp(-ad_mat^2/(2*maxad^2)) * (ad_mat > 0 & ad_mat <= maxad)
  }
  
  Matrix(ad_mat, sparse=TRUE)
}

gplot <- function(graph, signal=NULL, size=1, limits=range(V(graph)$z), gpalette=NULL, legend=TRUE) {
  
  if (is.null(gpalette)) 
    gpalette <- palette(c('#00008D', '#002AFF', '#00D4FF', '#2AFFD4', '#FFFF00', '#FF8D00', '#FF0000'))
  
  if (is.null(signal)) 
    signal <- V(graph)$z
  
  ncolours = length(gpalette)
  bpoints = seq(limits[1], limits[2], length=length(gpalette))
  
  x <- V(graph)$x; y <- V(graph)$y 
  edgelist <- as_edgelist(graph, names=FALSE)
  
  x1 <- x[edgelist[,1]]
  y1 <- y[edgelist[,1]]  
  x2 <- x[edgelist[,2]]
  y2 <- y[edgelist[,2]]
  
  gplot <- ggplot(data = data.frame(x = x, y = y), aes(x, y)) + xlab("") + ylab("")
  gplot <- gplot + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = "gray", 
                                data = data.frame(x1 = x1, y1 = y1, x2 = x2, y2 = y2)) 
  gplot <- gplot + geom_point(size = size, aes(colour = signal), show.legend=legend) + labs(colour = "")
  gplot <- gplot + scale_color_gradientn(limits = limits, colours = gpalette, 
                                         breaks = bpoints, labels = format(bpoints, digits=2)) 
  gplot <- gplot + 
    theme(
      legend.margin = margin(0.0, 0.0, 0.0, 0.0),
      plot.margin = unit(rep(0.02, 4), "cm"),
      axis.title=element_blank()
    )
  
  gplot
}

gextrema <- function(ad_mat, signal) {
  maxima_list <- minima_list <- NULL
  for (vertex in 1:length(signal)) {
    connected_vertices <- which(ad_mat[vertex,] != 0)
    
    if (all(signal[vertex] >= signal[connected_vertices])) {
      maxima_list <- c(maxima_list, vertex)
    } else if (all(signal[vertex] <= signal[connected_vertices])) {
      minima_list <- c(minima_list, vertex)
    }
  }
  list(maxima_list=maxima_list, minima_list=minima_list, n_extrema=c(length(maxima_list), length(minima_list)))
}

ginterpolating <- function(ad_mat, signal, vertices) {
  
  n <- length(signal) 
  Laplacian <- diag(rowSums(ad_mat)) - ad_mat 
  s.known <- vertices
  s.unknown <- setdiff(1:n, s.known)
  
  RT <- Laplacian[s.unknown, s.known]
  Lu <- Laplacian[s.unknown, s.unknown]
  sb <- signal[s.known] 

  su <- solve(Lu, -as.matrix(RT %*% sb))

  signal[s.unknown] <- su
  
  signal
}

gsmoothing <- function(ad_mat, signal) {
  
  L <- diag(rowSums(ad_mat)) - ad_mat
  eigenvectors <- eigen(L)$vectors
  
  GFT <- t(eigenvectors) %*% signal 
  GFT.Thresh <- ebayesthresh(GFT, verbose = TRUE, threshrule = 'soft') 
  
  GFT.inverse.smoothing <- eigenvectors %*% GFT.Thresh$muhat 
  
  c(GFT.inverse.smoothing)
}

sgemd <- function(graph, nimf, smoothing=FALSE, smlevels = c(1), boundary=FALSE, reflperc=0.3, 
                      reflaver=FALSE, connperc=0.05, connweight="boundary", tol=0.1^3, max.sift=50, verbose=FALSE) {
  
  tol_n_extrema <- 4; n_extrema <- NULL
  n_vertices <- vcount(graph)
  ad_mat <- as_adjacency_matrix(graph, attr="weight")
  
  if (boundary) {
    # define the reflected graph
    xy <- data.frame(x=V(graph)$x, y=V(graph)$y)

    connperc <- min(connperc, reflperc)
    
    x_refl <- quantile(xy$x, c(reflperc, 1-reflperc)) 
    y_refl <- quantile(xy$y, c(reflperc, 1-reflperc)) 
    
    x_conn <- quantile(xy$x, c(connperc, 1-connperc)) 
    y_conn <- quantile(xy$y, c(connperc, 1-connperc)) 
    x_conn <- c(max(x_conn[1], sort(unique(xy$x))[2]), min(x_conn[2], sort(unique(xy$x), decreasing=TRUE)[2]))
    y_conn <- c(max(y_conn[1], sort(unique(xy$y))[2]), min(y_conn[2], sort(unique(xy$y), decreasing=TRUE)[2]))
    
    x_end <- c(min(xy$x), max(xy$x)); y_end <- c(min(xy$y), max(xy$y)) 
    
    x_delta <- abs(x_end - x_conn)  
    y_delta <- abs(y_end - y_conn)  
      
    # indices of vertices to be reflected
    left_refl_idx <- which(xy$x > x_end[1] & xy$x <= x_refl[1]); n_left <- length(left_refl_idx)
    right_refl_idx <- which(xy$x < x_end[2] & xy$x >= x_refl[2]); n_right <- length(right_refl_idx)
    lower_refl_idx <- which(xy$y > y_end[1] & xy$y <= y_refl[1]); n_lower <- length(lower_refl_idx)
    upper_refl_idx <- which(xy$y < y_end[2] & xy$y >= y_refl[2]); n_upper <- length(upper_refl_idx)
    
    left_upper_refl_idx <- intersect(left_refl_idx, upper_refl_idx); n_left_upper <- length(left_upper_refl_idx)
    right_upper_refl_idx <- intersect(right_refl_idx, upper_refl_idx); n_right_upper <- length(right_upper_refl_idx)
    left_lower_refl_idx <- intersect(left_refl_idx, lower_refl_idx); n_left_lower <- length(left_lower_refl_idx)
    right_lower_refl_idx <- intersect(right_refl_idx, lower_refl_idx); n_right_lower <- length(right_lower_refl_idx)
    
    n_added <- cumsum(c(n_vertices, n_left, n_right, n_lower, n_upper, n_left_upper, n_right_upper, n_left_lower))
    
    refl_idx <- c(left_refl_idx, right_refl_idx, lower_refl_idx, upper_refl_idx, 
                  left_upper_refl_idx, right_upper_refl_idx, left_lower_refl_idx, right_lower_refl_idx)
    n_refl <- length(refl_idx)

    # define the coordinates of the reflected vertices and find vertices for connecting a given graph and the reflected parts 
    # conn_~ denotes the vertices which need to be connected with the neighbor reflected sectors (original graph, left, right, ...)
    
    left_refl <- right_refl <- lower_refl <- upper_refl <- 
      left_upper_refl <- right_upper_refl <- left_lower_refl <- right_lower_refl <- NULL
    conn_left_refl_idx <- conn_right_refl_idx <- conn_lower_refl_idx <- conn_upper_refl_idx <- 
      conn_left_upper_refl_idx <- conn_right_upper_refl_idx <- conn_left_lower_refl_idx <- conn_right_lower_refl_idx <- NULL
 
    tmpx1_idx <- which(xy$x <= x_conn[1])
    tmpx2_idx <- which(xy$x >= x_conn[2])
    tmpy1_idx <- which(xy$y <= y_conn[1])
    tmpy2_idx <- which(xy$y >= y_conn[2])
    
    conn_center_idx <- unique(c(tmpx1_idx, tmpx2_idx, tmpy1_idx, tmpy2_idx))
    n_center <- length(conn_center_idx)
    
    if (n_center == 0) 
      conn_center_idx <- NULL
    
    if (n_left > 0) {
      left_refl <- xy[left_refl_idx,]
      left_refl$x <- left_refl$x * (-1) + 2 * x_end[1]
      conn_left_refl_idx <- which(!is.na(match(left_refl_idx, conn_center_idx))) + n_added[1]
    } 
    if (n_right > 0) {
      right_refl <- xy[right_refl_idx,]
      right_refl$x <- right_refl$x * (-1) + 2 * x_end[2]
      conn_right_refl_idx <- which(!is.na(match(right_refl_idx, conn_center_idx))) + n_added[2]
    }
    if (n_lower > 0) {
      lower_refl <- xy[lower_refl_idx,]
      lower_refl$y <- lower_refl$y * (-1) + 2 * y_end[1]
      conn_lower_refl_idx <- which(!is.na(match(lower_refl_idx, conn_center_idx))) + n_added[3]
    }
    if (n_upper > 0) {
      upper_refl <- xy[upper_refl_idx,]
      upper_refl$y <- upper_refl$y * (-1) + 2 * y_end[2]
      conn_upper_refl_idx <- which(!is.na(match(upper_refl_idx, conn_center_idx))) + n_added[4]
    }
    if (n_left_upper > 0) {
      left_upper_refl <- xy[left_upper_refl_idx,]
      left_upper_refl$x <- left_upper_refl$x * (-1) + 2 * x_end[1]
      left_upper_refl$y <- left_upper_refl$y * (-1) + 2 * y_end[2]
      conn_left_upper_refl_idx <- which(!is.na(match(left_upper_refl_idx, conn_center_idx))) + n_added[5]
    }
    if (n_right_upper > 0) {
      right_upper_refl <- xy[right_upper_refl_idx,]
      right_upper_refl$x <- right_upper_refl$x * (-1) + 2 * x_end[2]
      right_upper_refl$y <- right_upper_refl$y * (-1) + 2 * y_end[2]
      conn_right_upper_refl_idx <- which(!is.na(match(right_upper_refl_idx, conn_center_idx))) + n_added[6]
    }
    if (n_left_lower > 0) {
      left_lower_refl <- xy[left_lower_refl_idx,]
      left_lower_refl$x <- left_lower_refl$x * (-1) + 2 * x_end[1]
      left_lower_refl$y <- left_lower_refl$y * (-1) + 2 * y_end[1]
      conn_left_lower_refl_idx <- which(!is.na(match(left_lower_refl_idx, conn_center_idx))) + n_added[7]
    } 
    if (n_right_lower > 0) {
      right_lower_refl <- xy[right_lower_refl_idx,]
      right_lower_refl$x <- right_lower_refl$x * (-1) + 2 * x_end[2]
      right_lower_refl$y <- right_lower_refl$y * (-1) + 2 * y_end[1]
      conn_right_lower_refl_idx <- which(!is.na(match(right_lower_refl_idx, conn_center_idx))) + n_added[8]
    }

    xy_refl <- rbind(xy, left_refl, right_refl, lower_refl, upper_refl, 
                     left_upper_refl, right_upper_refl, left_lower_refl, right_lower_refl)

    # define expanded adjacency matrix
    n_ad_mat_refl <- n_vertices + n_refl
    ad_mat_refl <- Matrix(0, nrow = n_ad_mat_refl, ncol = n_ad_mat_refl, sparse = TRUE)
    ad_mat_refl[1:n_vertices, 1:n_vertices] <- ad_mat
    
    # Keep the connectivity of original graph within each sector    
    if (n_left > 0) ad_mat_refl[(n_added[1]+1):n_added[2], (n_added[1]+1):n_added[2]] <- ad_mat[left_refl_idx, left_refl_idx]                        
    if (n_right > 0) ad_mat_refl[(n_added[2]+1):n_added[3], (n_added[2]+1):n_added[3]] <- ad_mat[right_refl_idx, right_refl_idx]
    if (n_lower > 0) ad_mat_refl[(n_added[3]+1):n_added[4], (n_added[3]+1):n_added[4]] <- ad_mat[lower_refl_idx, lower_refl_idx]
    if (n_upper > 0) ad_mat_refl[(n_added[4]+1):n_added[5], (n_added[4]+1):n_added[5]] <- ad_mat[upper_refl_idx, upper_refl_idx]
    if (n_left_upper > 0) ad_mat_refl[(n_added[5]+1):n_added[6], (n_added[5]+1):n_added[6]] <- ad_mat[left_upper_refl_idx, left_upper_refl_idx]
    if (n_right_upper > 0) ad_mat_refl[(n_added[6]+1):n_added[7], (n_added[6]+1):n_added[7]] <- ad_mat[right_upper_refl_idx, right_upper_refl_idx]
    if (n_left_lower > 0) ad_mat_refl[(n_added[7]+1):n_added[8], (n_added[7]+1):n_added[8]] <- ad_mat[left_lower_refl_idx, left_lower_refl_idx]
    if (n_right_lower > 0) ad_mat_refl[(n_added[8]+1):n_ad_mat_refl, (n_added[8]+1):n_ad_mat_refl] <- ad_mat[right_lower_refl_idx, right_lower_refl_idx]
  
    # define dist_max for calculating edge weights between a given graph and the reflected parts (conn_~ sets)
    if (connweight == "boundary") {
      dst_max <- max(x_delta, y_delta) 
    } else if (connweight == "graph")
      dst_max <- max(as.matrix(dist(xy, diag=TRUE, upper=TRUE))[as.logical(ad_mat != 0)]) 
    
    # connectivity between center and left
    tmpidx1 <- conn_center_idx
    tmpidx2 <- conn_left_refl_idx
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= x_delta[1] & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }   
    # connectivity between center and right
    tmpidx1 <- conn_center_idx
    tmpidx2 <- conn_right_refl_idx
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= x_delta[2] & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }    
    # connectivity between center and lower
    tmpidx1 <- conn_center_idx
    tmpidx2 <- conn_lower_refl_idx
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= y_delta[1] & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }       
    # connectivity between center and upper
    tmpidx1 <- conn_center_idx
    tmpidx2 <- conn_upper_refl_idx
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= y_delta[2] & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }      
    
    # connectivity between left upper and left, upper, center
    tmpidx1 <- conn_left_upper_refl_idx
    tmpidx2 <- c(conn_center_idx, conn_left_refl_idx, conn_upper_refl_idx)
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= max(x_delta[1], y_delta[2]) & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }
    
    # connectivity between right upper and right, upper, center
    tmpidx1 <- conn_right_upper_refl_idx
    tmpidx2 <- c(conn_center_idx, conn_right_refl_idx, conn_upper_refl_idx)
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= max(x_delta[2], y_delta[2]) & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }
    
    # connectivity between left lower and left, lower, center
    tmpidx1 <- conn_left_lower_refl_idx
    tmpidx2 <- c(conn_center_idx, conn_left_refl_idx, conn_lower_refl_idx)
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= max(x_delta[1], y_delta[1]) & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }
    
    # connectivity between right lower and right, lower, center
    tmpidx1 <- conn_right_lower_refl_idx
    tmpidx2 <- c(conn_center_idx, conn_right_refl_idx, conn_lower_refl_idx)
    if (length(tmpidx1) != 0 & length(tmpidx2) != 0) {
      tmpad <- matrix(0, nrow=length(tmpidx1), ncol=length(tmpidx2))
      tmpad <- outer(1:length(tmpidx1), 1:length(tmpidx2), 
                     FUN = Vectorize(function(i,j) dist(rbind(xy_refl[tmpidx1,][i,], xy_refl[tmpidx2,][j,]))))
      tmpad <- exp(-tmpad^2 / (2*dst_max^2)) * (tmpad <= max(x_delta[2], y_delta[1]) & tmpad != 0) 
      ad_mat_refl[tmpidx1, tmpidx2] <- tmpad
    }
    
    # Make the adjacency matrix symmetric
    tmpidx <- which(ad_mat_refl > t(ad_mat_refl), arr.ind=TRUE)
    ad_mat_refl[cbind(tmpidx[,2], tmpidx[,1])] <- ad_mat_refl[tmpidx]
    
    graph_refl <- graph_from_adjacency_matrix(ad_mat_refl, mode="undirected", weighted=TRUE)
    
    # obtain the connected graph and adjacency matrix
    graph_refl_connectivity <- components(graph_refl)
    graph_refl_mainindex <- which(graph_refl_connectivity$membership == which.max(graph_refl_connectivity$csize))
    ad_mat_refl <- ad_mat_refl[graph_refl_mainindex, graph_refl_mainindex]
    
    if (reflaver) {
      weights <- neighbor_vertices <- list()
      for (i in 1:n_refl) {
        neighbor_vertices[[i]] <- as.numeric(which(ad_mat[refl_idx[i], ] != 0))
        constants <- ad_mat[refl_idx[i], neighbor_vertices[[i]]]
        weights[[i]] <- constants / sum(constants)
      } 
    }
  } else 
    ad_mat_refl <- ad_mat

  imf <- list()
  residue <- V(graph)$z #signal
  
  for (i in 1:nimf) {
    input <- residue
    
    extrema <- gextrema(ad_mat, input)
    n_extrema <- rbind(n_extrema, extrema$n_extrema)
    if (extrema$n_extrema[1] < tol_n_extrema || extrema$n_extrema[2] < tol_n_extrema) { 
      residue <- input
      nimf <- i - 1
      warning("The number of extrema of remaining signal is not sufficient for interpolation to extract imf ", i, ".", "\n",
              "Thus resulting number of imf are ", nimf, ".")
      break  
    } 
    
    j <- 1
    repeat {
      if (verbose) print(paste0('imf ', i, ' : sifting = ', j))
      if (boundary) {
        if (reflaver) {
          # assign the remaining signal to reflected vertices
          input.refl <- NULL
          for (k in 1:n_refl) { 
            weighted_average <- sum(weights[[k]] * input[neighbor_vertices[[k]]]) 
            input.refl <- c(input.refl, weighted_average)
          }        
        } else 
          input.refl <- input[refl_idx]        
          input.refl <- c(input, input.refl)[graph_refl_mainindex]
      } else 
        input.refl <- input
      
      # extrema of the expanded graph signal by reflection
      extrema <- gextrema(ad_mat_refl, input.refl)
      maxima <- extrema$maxima_list; minima <- extrema$minima_list
      
      # Interpolation on expanded graph signal and Get signal of original graph
      uenvelope <- ginterpolating(ad_mat_refl, input.refl, maxima)[1:n_vertices]
      lenvelope <- ginterpolating(ad_mat_refl, input.refl, minima)[1:n_vertices]        
      if (any(i == smlevels) && smoothing) {
        # smoothing upper and lower envelopes
        uenvelope <- gsmoothing(ad_mat, uenvelope)
        lenvelope <- gsmoothing(ad_mat, lenvelope)
      }
      
      menvelope <- (uenvelope + lenvelope)  / 2 
      
      if (menvelope %*% menvelope < (input %*% input) * tol || j >= max.sift) 
        break 
      
      input <- input - menvelope
      
      j <- j + 1
    }
    
    imf[[i]] <- input 
    residue <- residue - input # remaining signal after extracting an IMF
    nimf <- i
  }
  #residue <- residue 
  
  if (extrema$n_extrema[1] >= tol_n_extrema && extrema$n_extrema[2] >= tol_n_extrema) { 
    extrema <- gextrema(ad_mat, residue)
    n_extrema <- rbind(n_extrema, extrema$n_extrema)
  }
  
  list(imf=imf, residue=residue, nimf=nimf, n_extrema=n_extrema)
}

gfdecomp <- function(graph, K) {
  
  signal <- V(graph)$z
  #n_vertices <- length(signal)
  ad_mat <- as_adjacency_matrix(graph, attr="weight")
  
  Laplacian_matrix <- diag(rowSums(ad_mat)) - ad_mat
  eigen.Laplacian <- eigen(Laplacian_matrix)
  eigenvalues <- eigen.Laplacian$values
  eigenvectors <- eigen.Laplacian$vectors
  
  GFT.signal <- t(eigenvectors) %*% signal
  sqrt_periodogram <- abs(GFT.signal)
  
  # Obtain the eigenvalues and eigenvectors corresponding to the K largest values of sqrt_periodogram
  idx.significant_component <- order(sqrt_periodogram, decreasing=TRUE)[1:K] 
  idx.significant_component <- idx.significant_component[order(eigenvalues[idx.significant_component])]
  
  fc <- list()
  # Extract the components from low frequency to high frequency
  for (i in 1:K)
    fc[[i]] <- GFT.signal[idx.significant_component[i]] * eigenvectors[, idx.significant_component[i]]
  
  residue <- signal - Reduce("+", fc)
  
  list(fc=fc, residue=residue)
}

gftplot <- function(graph, signal=NULL, K=NULL, size=1, plot=TRUE) {
  
  if (is.null(signal)) signal <- V(graph)$z
  if (is.null(K)) K <- length(signal)
  ad_mat <- as_adjacency_matrix(graph, attr="weight")
  
  Laplacian_matrix <- diag(rowSums(ad_mat)) - ad_mat
  eigen.Laplacian <- eigen(Laplacian_matrix)
  eigenvalues <- eigen.Laplacian$values
  eigenvectors <- eigen.Laplacian$vectors
  
  GFT.signal <- t(eigenvectors) %*% signal
  sqrt_periodogram <- abs(GFT.signal)
  
  # Obtain the eigenvalues and eigenvectors corresponding to the values of sqrt_periodogram
  idx.significant_component <- order(sqrt_periodogram, decreasing=TRUE)[1:K]
  x <- sqrt_periodogram[idx.significant_component]
  y <- eigenvalues[idx.significant_component]
  
  GFTcolor <- rep("red", K)
  GFTcolor[GFT.signal[idx.significant_component] < 0] <- "blue"
  
  gplot <- ggplot(data = data.frame(x = x, y = y), aes(x, y)) 
  gplot <- gplot + xlab("|graph Fourier Coefficients|") + ylab("eigenvalues")
  gplot <- gplot + geom_point(size = size, color=GFTcolor)
  
  if (plot) gplot
  else list(absgFCoeffs=x, eigenvalues=y)
}
