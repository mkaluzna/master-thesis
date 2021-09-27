library(matlab)
library(pheatmap)
library(dplyr)
library(igraph)
library(reshape2)
library(gplots)
library(ggplot2)
library(corrplot)
library(reticulate)
library(seriation)
library(ramify)
library(pracma)
library(stats)
library(rTensor)
library(RColorBrewer)
# library(dendextend)
library(aricode)
library(dbscan)

use_virtualenv("base")
source_python("Spinner/ProxG.py")
source_python("Spinner/Prox_h.py")
source_python("Spinner/prox_Fsvd.py")
source_python("Spinner/proxH_lasso.py")
source_python("Spinner/min_norm_estim.py")
source_python("Spinner/sol_options.py")
source_python("Spinner/spinner_lasso.py")
source_python("Spinner/spinner_nuclear.py")
source_python("Spinner/spinner_both.py")
source_python("Spinner/spinner.py")
source_python("Spinner/spinner_CV.py")


generate_AA_matrix <- function(n, p){
  AA = ramify::randn(p, p, n)
  AA = AA + aperm(AA, c(2, 1, 3))
  AA2 = Reshape(AA, p**2, n)
  AA2 = t(scale(t(AA2)))
  idxsDiag = as.logical(Reshape(diag(p), p**2, 1))
  AA2[idxsDiag,] = 0
  return(AA2)
}


generate_y_resp <- function(sigma, n, AA, B, X, beta){
  eps = sigma * randn(n, 1)
  y = array(t(AA) %*% as.vector(B) + X %*% beta + eps)
  return(y)
}


plot_heatmap <- function(M, title = "", paletteLength = 50, width = 5.5, height = 5.5,
                         col_names = FALSE, row_names = FALSE){
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(-max(M), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(M)/paletteLength, max(M), length.out = floor(paletteLength/2)))
  pheatmap(M, color = myColor, breaks = myBreaks, cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = col_names, show_colnames = row_names, cellwidth = width, cellheight = height,
           main = title, border_color = "NA")
}


plot_weighted_graph <- function(M){
  net = graph.adjacency(M, mode = "undirected", weighted = TRUE, diag = FALSE)
  plot.igraph(net, vertex.label=V(net)$name, layout=layout.fruchterman.reingold, 
              edge.color="red", edge.width=E(net)$weight, vertex.color="grey")
}


create_line_pattern_matrix <- function(n){
  M <- matrix(0, n, n)
  non_zeros <- sample(2:(n-1), 1)
  M[non_zeros,] <- runif(n)
  M[,non_zeros] <- M[non_zeros,]
  return(M)
}


create_bands_pattern_matrix <- function(n, bands_start = 2, bands_prop = 0.4){
  bands = c(bands_start:(floor(bands_prop * n) + bands_start))
  num_bands = length(bands)
  bMat <- matrix(runif(n * num_bands), n, num_bands, byrow=TRUE)
  M <- as.matrix(Matrix::bandSparse(n, n, k = bands, diag = bMat, symm=TRUE))
  return(M)
}


create_noise_matrix <- function(n, non_zero_prop = 0.5){
  M <- as.matrix(Matrix::rsparsematrix(n, n, non_zero_prop, symmetric = TRUE, 
                                       rand.x = function(n) runif(n)))
  return(M)
}


create_block_matrix <- function(p, dim_block, s, noise_nonzero_prop = 0.75, bands_nonzero_prop = 0.25){
  D1 <- create_noise_matrix(dim_block, non_zero_prop = noise_nonzero_prop)
  D1 <- D1 / max(eigen(D1)$values)
  D2 <- create_bands_pattern_matrix(dim_block, bands_prop = bands_nonzero_prop)
  D2 <- D2 / max(eigen(D2)$values)
  D3 <- create_line_pattern_matrix(dim_block)
  D3 <- D3 / max(eigen(D3)$values)
  dim_zeros <- p - 3 * dim_block
  D4 <- zeros(dim_zeros, dim_zeros)
  listMatrix <- list(D1, D2, D3, D4)
  B <- as.matrix(Matrix::.bdiag(listMatrix))
  B <- s * B
  rownames(B) <- 1:p
  colnames(B) <- 1:p
  return(B)
}


create_block_matrix_no_pattern <- function(p, dim_block, s, noise_nonzero_prop = 1){
  D1 <- create_noise_matrix(dim_block, noise_nonzero_prop)
  D1 <- D1 / max(eigen(D1)$values)
  D2 <- create_noise_matrix(dim_block - 1, noise_nonzero_prop)
  D2 <- D2 / max(eigen(D2)$values)
  D3 <- create_noise_matrix(dim_block - 2, noise_nonzero_prop)
  D3 <- D3 / max(eigen(D3)$values)
  dim_zeros <- p - (3 * dim_block - 3)
  D4 <- zeros(dim_zeros, dim_zeros)
  listMatrix <- list(D1, D2, D3, D4)
  B <- as.matrix(Matrix::.bdiag(listMatrix))
  B <- s * B
  rownames(B) <- 1:p
  colnames(B) <- 1:p
  return(B)
}


hierarchical_clustering <- function(B, k = 4){
  p = dim(B)[1]
  E = matrix(runif(p^2), p, p)
  E = (E + t(E)) / 2
  abs_B <- abs(B)
  dist_B <- max(abs_B) - abs_B
  dist_B <- dist_B + E
  diag(dist_B) <- 0
  dist_B <- as.dist(dist_B)
  hc_obj <- hclust(dist_B, method = "complete")
  hc_clusters <- cutree(hc_obj, k = k)
  names(hc_clusters) <- 1:length(hc_clusters)
  return(hc_clusters)
}

mutual_knn_graph <- function(X, nn){
  D <- as.matrix(as.dist(X))
  knn_mat <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  # find the nearest neighbors for each point
  for (i in 1: nrow(X)){
    neighbor_index <- order(D[i,])[2:(nn + 1)]
    knn_mat[i,][neighbor_index] <- 1 
  }
  
  # i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1 
  knn_mat <- knn_mat + t(knn_mat) # find mutual knn
  knn_mat[knn_mat == 2] = 1
  
  return(knn_mat)
}

spectral_clustering <- function(B,
                                n_eig = 3, # number of eigenvectors to keep
                                k = 4 # number of clusters
                                ){
  D <- max(abs(B)) - abs(B)  # 1. Compute dissimilarity matrix D
  diag(D) <- 0
  D_square <- D**2
  c <- max(D_square)
  S = exp(-(D_square)/c) # 2. Compute similarity matrix S
  nn = ceil(sqrt(dim(B)[1]))
  G = mutual_knn_graph(D, nn)  # 3. Compute similarity graph G
  S[G == 0] = 0
  diag(S) <- 0
  L = diag(colSums(S)) - S
  ei = eigen(L, symmetric = TRUE) 
  n = nrow(L) 
  sc <- ei$vectors[,(n - n_eig):(n - 1)] 
  B_kmeans <- kmeans(sc, k)
  sc_clusters <- B_kmeans$cluster
  names(sc_clusters) <- 1:length(sc_clusters)
  return(sc_clusters)
}


dbscan_clustering <- function(B, eps_val, minPts_val){
  abs_B <- abs(B)  # Compute dissimilarity matrix
  dist_B <- max(abs_B) - abs_B
  diag(dist_B) <- 0
  dbscan_clusters <- fpc::dbscan(dist_B, eps = eps_val, MinPts = minPts_val, method = "dist")$cluster
  # dbscan_clusters <- dbscan(B, eps=eps_val, minPts = minPts_val)$cluster
  names(dbscan_clusters) <- 1:length(dbscan_clusters)
  return(dbscan_clusters)
}


## df: Mutual Information
create_df_MI <- function(results_list, param, param_val){
  avg_MI <- lapply(results_list, colMeans)
  df_avg_MI <- data.frame(do.call(rbind, avg_MI))
  
  df_avg_MI[,param] <- param_val
  df_avg_MI_long <- df_avg_MI %>% 
    melt(id.vars = param, variable.name = "method") %>% 
    rename(MI = "value")
  
  return(df_avg_MI_long)
}

## df: Rand Index
create_df_RI <- function(results_list, param, param_val){
  avg_RI <- lapply(results_list, colMeans)
  df_avg_RI <- data.frame(do.call(rbind, avg_RI))
  
  df_avg_RI[,param] <- param_val
  df_avg_RI_long <- df_avg_RI %>% 
    melt(id.vars = param, variable.name = "method") %>% 
    rename(RI = "value")
  
  return(df_avg_RI_long)
}

## plot results
plot_metric <- function(df_long, param, metric){
  df_long  %>% 
    ggplot(aes(x = get(param), y = get(metric), group = method, color = method, shape = method)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    xlab(param) +
    ylab(metric)
}