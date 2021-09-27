source("final_functions.R")


### Fixed parameters
sigma = 0.1
m = 1
dim_block = 10
noise_nonzero_prop = 0.7
bands_nonzero_prop = 0.25


# Pattern matrix ----------------------------------------------------------
### Simulations for n
# Setting parameters
iter = 1
n_param = c(250, 300, 350, 400, 450)
p = 60
s = 25
W = ones(p, p) - diag(p)
beta <- runif(m)
B <- create_block_matrix(p, dim_block, s, noise_nonzero_prop = noise_nonzero_prop,
                         bands_nonzero_prop = bands_nonzero_prop)
real_clusters <- c(rep(1, dim_block), rep(2, dim_block), rep(3, dim_block), rep(4, p - 3 * dim_block))
names(real_clusters) <- 1:length(real_clusters)
plot_heatmap(B, "True signal")
sim_results_n <- list()
for(n in n_param){
  param <- paste0("n_", n)
  sim_results_n[[param]] = list(B = list(),
                            B_hat = list(),
                            real_clusters = list()
  )
}
for(n in n_param){
  X <- matrix(1, n, m)
  param <- paste0("n_", n)
  for(i in 1:iter){
    AA <- generate_AA_matrix(n, p)
    AA_vec <- array(as.vector(AA))
    y <- generate_y_resp(sigma, n, AA, B, X, beta)
    spinner_fit = calculate_CV(y = y, 
                               AA = AA_vec, 
                               W = W, 
                               X = X)
    B_hat = spinner_fit$B
    rownames(B_hat) <- 1:(dim(B_hat)[1])
    colnames(B_hat) <- 1:(dim(B_hat)[2])
    sim_results_n[[param]]$B[[(length(sim_results_n[[param]]$B) + 1)]] <- B
    sim_results_n[[param]]$B_hat[[(length(sim_results_n[[param]]$B_hat) + 1)]] <- B_hat
    sim_results_n[[param]]$real_clusters[[(length(sim_results_n[[param]]$real_clusters) + 1)]] <- real_clusters
    plot_heatmap(B_hat, "Estimated B")
  }
}


### Simulations for p
# Setting parameters
iter = 1
n = 350
p_param = c(40, 60, 80, 100)
s = 25
X <- matrix(1, n, m)
beta <- runif(m)
# plot_heatmap(B, "True signal")
sim_results_p <- list()
for(p in p_param){
  param <- paste0("p_", p)
  sim_results_p[[param]] = list(B = list(),
                          B_hat = list(),
                          real_clusters = list()
  )
}
for(p in p_param){
  W = ones(p, p) - diag(p)
  B <- create_block_matrix(p, dim_block, s, noise_nonzero_prop = noise_nonzero_prop,
                           bands_nonzero_prop = bands_nonzero_prop)
  plot_heatmap(B, "True signal")
  real_clusters <- c(rep(1, dim_block), rep(2, dim_block), rep(3, dim_block), rep(4, p - 3 * dim_block))
  names(real_clusters) <- 1:length(real_clusters)
  param <- paste0("p_", p)
  for(i in 1:iter){
    AA <- generate_AA_matrix(n, p)
    AA_vec <- array(as.vector(AA))
    y <- generate_y_resp(sigma, n, AA, B, X, beta)
    spinner_fit = calculate_CV(y = y, 
                               AA = AA_vec, 
                               W = W, 
                               X = X)
    B_hat = spinner_fit$B
    rownames(B_hat) <- 1:(dim(B_hat)[1])
    colnames(B_hat) <- 1:(dim(B_hat)[2])
    sim_results_p[[param]]$B[[(length(sim_results_n[[param]]$B) + 1)]] <- B
    sim_results_p[[param]]$B_hat[[(length(sim_results_n[[param]]$B_hat) + 1)]] <- B_hat
    sim_results_p[[param]]$real_clusters[[(length(sim_results_n[[param]]$real_clusters) + 1)]] <- real_clusters
    plot_heatmap(B_hat, "Estimated B")
  }
}

### Simulations for s
# Setting parameters
iter = 1
n = 350
p = 60
s_param = c(16, 25, 36, 49, 64)
W = ones(p, p) - diag(p)
X <- matrix(1, n, m)
beta <- runif(m)
real_clusters <- c(rep(1, dim_block), rep(2, dim_block), rep(3, dim_block), rep(4, p - 3 * dim_block))
names(real_clusters) <- 1:length(real_clusters)
sim_results_s <- list()
for(s in s_param){
  param <- paste0("s_", s)
  sim_results_s[[param]] = list(B = list(),
                          B_hat = list(),
                          real_clusters = list()
  )
}
for(s in s_param){
  B <- create_block_matrix(p, dim_block, s, noise_nonzero_prop = noise_nonzero_prop,
                           bands_nonzero_prop = bands_nonzero_prop)
  plot_heatmap(B, "True signal")
  param <- paste0("s_", s)
  for(i in 1:iter){
    AA <- generate_AA_matrix(n, p)
    AA_vec <- array(as.vector(AA))
    y <- generate_y_resp(sigma, n, AA, B, X, beta)
    spinner_fit = calculate_CV(y = y, 
                               AA = AA_vec, 
                               W = W, 
                               X = X)
    B_hat = spinner_fit$B
    rownames(B_hat) <- 1:(dim(B_hat)[1])
    colnames(B_hat) <- 1:(dim(B_hat)[2])
    sim_results_s[[param]]$B[[(length(sim_results_s[[param]]$B) + 1)]] <- B
    sim_results_s[[param]]$B_hat[[(length(sim_results_s[[param]]$B_hat) + 1)]] <- B_hat
    sim_results_s[[param]]$real_clusters[[(length(sim_results_s[[param]]$real_clusters) + 1)]] <- real_clusters
    plot_heatmap(B_hat, "Estimated B")
  }
}




# Create lists to store data ----------------------------------------------

results_n <- list()
results_s <- list()
results_p <- list()
metrics <- c("MI", "RI", "Q", "FP")
algorithms <- c("HC", "SC", "DBSCAN")
num_alg <- length(algorithms)
for(metric_name in metrics){
  for(n in n_param){
    param <- paste0("n_", n)
    iter <- length(sim_results_n[[param]]$real_clusters)
    df_results <- data.frame(matrix(0, iter, num_alg))
    colnames(df_results) <- algorithms
    results_n[[metric_name]][[param]] <- df_results
  }
  for(p in p_param){
    param <- paste0("p_", p)
    iter <- length(sim_results_p[[param]]$real_clusters)
    df_results <- data.frame(matrix(0, iter, num_alg))
    colnames(df_results) <- algorithms
    results_p[[metric_name]][[param]] <- df_results
  }
  for(s in s_param){
    param <- paste0("s_", s)
    iter <- length(sim_results_s[[param]]$real_clusters)
    df_results <- data.frame(matrix(0, iter, num_alg))
    colnames(df_results) <- algorithms
    results_s[[metric_name]][[param]] <- df_results
  }
}

clusters_n <- list()
clusters_p <- list()
clusters_s <- list()
for(alg_name in algorithms){
  for(n in n_param){
    param <- paste0("n_", n)
    real_clusters <- sim_results_n[[param]]$real_clusters[[1]]
    iter <- length(sim_results_n[[param]]$real_cluster)
    clusters_results <- data.frame(matrix(0, iter + 1, length(real_clusters)))
    clusters_results[1,] <- real_clusters
    colnames(clusters_results) <- 1:length(real_clusters)
    clusters_n[[alg_name]][[param]] <- clusters_results
  }
  for(p in p_param){
    param <- paste0("p_", p)
    real_clusters <- sim_results_p[[param]]$real_clusters[[1]]
    iter <- length(sim_results_p[[param]]$real_cluster)
    clusters_results <- data.frame(matrix(0, iter + 1, length(real_clusters)))
    clusters_results[1,] <- real_clusters
    colnames(clusters_results) <- 1:length(real_clusters)
    clusters_p[[alg_name]][[param]] <- clusters_results
  }
  for(s in s_param){
    param <- paste0("s_", s)
    real_clusters <- sim_results_s[[param]]$real_clusters[[1]]
    iter <- length(sim_results_s[[param]]$real_cluster)
    clusters_results <- data.frame(matrix(0, iter + 1, length(real_clusters)))
    clusters_results[1,] <- real_clusters
    colnames(clusters_results) <- 1:length(real_clusters)
    clusters_s[[alg_name]][[param]] <- clusters_results
  }
}


# Clustering --------------------------------------------------------------

for(n in n_param){
  param <- paste0("n_", n)
  real_clusters <- sim_results_n[[param]]$real_clusters[[1]]
  iter <- length(sim_results_n[[param]]$real_clusters)
  for(i in 1:iter){
    B_hat = sim_results_n[[param]]$B_hat[[i]]
    # hierarchical clustering
    hc_clusters <- hierarchical_clustering(B_hat)
    clusters_n$HC[[param]][i + 1, ] = hc_clusters
    results_n$MI[[param]][i, "HC"] = NMI(real_clusters, hc_clusters)
    results_n$RI[[param]][i, "HC"] = ARI(real_clusters, hc_clusters)
    # spectral clustering
    sc_clusters <- spectral_clustering(B_hat)
    clusters_n$SC[[param]][i + 1, ] = sc_clusters
    results_n$MI[[param]][i, "SC"] = NMI(real_clusters, sc_clusters)
    results_n$RI[[param]][i, "SC"] = ARI(real_clusters, sc_clusters)
  }
}

# dbscan
param = "n_450"
B_hat <- sim_results_n[[param]]$B_hat[[1]]
real_clusters <- sim_results_n[[param]]$real_clusters[[1]]
abs_B <- abs(B_hat)
dist_B <- max(abs_B) - abs_B
diag(dist_B) <- 0
minPts = 4
k = minPts - 1
kNNdistplot(B_hat, k)
eps = 5.25
abline(h=eps, col="red")
dbscan_clustering(B_hat, eps, minPts)
db_clusters <- dbscan_clustering(B_hat, eps, minPts)
clusters_n$DBSCAN[[param]][2, ] = db_clusters
results_n$MI[[param]][i, "DBSCAN"] = NMI(real_clusters, db_clusters)
results_n$RI[[param]][i, "DBSCAN"] = ARI(real_clusters, db_clusters)

for(p in p_param){
  param <- paste0("p_", p)
  real_clusters <- sim_results_p[[param]]$real_clusters[[1]]
  iter <- length(sim_results_p[[param]]$real_clusters)
  for(i in 1:iter){
    B_hat = sim_results_p[[param]]$B_hat[[i]]
    # hierarchical clustering
    hc_clusters <- hierarchical_clustering(B_hat)
    clusters_p$HC[[param]][i + 1, ] = hc_clusters
    results_p$MI[[param]][i, "HC"] = NMI(real_clusters, hc_clusters)
    results_p$RI[[param]][i, "HC"] = ARI(real_clusters, hc_clusters)
    # spectral clusteriing
    sc_clusters <- spectral_clustering(B_hat)
    clusters_p$SC[[param]][i + 1, ] = sc_clusters 
    results_p$MI[[param]][i, "SC"] = NMI(real_clusters, sc_clusters)
    results_p$RI[[param]][i, "SC"] = ARI(real_clusters, sc_clusters)
  }
}

# dbscan
param = "p_80"
B_hat <- sim_results_p[[param]]$B_hat[[1]]
real_clusters <- sim_results_p[[param]]$real_clusters[[1]]
abs_B <- abs(B_hat)
dist_B <- max(abs_B) - abs_B
diag(dist_B) <- 0
minPts = 4
k = minPts - 1
kNNdistplot(B_hat, k)
eps = 4.75
abline(h=eps, col="red")
dbscan_clustering(B_hat, eps, minPts)
db_clusters <- dbscan_clustering(B_hat, eps, minPts)
clusters_p$DBSCAN[[param]][2, ] = db_clusters
results_p$MI[[param]][i, "DBSCAN"] = NMI(real_clusters, db_clusters)
results_p$RI[[param]][i, "DBSCAN"] = ARI(real_clusters, db_clusters)

for(s in s_param){
  param <- paste0("s_", s)
  real_clusters <- sim_results_s[[param]]$real_clusters[[1]]
  iter <- length(sim_results_s[[param]]$real_clusters)
  for(i in 1:iter){
    B_hat = sim_results_s[[param]]$B_hat[[i]]
    # hierarchical clustering
    hc_clusters <- hierarchical_clustering(B_hat)
    clusters_s$HC[[param]][i + 1, ] = hc_clusters
    results_s$MI[[param]][i, "HC"] = NMI(real_clusters, hc_clusters)
    results_s$RI[[param]][i, "HC"] = ARI(real_clusters, hc_clusters)
    # spectral clustering
    sc_clusters <- spectral_clustering(B_hat)
    clusters_s$SC[[param]][i + 1, ] = sc_clusters
    results_s$MI[[param]][i, "SC"] = NMI(real_clusters, sc_clusters)
    results_s$RI[[param]][i, "SC"] = ARI(real_clusters, sc_clusters)
  }
}


# dbscan
param = "s_49"
B_hat <- sim_results_s[[param]]$B_hat[[1]]
real_clusters <- sim_results_s[[param]]$real_clusters[[1]]
abs_B <- abs(B_hat)
dist_B <- max(abs_B) - abs_B
diag(dist_B) <- 0
minPts = 4
k = minPts - 1
kNNdistplot(B_hat, k)
eps = 9
abline(h=eps, col="red")
dbscan_clustering(B_hat, eps, minPts)
db_clusters <- dbscan_clustering(B_hat, eps, minPts)
clusters_s$DBSCAN[[param]][2, ] = db_clusters
results_s$MI[[param]][i, "DBSCAN"] = NMI(real_clusters, db_clusters)
results_s$RI[[param]][i, "DBSCAN"] = ARI(real_clusters, db_clusters)



# Plots -------------------------------------------------------------------

df_MI_n <- create_df_MI(results_n$MI, "n", n_param)
plot_metric(df_MI_n, "n", "MI")  + ylim(0, 1) + theme_bw()
df_MI_p <- create_df_MI(results_p$MI, "p", p_param)
plot_metric(df_MI_p, "p", "MI") + ylim(0, 1) + theme_bw()
df_MI_s <- create_df_MI(results_s$MI, "s", s_param)
plot_metric(df_MI_s, "s", "MI") + ylim(0, 1) + theme_bw()

df_RI_n <- create_df_RI(results_n$RI, "n", n_param)
plot_metric(df_RI_n, "n", "RI") + ylim(0, 1) + theme_bw()
df_RI_p <- create_df_RI(results_p$RI, "p", p_param)
plot_metric(df_RI_p, "p", "RI") + ylim(NA, 1) + theme_bw()
df_RI_s <- create_df_RI(results_s$RI, "s", s_param)
plot_metric(df_RI_s, "s", "RI") + ylim(0, 1) + theme_bw()
