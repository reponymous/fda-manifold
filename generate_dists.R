# compute dists
source("utils/setup.R")

load("data/sim_data.RData")


# computing fun space dist l2 and dtw -----------------------------------------
d_l2 <-
  mclapply(data_list, function(dat)
    compute_fun_dist(dat, "l2"), mc.cores = 12)

# save(d_l2, file = "distance_matrices.RData")



# computing geodesic distance and minimal neighborhood sizes ------------------

# help function ---------------------------------------------------------------

compute_geo_dist <- function(dists) {
  min_ks_i <- min_ks_u <- vector("integer", length(dists))
  d_geo <- vector("list", length(dists))
  names(d_geo) <- names(min_ks_i) <- names(min_ks_u) <- names(dists)
  
  for (dis in names(dists)) {
    for (min_k_i in 3:25) {
      val <- 
        tryCatch(embed(dist_obj = dists[[dis]], 
                       method = "isomap", 
                       k = min_k_i, 
                       k2 = 20), 
                 error = function(e) {"k-too-small"})
      if (!is.character(val)) break
    }
    min_ks_i[[dis]] <- min_k_i
    
    for (min_k_u in 3:25) {
      val <- 
        tryCatch(embed(dist_obj = dists[[dis]], 
                       method = "umap", 
                       k = min_k_u, 
                       k2 = 20), 
                 error = function(e) {"k-too-small"})
      if (!is.character(val)) break
    }
    min_ks_u[[dis]] <- min_k_u
    
    d_geo[[dis]] <- list(
      dists = isomapdist(dists[[dis]]$dists,
                         k = min_k_i), # use minimal k of isomap, since geo dist is computed via isomap routine
      metric = "geo",
      data = dists[[dis]]$data
    )
    class(d_geo[[dis]]) <- "hd_dist"
  }
  list(dist_geo = d_geo, min_ks_u = min_ks_u, min_ks_i = min_ks_i)
}

# constructing the neighborhood graph can fail if parameter k is choosen to small
# To avoid failures in following optimization routine be determine the minimal 
# working k for each data set based on the two distance measrues.
# In addition the geodesic distance in function space is computed.

out_l2 <- compute_geo_dist(d_l2)

d_geo_l2 <- out_l2$dist_geo

dist_list <- c(dist_list, d_geo_l2)
names(dist_list) <-
  paste(rep(c("l2", "geo_l2"), each = 12), 
        names(data_list), 
        sep = "_")

min_ks <- c(out_l2[2:3], out_dtw[2:3])
min_ks_i <- c(min_ks[[2]], min_ks[[4]])
# save(dist_list, file = "distance_matrices.RData")
# save(min_ks_i, file = "min_ks_i.RData")
