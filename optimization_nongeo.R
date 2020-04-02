source("utils/setup.R")
source("utils/optim_embedding.R")
source("utils/optim_emb2.R")

load("data/distance_matrices.RData")
load("data/min_ks_i.RData")

methods <- c("isomap", "umap", "diffmap", "tsne")

# computations are very time demanding about 6 on system with 30 kernels
# caution !!! set core according to system

l_dists <- dist_list[1:11] 
opt_embeddings <- vector("list", length = length(l_dists))
names(opt_embeddings) <- names(l_dists)

for (dat_dist in names(l_dists)) {
  print(paste0("here we are: ", dat_dist))
  
  params <- list(
    "isomap" = seq(min_ks_i[[dat_dist]], 975L, by = 3), # by = 3
    "umap" = seq(3L, 975L, by = 3), # by = 3
    "diffmap" = eps_comp(l_dists[[dat_dist]], ratio = 0.9, length = 300), # l = 300
    "tsne" = 3L:floor(999/3) # 999/3
  )
  
  opt_embs <- vector("list", length = length(methods))
  names(opt_embs) <- methods
  for (meth in methods) {
    pars <- params[meth]
    opt_embs[[meth]] <- optim_emb2(l_dists[dat_dist], 
                                   method = meth, 
                                   params = params[[meth]], 
                                   cores = 30L)  #  <----------------------- set cores 
  }
  opt_embeddings[[dat_dist]] <- do.call(rbind, opt_embs)
}

save(opt_embeddings, file = "opt_embs_nongeo.RData")