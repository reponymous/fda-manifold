optim_emb2 <- function(dist, method, params, cores = 2, ...) {
  quals <- 
    mclapply(params, 
             function(par) compute_quals(dist, method, par = par),
             mc.cores = cores)
    # umap has special settings: output space fixed and random init
  
  df_quals <- do.call(rbind, quals)
  
  # out <- vector("list", 2)
  # names(out) <- c("ps", "fs")
  # 
  # for (sp in names(out)) {
  #   qual <- quals[[paste0("quality_", sp)]]
  #   
  #   k_inds <- apply(qual, 1, which.max)
  # 
  #   max <- list(meas = apply(qual, 1, max),
  #               params = vapply(k_inds, 
  #                               function(ind) params[ind], 
  #                               FUN.VALUE = numeric(1)))
  # 
  #   emb <- lapply(k_inds, function(ind) (embs[[ind]]))
  #     
  #   out[[sp]] <- list(embs = emb,
  #                     quality = qual,
  #                     max = max)        
  # }
  # out
}

compute_quals <- function(dist, method, par, ...) {
  # TODO check seeding and ggf store it for reproducability
  embs <- embed(dist[[1]], method, k = par, ...)
  
  qual_ps <- quality(embs, p_space = TRUE)
  qual_fs <- quality(embs, p_space = FALSE)
  
  data.table(data = rep(names(dist), 6),
             method = rep(method, 6),
             k = rep(par, 6),
             space = rep(c("ps", "fs"), each = 3), 
             meas = names(c(qual_ps, qual_fs)), 
             value = unname(c(qual_ps, qual_fs)))
}

quality <- function(embs, p_space) {
  auc <- auc_rnx(embs, p_space = p_space, ndim = 3, weight = "log10")
  q_local <- local_q(embs, p_space = p_space, ndim = 3)  
  q_global <- global_q(embs, p_space = p_space, ndim = 3)
  
  c(auc_rnx = auc, q_local = q_local, q_global = q_global)
}

get_opt_embmeas <- function(opt_emb, method) {
  fs <- opt_emb[[method]]$fs$max
  ps <- opt_emb[[method]]$ps$max
  
  df_res <- data.frame(space = rep(c("fs", "ps"), each = 3),
                       meas = c(names(fs$meas), names(ps$meas)),
                       val = c(fs$meas, ps$meas),
                       params = c(fs$params, ps$params))
  df_res
}

