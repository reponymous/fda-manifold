optim_emb <- function(dist, params = NULL, method, meas = c("auc_rnx", "local_q"), p_space = FALSE, cores = c(2, 2), emb_out = TRUE, ...) {
  embs <- 
    if (method == "diffmap") {
      # epsval <- eps_comp(dist)
      mclapply(params, 
               function(par) embed(dist, method, eps.val = par, ...), 
               mc.cores = cores[1])
    } else if (method == "tsne") {
      mclapply(params, 
               function(par) embed(dist, method, perplexity = par, ...),
               mc.cores = cores[1])
    } else {
      mclapply(params, 
               function(par) embed(dist, method, k = par, ...),
               mc.cores = cores[1])
    }
  
  meas <- match.arg(meas)
  quality <- simplify2array(mclapply(embs, get(meas), p_space = p_space, mc.cores = cores[2]), FALSE)
  # browser()
  if (emb_out) {
    max <- list(qual = quality[which.max(quality)])
    if (!is.null(params)) max$param = params[which.max(quality)]
    list(embs = embs[[which.max(quality)]], 
         quality = quality, 
         max = max)
  } else {
    quality
  }
}

eps_comp <- function(d, ratio = 0.9, length = 75) {
  knn.m <- epsilonCompute(d$dists)
  seq(max(1, knn.m - ratio * knn.m),
      knn.m + ratio * knn.m,
      l = length)
}

local_q <- 
  function (object, ndim = 3, p_space) {
    ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
    
    if(p_space) {
      Q <- coRanking::coranking(object$p_dist,
                                ld_dist,
                                input_Xi = "dist")
    } else {
      Q <- coRanking::coranking(object$f_dist,
                                ld_dist,
                                input_Xi = "dist")
    }
    
    nQ <- nrow(Q)
    N <- nQ + 1
    
    Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) / seq_len(nQ) / N
    lcmc <- Qnx - seq_len(nQ) / nQ
    
    Kmax <- which.max(lcmc)
    
    Qlocal <- sum(lcmc[1:Kmax]) / Kmax
    return(as.vector(Qlocal))
  }

global_q <- 
  function(object, ndim = 3, p_space){
    ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
    
    if(p_space) {
      Q <- coRanking::coranking(object$p_dist,
                                ld_dist,
                                input_Xi = "dist")
    } else {
      Q <- coRanking::coranking(object$f_dist,
                                ld_dist,
                                input_Xi = "dist")
    }
    
    nQ <- nrow(Q)
    N <- nQ + 1
  
    Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) / seq_len(nQ) / N
    lcmc <- Qnx - seq_len(nQ) / nQ
  
    Kmax <- which.max(lcmc)
  
    Qglobal <- sum(lcmc[(Kmax + 1):nQ]) / (N - Kmax)
    return(as.vector(Qglobal))
}

auc_rnx <- function(object, ndim = 2, p_space = FALSE, weight = "inv") {
  rnx <- R_nx2(object, ndim = ndim, p_space = p_space)
  
  weight <- match.arg(weight, c("inv", "ln", "log", "log10"))
  switch(
    weight,
    inv   = auc_ln_k_inv(rnx),
    log   = auc_log_k(rnx),
    ln    = auc_log_k(rnx),
    log10 = auc_log10_k(rnx),
    stop("wrong parameter for weight")
  )
}

R_nx2 <- function(object, ndim = 2, p_space = FALSE) {
  # chckpkg("coRanking")
  # if (!object@has.org.data) stop("object requires original data")
  ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
  
  if(p_space) {
    Q <- coRanking::coranking(object$p_dist,
                              ld_dist,
                              input_Xi = "dist")
  } else {
    Q <- coRanking::coranking(object$f_dist,
                              ld_dist,
                              input_Xi = "dist")
  }
  
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) /
    seq_len(nQ) / N
  
  Rnx <- ((N - 1) * Qnx - seq_len(nQ)) /
    (N - 1 - seq_len(nQ))
  Rnx[-nQ]
}

auc_ln_k_inv <- function(rnx) {
  Ks <- seq_along(rnx)
  return (sum(rnx / Ks) / sum(1 / Ks))
}

auc_log_k <- function(rnx) {
  Ks <- seq_along(rnx)
  return (sum(rnx * log(Ks)) / sum(log(Ks)))
}

auc_log10_k <- function(rnx) {
  Ks <- seq_along(rnx)
  return (sum(rnx * log10(Ks)) / sum(log10(Ks)))
}


ld_dist <- function(embed) {
  ld_pts <- embed$points[, 1:3]
  ld_dist <- dist(ld_pts)
  as.matrix(ld_dist)
}

ps_dist <- function(fn_dat, k) {
  as.matrix(
    isomapdist(
      dist(
        do.call(cbind, fn_dat[[2]])
      ),
      k = k
    )
  )
}

distor_var <- function(dist) {
  var(as.vector(dist[upper.tri(dist)]))
}

distor_heat <- function(dist) {
  d1 <- d2 <- 1:ncol(dist)
  dd <- expand.grid(X = d1, Y = d2)
  dd$scores <- as.vector(dist)
  dd$scores[is.nan(dd$scores)] <- 0
  
  ggplot(dd) +
    geom_tile(aes(X, Y, fill = scores)) +
    scale_fill_viridis_c() +
    ggtitle("heatmap of distance distortion") 
  # theme(axis.text.x = element_text(angle=90))
}