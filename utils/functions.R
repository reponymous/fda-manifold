# functions

# utils 

plot_emb2 <- function(embedding, grouping = NULL, labels = FALSE, title = "2d-embedding", size = 1) {
  if (inherits(embedding, "embedding")) {
    emb <- embedding$emb
  } else {
    emb <- embedding
  }
  
  if (class(emb) == "matrix") {
    pts <- emb
  } else {
    pts <- extract_points(emb, 2)
  }
  
  dat <- data.frame(dim1 = pts[, 1],
                    dim2 = pts[, 2],
                    label = as.factor(1:nrow(pts)))
  
  if (!is.null(grouping)) {
    dat$group <- grouping
  } else {   
    if (inherits(embedding, "embedding")) dat$group <- embedding$data$param[[1]]
  }
  
  p <- ggplot(dat) +
    geom_point(aes(x = dim1, 
                   y = dim2, 
                   colour = if(is.numeric(dat$group)) {group} else {label}),
               size = size) +
    theme(legend.position = "Non") +
    ggtitle(label = title)
  if (labels) p <- p + geom_text_repel(aes(x = dim1, y = dim2, label = label))
  p
}

# sobolev distance
sobo_dist <- function(data, cores = 4, ...) {
  if (is.list(data)) {
    data <- data$funs
    grid <- data$grid
  }

  n_funs <- nrow(data)
  combs <- combn(1:n_funs, 2)
  m_dist <- matrix(nrow = n_funs, ncol = n_funs)
  diag(m_dist) <- 0
    
  vals <- unlist(mclapply(1:ncol(combs), function(comb, cores) {
      sobo_metric(data[combs[1, comb], ], data[combs[2, comb], ], ...)
      },
      mc.cores = cores)
  )
    
  m_dist[upper.tri(m_dist)] <- vals
  m_dist <- t(m_dist)
  m_dist[upper.tri(m_dist)] <- vals
  m_dist
}
  
sobo_metric <- function(fun1, fun2, grid = NULL, k = 2, p = 2, a = 1) { # default L2
  if (!is.null(grid)) {
    fun1 <- tibble("funs" = tfd(fun1, arg = grid))
    fun2 <- tibble("funs" = tfd(fun2, arg = grid))
  }
    
  g = fun1$funs - fun2$funs
  p_inner <- function(f1, f2, p = 2) {
    sum(abs(tf_evaluate(f1)[[1]] * tf_evaluate(f2)[[1]])^p)^(1/p)
  }
  
  sobo_inner <- function(f1, f2, a = 1, k = k, p = 2) { # p = 2
    p_inner(f1, f2, p = p) + a * p_inner(tf_derive(f1, k), tf_derive(f2, k), p = p)
  }
    
  sqrt(sobo_inner(g, g, a = a, k = k, p = p))
    
  # derivatives <-
    #   lapply(0:k, function(k)
    #     as.matrix(tf_derive(
    #       f = x$sim_fun, order = k
    #     )))
    # 
    # dists <- lapply(derivatives, function(dat) metric.lp(dat, lp = p)^2)
    # 
    # sob_dist <- sqrt(Reduce("+", dists)) # sums over elementswise
    # sob_dist
}
  
srvf_dist <- function(data, grid, type = 2, cores = 4, ...) {
  n_funs <- nrow(data)
  combs <- combn(1:n_funs, 2)
  m_dist <- matrix(nrow = n_funs, ncol = n_funs)
  diag(m_dist) <- 0
  
  vals <- unlist(mclapply(1:ncol(combs), function(comb, ...) {
    fdasrvf::elastic.distance(data[combs[1, comb], ], data[combs[2, comb], ], time = grid, ...)[[type]]
  },
  mc.cores = cores)
  )
    
  m_dist[upper.tri(m_dist)] <- vals
  m_dist <- t(m_dist)
  m_dist[upper.tri(m_dist)] <- vals
  m_dist
}  
  
# help fun for plot embedding - S3 class to extract embedding coordinates
extract_points <- function(x, dim = 2) {
  UseMethod("extract_points")
}

# S3 method for isomap
extract_points.isomap <- function(embedding, ndim = dim(embedding$points)[2]) {
  embedding$points[, 1:ndim]
}

# S3 method for umap
extract_points.umap <- function(embedding, ndim = dim(embedding$layout)[2]) {
  embedding$layout[, 1:ndim]  
}

# S3 method for diffusionMap
extract_points.diffuse <- function(embedding, ndim = dim(embedding$X)[2]) {
  embedding$X[, 1:ndim]
}

# S3 method for mds
extract_points.matrix <- function(embedding, ndim = dim(embedding)[2]) {
  embedding[, 1:ndim]
}

# S3 method for tsne
extract_points.tsne <- function(embedding, ndim = dim(embedding$Y)[2]) {
  embedding$Y[, 1:ndim]
}


# performance evaluation ------------------------------------------------------

# computing graph laplacian for knn graph
# graph_laplac <- function(dists, n, k, con = TRUE) {
#   nn_m <-  matrix(0, n, k)
#   W <-matrix(0, ncol = n, nrow = n)
#   
#   if (class(dists) == "dist") dists <- as.matrix(dists)
#   
#   for (i in 1:n) nn_m[i, ] = FastKNN::k.nearest.neighbors(i, dists, k = k)
#   
#   for (i in 1:n) {
#     knns <- nn_m[i, ]
#     W[i, knns] <- 
#       if(con) {1} else {exp(-dists[i, knns]/(2*100))} 
#   }
#   
#   L <- diag(apply(W, 1, sum)) - W
# }
# 
# # computing embedding difference based on graph laplacian
# emb_diff <- function(dist_1, dist_2, k, con = TRUE) {
#   n <- if (class(dist_1) == "dist") {attributes(dist_1)$size} else {dim(dist_1)[1]}
#   L1 <- graph_laplac(dist_1, n, k, con = con)
#   L2 <- graph_laplac(dist_2, n, k, con = con)
#   
#   eigen1 <- eigen(L1, symmetric = TRUE)
#   eigen2 <- eigen(L2, symmetric = TRUE)
#   sqrt(sum((eigen1$values - eigen2$values)^2)) # sqrt according to Wilson and Zhu (2008), (dived by n also reasonalbe?)
# }

# fast wrapper for as.matrix
# m <- function(x) {
#   as.matrix(x)
# }
# 
# # differences based on spectra
# eigen_diff <- function(m1, m2) {
#   eigen1 <- eigen(m1, symmetric = TRUE)
#   eigen2 <- eigen(m2, symmetric = TRUE)
#   (sum(eigen1$values - eigen2$values))^2
# }
# 
# # compare several embeddings based at once (expects list of embeddings 'embs',
# # number of nearest neighbors to be regarded 'k',
# # combinations to consider "combs" (defaults to all combinations))
# compare_embs <- function(embs, k, combs = NULL) {
#   emb_nams <- names(embs)
#   if (is.null(combs)) combs <- combn(seq_along(embs), 2)
#   n_combs <- ncol(combs)
#   
#   comps <- matrix(0, nrow = n_combs, ncol = 2) 
#   for (i in 1:n_combs) {
#     arg <- c(unname(embs[combs[, i]]), k = k)
#     comps[i, ] <- c(do.call(k_tau, arg), do.call(lcmc, arg))
#   }
#   
#   colnames(comps) <- c("kendall", "lcmc")
#   if (is.null(emb_nams)) {
#     rownames(comps) <- paste0("emb-", combs[1, ], "vs", combs[2, ])
#   } else {
#     rownames(comps) <- apply(combn(emb_nams, 2), 2, paste, collapse = "_")
#   }
# 
#   comps
# }


# adapting exsting funs

diffuse2 <- function (D, eps.val = epsilonCompute(D), neigen = NULL, t = 0, 
          maxdim = 50, delta = 10^-5, maxiter = 100000) 
{
  start = proc.time()[3]
  D = as.matrix(D)
  n = dim(D)[1]
  K = exp(-D^2/(eps.val))
  v = sqrt(apply(K, 1, sum))
  A = K/(v %*% t(v))
  ind = which(A > delta, arr.ind = TRUE)
  Asp = sparseMatrix(i = ind[, 1], j = ind[, 2], x = A[ind], 
                     dims = c(n, n))
  f = function(x, A = NULL) {
    as.matrix(A %*% x)
  }
  cat("Performing eigendecomposition\n")
  if (is.null(neigen)) {
    neff = min(maxdim + 1, n)
  }
  else {
    neff = min(neigen + 1, n)
  }
  decomp = arpack(f, extra = Asp, sym = TRUE, options = list(which = "LA", maxiter = maxiter,
                                                             nev = neff, n = n, ncv = max(min(c(n, 4 * neff)))))
  psi = decomp$vectors/(decomp$vectors[, 1] %*% matrix(1, 1, 
                                                       neff))
  phi = decomp$vectors * (decomp$vectors[, 1] %*% matrix(1, 
                                                         1, neff))
  eigenvals = decomp$values
  cat("Computing Diffusion Coordinates\n")
  if (t <= 0) {
    lambda = eigenvals[-1]/(1 - eigenvals[-1])
    lambda = rep(1, n) %*% t(lambda)
    if (is.null(neigen)) {
      lam = lambda[1, ]/lambda[1, 1]
      neigen = min(which(lam < 0.05))
      neigen = min(neigen, maxdim)
      eigenvals = eigenvals[1:(neigen + 1)]
      cat("Used default value:", neigen, "dimensions\n")
    }
    X = psi[, 2:(neigen + 1)] * lambda[, 1:neigen]
  }
  else {
    lambda = eigenvals[-1]^t
    lambda = rep(1, n) %*% t(lambda)
    if (is.null(neigen)) {
      lam = lambda[1, ]/lambda[1, 1]
      neigen = min(which(lam < 0.05))
      neigen = min(neigen, maxdim)
      eigenvals = eigenvals[1:(neigen + 1)]
      cat("Used default value:", neigen, "dimensions\n")
    }
    X = psi[, 2:(neigen + 1)] * lambda[, 1:neigen]
  }
  cat("Elapsed time:", signif(proc.time()[3] - start, digits = 4), 
      "seconds\n")
  y = list(X = X, phi0 = phi[, 1], eigenvals = eigenvals[-1], 
           eigenmult = lambda[1, 1:neigen], psi = psi, phi = phi, 
           neigen = neigen, epsilon = eps.val)
  class(y) = "diffuse"
  return(y)
}


lcmc2 <- function(x, dist = NULL, k) {
  UseMethod("lcmc2")
}

lcmc2.embedding <- function(emb, k) {
  dist1 <- as.matrix(emb$ld_dist)
  dist2 <- emb$hd_dist
  
  lcmc2(dist1, dist2, k)
}

lcmc2.default <- function(dist1, dist2, k) {
  # dist1,2 must be matrices
  cork <- coranking(dist1, dist2, input = "dist")
  LCMC(cork, K = as.integer(k))
}  

lcmc_optim <- function(dist1, dist2, ks) {
  lcmcs <- unlist(mclapply(ks, function(k) {
    lcmc2(dist1, dist2, k = k)
  }))
  lcmcs[which.max(lcmcs)]
}

embeding_cor <- function(x, dist = NULL, method = "spearman", ...) { # set k w.r.t to n_funs/obs
  UseMethod("embeding_cor")
}

embeding_cor.embedding <- function(emb, method = "spearman", ...) {
  dist1 <- as.matrix(emb$ld_dist)
  dist2 <- emb$hd_dist
  
  embeding_cor(dist1, dist2, method = method, ...)
}

embeding_cor.default <- function(dist1, dist2, method = "spearman", ...) {
  # dist1,2 must be matrices
  # vectorize dist so comparison for kendall tau
  if (method == "spearman") {
    v_d1 <- c(dist1[upper.tri(dist1)])
    v_d2 <- c(dist2[upper.tri(dist1)])
    
    cor(v_d1, v_d2, method = "spearman")
  } else {
    n_obs <- nrow(dist1)
    taus <-  unlist(lapply(c(1:1000), function(x) {
        knn1 <- k.nearest.neighbors(x, dist1, ...)
        knn2 <- k.nearest.neighbors(x, dist2, ...)
        cor(knn1, knn2, method = "kendall")
      }
    ))
    sum(taus)/n_obs
  }
}

# 
# 
# sobo_dist <- function(data, grid, ...) {
#   # if (is.list(data)) {
#   #   data <- data$funs
#   #   grid <- data$grid
#   # }
#   
#   n_funs <- nrow(data)
#   combs <- combn(1:n_funs, 2)
#   m_dist <- matrix(nrow = n_funs, ncol = n_funs)
#   diag(m_dist) <- 0
#   
#   vals <- unlist(lapply(c(1:ncol(combs)), function(comb, ...) {
#     sobo_metric(data[combs[1, comb], ], data[combs[2, comb], ], ...)
#   })
#   )
#   
#   m_dist[upper.tri(m_dist)] <- vals
#   m_dist <- t(m_dist)
#   m_dist[upper.tri(m_dist)] <- vals
#   m_dist
# }

# dist of fun and derivative
# d_metric <- function(f1, f2, a = 0.5, grid) {
#   # f1 <- tibble("funs" = tfd(fun1, arg = grid))
#   # f2 <- tibble("funs" = tfd(fun2, arg = grid))
#   
#   df1 <- my_deriv(grid, f1)
#   df2 <- my_deriv(grid, f2)
# 
#   f_dist <- dist(f1, f2)
#   df_dist <- dist(df2, df2)
#   
#   f_dist + a * df_dist
# }

d_dist <- function(mat, a = 0.5, grid, d_only = TRUE) {
  
  derivs <- t(apply(mat, 1, function(f) my_deriv(grid = grid, vals = f)))
  
  dfuns <- dist(mat)
  dders <- dist(derivs)
  
  a * as.matrix(dfuns) + (1 - a) * as.matrix(dders)
}

my_deriv = function(grid, vals){
  diff(vals)/diff(grid)
}


# getter funs

get_funs <- function(data) {
  data[[1]][[1]]
}

get_grid <- function(data) {
  data[[1]][[2]]
}


# measures 

R_nx <- function(object) {
  # chckpkg("coRanking")
  # if (!object@has.org.data) stop("object requires original data")
  Q <- coRanking::coranking(object$hd_dist,
                            as.matrix(object$ld_dist),
                            input = "dist")
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) /
    seq_len(nQ) / N
  
  Rnx <- ((N - 1) * Qnx - seq_len(nQ)) /
    (N - 1 - seq_len(nQ))
  Rnx[-nQ]
}

R_nx2 <- function(object, ndim = 2) {
  # chckpkg("coRanking")
  # if (!object@has.org.data) stop("object requires original data")
  ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
  
  Q <- coRanking::coranking(object$hd_dist,
                            ld_dist,
                            input = "dist")
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) /
    seq_len(nQ) / N
  
  Rnx <- ((N - 1) * Qnx - seq_len(nQ)) /
    (N - 1 - seq_len(nQ))
  Rnx[-nQ]
}


auc_rnx <- function(object, ndim = 2, weight = "inv") {
  rnx <- R_nx2(object, ndim = ndim)
  
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

local_q <- function (object, ndim = 2) {
  # if (!object@has.org.data) stop("object requires original data")
  # chckpkg("coRanking")
  
  # Q <- coRanking::coranking(object@org.data,
  #                           object@data@data[, seq_len(ndim), drop = FALSE])
  ld_dist <- as.matrix(dist(object$points[, 1:ndim]))
  
  Q <- coRanking::coranking(object$hd_dist,
                            ld_dist,
                            input = "dist")
  nQ <- nrow(Q)
  N <- nQ + 1
  
  Qnx <- diag(apply(apply(Q, 2, cumsum), 1, cumsum)) / seq_len(nQ) / N
  lcmc <- Qnx - seq_len(nQ) / nQ
  
  Kmax <- which.max(lcmc)
  
  Qlocal <- sum(lcmc[1:Kmax]) / Kmax
  return(as.vector(Qlocal))
}

plotly_viz <- function(emb, ..., size = 0.1) {
  plotly::plot_ly(x = emb$points[, 1], y = emb$points[, 2], z = emb$points[, 3],
                  size = size,
                  type = "scatter3d", ...)
}
