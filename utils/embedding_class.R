compute_fun_dist <- function(dat, metr = c("l2", "dtw", "deriv", "sobo", "srvf"), ...) {
  if(is.list(dat)) {
    fun_mat <- dat$dat$funs
  } else {
    fun_mat <- dat
  }
  
  dists <- compute_dist(fun_mat, metr, ...)
  
  res <- list("dists" = dists, "metric" = metr, "data" = dat)
  class(res) <- "hd_dist"
  res
}

compute_dist <- function(dat, method = c("l2", "dtw", "deriv", "sobo", "srvf", "geo"), ...) {
  assert_matrix(dat, mode = "numeric") 
  
  dist <- 
    switch(method,
           "l2" = fda.usc::metric.lp(dat, ...), # slow use deriv with a = 0 instead
           "dtw" = proxy::dist(dat, method = "DTW", ...),
           "deriv" = d_dist(dat, ...),
           "srvf" = srvf_dist(dat, ...),
           "geo" = vegan::isomapdist(dat, k),
           "sobo" = sobo_dist(dat, ...))
  dist
}

embed <- function(dist_obj, method = c("isomap", "umap", "diffmap", "mds", "tsne"), k = 10, k2 = 10, ...) {
  if (class(dist_obj) == "hd_dist") {
    dist_mat <- if (is.matrix(dist_obj$dists)) {dist_obj$dists} else {as.matrix(dist_obj$dists)}
    dat <- dist_obj$data
    met <- dist_obj$metric
  } else {
    dist_mat <- dist_obj
    dat <- NULL
    met <- NULL
  }
  
  if (method == "umap" && !exists("umap_settings")) {
    umap_settings <- umap.defaults
    umap_settings$n_neighbors <- k
    umap_settings$input <- "dist"
    umap_settings$n_components <- 10 # Output space fixed to 10! 
    umap_settings$init <- "random" # only as a test: random init not recommended
  }
  
  emb <- switch(method,
                "isomap" = vegan::isomap(dist_mat, ndim = nrow(dist_mat) - 1, k = k, ...),
                "umap" = umap::umap(dist_mat, config = umap_settings, ...),
                "diffmap" =  diffuse2(dist_mat, eps.val = k, ...),
                "mds" = cmdscale(dist_mat, k = 10), # Output space fixed to 10!
                "tsne" = Rtsne(dist_mat, perplexity = k, ..., is_distance = TRUE, dims = 3)) # Output space fixed to 3!

  if (method == "tsne") class(emb) <- "tsne"
  
  points <- extract_points(emb)
  layout <- points[, 1:2]
  ld_dists <- dist(points)
  
  l_emb <-
    list(
      "emb" = emb,
      "points" = points,
      "layout" = layout, # delete and add layout-getter instead
      "f_dist" = dist_mat,
      "p_dist" = ps_dist(dist_obj$data, k = k2),
      "ld_dist" = ld_dists,
      "metric" = met,
      "method" = method,
      "k" = k,
      "data" = dat
    )
  class(l_emb) <- "embedding"
  l_emb
}

mul_emb <- function(dists, methods, k = 10, cores = 3, ...) {
  embs <- mclapply(methods, function(meth) embed(dists, meth, k = k), mc.cores = cores, ...)
  names(embs) <- methods
  embs
}


s_rho <- function(dist1, dist2) {
  # dist1,2 must be matrices
  # vectorize dist so comparison for kendall tau
  v_d1 <- c(dist1[upper.tri(dist1)])
  v_d2 <- c(dist2[upper.tri(dist1)])
  
  cor(v_d1, v_d2, method = "spearman") 
}


# shiny app vis ------------------------------------------------------------
# TODO arg check
# TODO add groups
shiny_viz <- function(l_embs, grouping = NULL) {
  methods <- vapply(l_embs, 
                    function(emb) emb[["method"]], 
                    FUN.VALUE = character(length(1)))
  
  if (length(unique(methods)) != length(methods)) {
    methods <- paste0(methods, "-", 1:length(methods))
  } 
  
  # TODO check metric slot to be non-NULL
  metrics <- vapply(l_embs,
                    function(emb) emb[["metric"]],
                    FUN.VALUE = character(length(1)))

  # TODO check datasets to be equal!!!
  
  # make sure we have unique names for every embedding/plot
  unique_nams <- paste(methods, metrics, sep = "-")
  names(l_embs) <- unique_nams
  
  funs <- get_funs(l_embs[[1]]$data) # funs in rows
  n_funs <- nrow(funs)
  grid <- get_grid(l_embs[[1]]$data)
  grid_size <- length(grid)
  
  # df_funs <- data.frame(args = rep(grid, n_funs),
  #                       vals = c(t(funs)), # transposed funs!!!
  #                       id = as.factor(rep(1:n_funs, each = grid_size)))
  # 
  # if (!is.null(grouping)) df_funs$group <- rep(grouping, grid_size)
  # 
  tf_funs <- tibble("funs" = tfd(funs, arg = grid),
                    "id" = as.factor(1:n_funs))
  if (!is.null(grouping)) tf_funs$group <- grouping
  
  shinyApp(
    ui <- fluidPage(
      fluidRow(
        lapply(unique_nams, function(meth) { 
            column(3, 
                   plotOutput(paste0("plot1_", meth), 
                              brush = paste0("plot_brush_", meth),
                              height = 500),
                   verbatimTextOutput(paste0("info_", meth)),
                   plotOutput(paste0("plot2_", meth)))
        })
      )
    ),
    server <- function(input, output) {
      lapply(unique_nams, function(meth) { # shorten by using plot_emb2!
        pts <- l_embs[[meth]]$points
        # browser()
        dat <- data.frame(dim1 = pts[, 1],
                          dim2 = pts[, 2],
                          label = as.factor(1:nrow(pts)))
        
        if (!is.null(grouping)) dat$group <- grouping
        
        output[[paste0("plot1_", meth)]] <- renderPlot({
          ggplot(dat) +
            geom_point(aes(x = dim1, y = dim2, colour = if (!is.null(grouping)) {group} else {label})) +
            theme(legend.position = "Non") +
            ggtitle(label = meth)
        })
        
        output[[paste0("info_", meth)]] <- renderPrint({
          ids <- brushedPoints(dat, input[[paste0("plot_brush_", meth)]]) 
          as.integer(ids$label)
        })

        output[[paste0("plot2_", meth)]] <- renderPlot({
          ids <- brushedPoints(dat, input[[paste0("plot_brush_", meth)]])
          ind <- as.integer(ids$label)
          mean_fun <- tf_funs %>% summarise(mean = mean(funs))
          p <- ggplot(tf_funs[ind, ]) +
               geom_spaghetti(aes(y = mean), 
                              data = tf_funs[ind, ] %>% summarise(mean = mean(funs)), 
                              size = 2) +
               geom_spaghetti(aes(y = mean), data = mean_fun)
          if (is.null(grouping)) {
            p + 
              geom_spaghetti(aes(y = funs, 
                                 colour = id)) +
              theme(legend.position = "Non")
          } else {
              p +
              geom_spaghetti(aes(y = funs, 
                                 colour = group))
          }
        })
      })
    }
  )
}

  
  
  


