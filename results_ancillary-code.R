# help funs
subset_dat <-
  function(fun_data, smpl = sample(1:1000, 10)) {
    list(funs = fun_data$dat$funs[smpl, ],
         grid = fun_data$dat$grid)
  }


list_plots <- function(opt_k, seed = 42, size = 0.5) {
  plt_lst <- list()
  for (i in 1:nrow(opt_k)) {
    dat <- opt_k[[i, 1]]
    meth <- opt_k[[i, 2]]
    mea <- opt_k[[i, 5]]
    val <- opt_k[[i, 6]]
    k <- opt_k[[i, 3]]
    if (meth == "umap") {
      temp <- embed(dist_list[[dat]], method = meth, k = k, random_state = seed)
    } else {
      set.seed(seed)
      temp <- embed(dist_list[[dat]], method = meth, k = k)
    }
    plt_lst[[i]] <- 
      plot_emb2(temp, 
                title = paste(dat, meth, paste0("k: ", round(k)), sep = "-"),
                size = size) +
      ggtitle(paste(paste(meth, round(k, 2), sep = "-"), 
                    paste0(mea, " = ", round(val, 2)), 
                    sep = ": ")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  plt_lst
}

plot_row <- function(rows, labels, plot_list) {
  do.call(plot_grid, 
          list(plotlist = plot_list[rows], nrow = 1, ncol = length(rows), labels = labels)
  )
}

get_optvals <- function(dats, meas = "auc_rnx", space = NULL, method = NULL) {
  # browser()
  temp <- do.call(rbind, lapply(dats, function(dat) opt_vals(dat, meas)))
  if (!is.null(space)) temp <- temp %>% filter(space == !!space)
  if (!is.null(method)) temp <- temp %>% filter(method == !!method)
  temp
}

opt_vals <- function(dat, mea) {
  df_res %>% 
    filter(data == dat, meas == mea) %>%
    group_by(space, method) %>%
    filter(value == max(value))
}


nice_plts <- function(p_list, lbls = NULL, thm = NULL) {
  for (i in seq_along(p_list)) {
    ttl <- p_list[[i]]$labels$title
    
    temp <- str_split(ttl, " ", simplify = TRUE)
    auc_val <- temp[4]
    meth <- str_split(temp[1], "-", simplify = TRUE)
    k <- str_split(meth[2], ":", simplify = TRUE)[1]
    
    tex_title <- TeX(paste0(str_to_upper(meth[1]), 
                            "_",
                            paste0("{", k, "}"),
                            ": ",
                            "AUC_{R_{nx}} = ",
                            auc_val))
    
    p_list[[i]] <- 
      p_list[[i]] + ggtitle(tex_title)
    
    if (!is.null(lbls)) p_list[[i]] + lbls
    if (!is.null(thm)) p_list[[i]] + thm
  }
  p_list
}

plot_funs2 <- function(data, col = NULL) {
  n <- nrow(data$funs)
  # browser()
  df_dat <- data.frame(
    args = rep(data$grid, n),
    vals = c(t(data$funs)),
    id = as.factor(rep(1:n, each = length(data$grid))),
    col = rep(col, each = length(data$grid))
  )
  
  if (is.null(col)) col <- df_dat$id
  ggplot(df_dat) +
    geom_line(aes(x = args, y = vals, group = id, colour = col)) +
    theme(legend.position = "None")
}
