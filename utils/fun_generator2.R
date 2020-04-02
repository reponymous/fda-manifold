library(tibble)
library(tidyfun)
library(dplyr)
library(ggplot2)
library(splines)

fun_generator2 <- function(n_funs, 
                           grid = seq(0, 1, l = 100), 
                           process, 
                           proc_param = NULL,  
                           warp_fun = function(t) {t}, # id warping
                           warp_param = NULL,
                           tfd = FALSE,
                           ...) {
  
  vals <- matrix(nrow = n_funs, ncol = length(grid))
  
  for (i in 1:n_funs) {
    wp <- lapply(warp_param, `[[`, i)
    pp <- lapply(proc_param, `[[`, i)
    
    warp_grid <- do.call(warp_fun, c(list(grid), wp))
    
    vals[i, ] <- do.call(process, c(list(warp_grid), pp))
  }
  
  if (tfd) {
    tibble("funs" = tfd(vals, arg = grid))
  } else {
    list(funs = vals, grid = grid) 
    # fruther possible outputs, warping fun, warpings, dat name, base fun
  }
}


plot_funs <- function(data) {
  n <- nrow(data$funs)
  df_dat <- data.frame(
    args = rep(data$grid, n),
    vals = c(t(data$funs)),
    id = as.factor(rep(1:n, each = length(data$grid)))
  )
  ggplot(df_dat) +
    geom_line(aes(x = args, y = vals, group = id, colour = id)) +
    theme(legend.position = "None")
}
