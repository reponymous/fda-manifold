# generating data
source("utils/setup.R")

n_funs <- 1000
grid_size <- 200
grid <- seq(0, 1, l = grid_size)

# base_fun --------------------------------------------------------------------

gauss1 <- function(s) {
  dnorm(s, mean = 0.25, sd = 0.1) 
}

gauss2 <- function(s) {
  dnorm(s, mean = 0.75, sd = 0.1)
}

bas_fun1 <- function(s, a) {a * (gauss1(s) + gauss2(s))}
bas_fun2 <- function(s, a1, a2) {a1 * gauss1(s) + a2 * gauss2(s)}
bas_fun3 <- function(s, a1, a2, a3) {(a1 * gauss1(s) + a2 * gauss2(s)) + a3}

set.seed(28)

# parameters for linear space -------------------------------------------------
warp1 <- runif(n_funs, min = 0.5, max = 3)
warp2 <- runif(n_funs, min = 0.5, max = 3)
amp1 <- runif(n_funs, min = 0.5, max = 3)
amp2 <- runif(n_funs, min = 0.5, max = 3)

#  nonlinear parameter spaces/manifolds ---------------------------------------
library(dimRed)

sr <- loadDataSet("Swiss Roll", n = n_funs)
sr2 <- loadDataSet("Swiss Roll", n = n_funs)
hx <- loadDataSet("Helix", n = n_funs)
sc <- loadDataSet("3D S Curve", n = n_funs)
tp <- loadDataSet("Twin Peaks", n = n_funs)


# 1. linear parameter space settings

# linear warping function 
lin_warp <- function(t, s, a1 = NULL) {
  # stopifnot(s %in% c(0:1))
  grid_length <- length(t)
  if (is.null(a1)) a1 <- runif(1, min = 0, max = 1/s)
  
  p1 <- t[1:(floor(length(t) * s))] * a1
  p2 <- (t[(floor(length(t) * s) + 1):grid_length] - 1) * ((s * a1 - 1)/(s - 1)) +  1
  c(p1, p2)
} 

data_list <- vector("list", 12)
names(data_list) <- c(
  "df1_a",
  "df1_p",
  "df1_dep",
  "df2_a",
  "df2_p",
  "df2_ind",
  "sr_df2_a",
  "sr_df2_dep",
  "hx_df3_a",
  "sr_df3_a",
  "sc_df3_a",
  "tp_df3_a"
)

# L1: amplitude variation in 1 df ---------------------------------------------
params_df1 <- list(amp1)

data_list[["df1_a"]] <- 
  list(dat = fun_generator2(n_funs, grid, 
                            bas_fun1,
                            # proc_param (for data generating process parameters) 
                            # is the slot amplitude variation parameters 
                            # (bad naming to be changed)
                            proc_param = params_df1), 
       param = params_df1)

# L2: phase variation in 1 df with linear warping -----------------------------
# linear warping with break point 0.5
lin_warp_param <- list(rep(0.5, n_funs),
                       runif(n_funs, min = 0.01, max = 0.99))

data_list[["df1_p"]] <- 
  list(
    dat = fun_generator2(
      n_funs,
      grid,
      bas_fun1,
      proc_param = list(rep(1, n_funs)),
      warp_fun = lin_warp,
      warp_param = lin_warp_param
    ),
    param = lin_warp_param
  )

# L3: dependent amp and phase variation in 1 df with warping via power fun ----
data_list[["df1_dep"]] <-
  list(
    dat = fun_generator2(
      n_funs,
      grid,
      bas_fun1,
      proc_param = params_df1,
      warp_fun = function(t, w) {
        t ^ w
      },
      warp_param = params_df1
    ),
    param = params_df1
  )

# L4: amplitude variation in 2 df ---------------------------------------------
data_list[["df2_a"]] <- 
  list(
    dat = fun_generator2(n_funs, grid,
                         bas_fun2,
                         proc_param = list(amp1, amp2)),
    param = list(amp1, amp2)
  )

# L5: phase variation in 2 df with warping via beta cdf -----------------------

data_list[["df2_p"]] <-
  list(
    dat = fun_generator2(
      n_funs,
      grid,
      bas_fun1,
      proc_param = list(rep(1, n_funs)), 
      warp_fun = pbeta,
      warp_param = list(warp1, warp2)
    ),
    param = list(warp1, warp2)
  )


# L6: independent amp and phase variation (2 df) with linear warping ----------

data_list[["df2_ind"]] <-
  list(
    dat = fun_generator2(
      n_funs,
      grid,
      bas_fun1,
      proc_param = params_df1,
      warp_fun = lin_warp,
      warp_param = lin_warp_param
    ),
    param = list(amp = params_df1[[1]],
                 warp = lin_warp_param[[2]])
  )

# 2. nonlinear parameter space settings

# N1: amplitude variation in 2 df, swiss role 1d ------------------------------
sr_param <- list(sr@data[, 1], # 1d swiss role embedded in 2d space
                 sr@data[, 3])

data_list[["sr_df2_a"]] <-
  list(
    dat = fun_generator2(
    n_funs,
    grid = grid,
    bas_fun2,
    proc_param = sr_param
  ),
  param = sr_param
)

# N3: amplitude variation in 3df, helix 1d ------------------------------------
hx_a_params <- list(hx@data[, 1], # 1d helix embedded in 3d space
                    hx@data[, 2],
                    hx@data[, 3])

data_list[["hx_df3_a"]] <- 
  list(
    dat = fun_generator2(n_funs,
                         grid = grid,
                         bas_fun3,
                         proc_param = hx_a_params),
    param = hx_a_params
  )


# N4: amplitude variation in 3 df, swiss role 2d ------------------------------
sr2_params <- list(sr2@data[, 1], # 2d swiss role embedded in 3d space
                   sr2@data[, 2],
                   sr2@data[, 3])

data_list[["sr_df3_a"]] <- 
  list(
    dat = fun_generator2(n_funs,
                         grid = grid,
                         bas_fun3,
                         proc_param = sr2_params),
    param = sr2_params
  )

# N5: amplitude variation in 3df, S-curve 2d ----------------------------------
sc_params <- list(sc@data[, 1], # 2d s-curve embedded in 3d space
                    sc@data[, 2],
                    sc@data[, 3])

data_list[["sc_df3_a"]] <-
  list(
    dat = fun_generator2(n_funs,
                         grid = grid,
                         bas_fun3,
                         proc_param = sc_params),
    param = sc_params
  )

# N6: amplitude variation in 3df, twin peaks 2d -------------------------------
tp_params <- list(tp@data[, 1], # 2d two-peak surface embedded in 3d space
                  tp@data[, 2],
                  tp@data[, 3])

data_list[["tp_df3_a"]] <-
  list(
    dat = fun_generator2(n_funs,
                         grid = grid,
                         bas_fun3,
                         proc_param = tp_params),
    param = tp_params
  )

# save(data_list, file = "paper/sim_data.RData")
