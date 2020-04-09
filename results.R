# setup -----------------------------------------------------------------------

library(cowplot)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(stringr)

source("utils/setup.R")
source("utils/embedding_class.R")
source("utils/optim_embedding.R")
source("utils/results_ancillary-code.R")

# load raw data, dist matrices and tuning results and merge results

load("data/distance_matrices1.RData")
load("data/distance_matrices2.RData")
load("data/sim_data.RData")
load("data/optimized_embeddings_nongeo.RData")
load("data/optimized_embeddings_geo.RData")

# merge results
dist_list <- c(dist_list1, dist_list2) # splitted to fit in GitHub

df_res_nongeo <- do.call(rbind, opt_emb_nongeo)
df_res_geo <- do.call(rbind, opt_emb_geo)

df_res <- rbind(df_res_nongeo, df_res_geo)


l_data_l2 <- c("l2_df1_a", "l2_df1_p", "l2_df1_dep", 
               "l2_df2_a", "l2_df2_p", "l2_df2_ind")

nl_dats_l2 <- c("l2_sr_df2_a", "l2_hx_df3_a", "l2_sr_df3_a", 
                "l2_sc_df3_a", "l2_tp_df3_a")

nl_geo_dats_l2 <- c("geo_l2_sr_df2_a",  "geo_l2_hx_df3_a", "geo_l2_sr_df3_a",
                    "geo_l2_sc_df3_a",  "geo_l2_tp_df3_a")

ex_labs <- list(xlab("t"), ylab("x(t)"))
emb_labs <- list(xlab(TeX("y_1")), ylab(TeX("y_2")))


# for exact reproduction of the results
seed <- 42 


# Section 3: Study design ------------------------------------------------------

# Fig 1 -----

thm <- theme(axis.title = element_text(size = rel(1.5)),
             axis.text = element_text(size = rel(1.1)),
             axis.text.x = element_text(angle = -45)) 

smpl <- sample(1:1000, 10)

ex_data_p2 <- 
  plot_funs2(subset_dat(data_list$df2_p, smpl),
             col = data_list$df2_p$param[[1]][smpl]) +
  ex_labs + thm
ex_data_ap <- 
  plot_funs2(subset_dat(data_list$df1_dep, smpl),
             col = data_list$df1_dep$param[[1]][smpl]) +
  ex_labs + thm
ex_data_sr2 <- plot_funs2(subset_dat(data_list$sr_df2_a, smpl),
                          col = data_list$sr_df2_a$param[[1]][smpl]) +
  ex_labs + thm
ex_data_hx3 <- plot_funs2(subset_dat(data_list$hx_df3_a, smpl),
                          col = data_list$hx_df3_a$param[[1]][smpl]) +
  ex_labs + thm

# Fig 1
plot_grid(plotlist = list(ex_data_ap, ex_data_p2, ex_data_sr2, ex_data_hx3), 
          labels = c("c1-l", "p1-l", "a2-sr", "a3-hx"), 
          hjust = 0, vjust = 1, nrow = 1)

# -----

# Section 4: Results ----------------------------------------------------------

labels_l <- c("a1-l", "p1-l", "c1-l", "a2-l", "p2-l", "i2-l")
labels_n <- c("a2-sr", "a3-hx", "a3-sr", "a3-sc", "a3-tp")

# Fig 2 -----

a1_l <- embed(dist_list$l2_df1_a, "isomap", k = 8)
p1_l <- embed(dist_list$l2_df1_p, "isomap", k = 8)
c1_l <- embed(dist_list$l2_df1_dep, "isomap", k = 8)
a2_l <- embed(dist_list$l2_df2_a, "isomap", k = 8)
p2_l <- embed(dist_list$l2_df2_p, "isomap", k = 8)
i2_l <- embed(dist_list$l2_df2_ind, "isomap", k = 8)

imap_embs_l <- list(a1_l, p1_l, c1_l, a2_l, p2_l, i2_l)

fig2_plts <- vector("list", 6)
for (i in 1:6) {
  fig2_plts[[i]] <- 
    plot_emb2(imap_embs_l[[i]], size = 0.5) +
    ggtitle(TeX("ISOMAP_{8}")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    emb_labs
}

fig2_plts[[1]] <- fig2_plts[[1]] + ylim(c(-25, 25))
fig2_plts[[2]] <- fig2_plts[[2]] + ylim(c(-20, 20))
fig2_plts[[3]] <- fig2_plts[[3]] + ylim(c(-50, 50))

# Fig 2
plot_grid(plot_row(rows = 1:3, labels_l[1:3], fig2_plts[1:3]),
          plot_row(rows = 1:3, labels_l[4:6], fig2_plts[4:6]),
          nrow = 2) 

# ----

# Fig 3 -----------------------------------------------------------------------
nl_optl2_fs <- get_optvals(nl_dats_l2, meas = "auc_rnx", space = "fs") %>% arrange(method)
nl_optl2_ps <- get_optvals(nl_dats_l2, meas = "auc_rnx", space = "ps") %>% arrange(method)

fig3_r1 <- nl_optl2_fs %>% filter(method == "isomap")
fig3_r2 <- nl_optl2_ps %>% filter(method == "isomap")

fig3_r1_plts <- list_plots(fig3_r1, seed = seed, size = 0.5)
fig3_r2_plts <- list_plots(fig3_r2, seed = seed, size = 0.5)

fig3_r1_plt <- nice_plts(fig3_r1_plts)
fig3_r2_plt <- nice_plts(fig3_r2_plts)

for (i in seq_along(fig3_r1_plt)) {
  fig3_r1_plt[[i]] <- fig3_r1_plt[[i]] + emb_labs
  fig3_r2_plt[[i]] <- fig3_r2_plt[[i]] + emb_labs
}

fig3_r2_plt[[1]] <- fig3_r2_plt[[1]] + ylim(c(-500, 1000))

row1 <-
  do.call(plot_grid,
          list(
            plotlist = fig3_r1_plt,
            nrow = 1,
            ncol = 5,
            labels = c("fs", rep("", 4))
          ))

row2 <- 
  do.call(plot_grid,
          list(
            plotlist = fig3_r2_plt,
            nrow = 1,
            ncol = 5,
            labels = c("fs", rep("", 4))
          ))

# Fig 3
plot_grid(row1, row2, nrow = 2) 
# -----

# Fig 4 -----------------------------------------------------------------------
nl_optl2_ps <- get_optvals(nl_dats_l2, meas = "auc_rnx", space = "ps") %>% arrange(method)

fig4_plts <- list_plots(nl_optl2_ps %>% filter(method != "isomap"), 
                        seed = seed, 
                        size = 0.5)

fig4_plt <- nice_plts(fig4_plts)
for (i in seq_along(fig4_plt)) {
  fig4_plt[[i]] <- fig4_plt[[i]] + emb_labs
}

# Fig 4
plot_grid(plotlist = fig4_plt, nrow = 3, ncol = 5)
# -----

# Tab 4 -----------------------------------------------------------------------
tt1 <- get_optvals(nl_dats_l2, space = "fs")  %>% filter(method != "isomap")
tt2 <- get_optvals(nl_dats_l2, space = "ps")  %>% filter(method != "isomap")
tt3 <- get_optvals(nl_geo_dats_l2, space = "fs")  %>% filter(method != "isomap")
tt4 <- get_optvals(nl_geo_dats_l2, space = "ps")  %>% filter(method != "isomap")

df_geo <- data.frame(c("diffmap", "tsne", "umap"),
                     c("eps.val", "perplexity", "k"))
for (dt in paste0("tt", rep(1:4))) {
  if (dt %in% c("tt1", "tt2")) {
    excl <- "l2_sr_df3_a"
  } else {
    excl <- "geo_l2_sr_df3_a"
  }
  temp <- 
    get(dt) %>%
    filter(data != !!excl) %>%
    group_by(method) %>%
    summarize(mean = mean(k)) %>%
    pull(mean) %>%
    round(digits = 2)
  
  df_geo <- cbind(df_geo, temp)
}

names(df_geo) <- c("method", "lp", "l2_fs", "l2_ps", "geo_fs", "geo_ps")
df_diffs <- df_geo %>% mutate(l2 = abs(l2_fs - l2_ps), geo = abs(geo_fs - geo_ps))

# Tab 4
kableExtra::kable(df_diffs[, c(1, 2, 7, 8)], "latex", digits = 2, align = "l")
# -----

# Supplement ------------------------------------------------------------------

# Sup Tab 1 ------
get_optvals(nl_dats_l2, meas = "auc_rnx", space = "fs") %>% arrange(method) %>%
  summarize(mean(k))
get_optvals(nl_dats_l2, meas = "q_global", space = "fs") %>% arrange(method) %>%
  summarize(mean(k))
# ----

# Sub Fig 1 ------
sup_fig1_dat <- get_optvals(l_data_l2, space = "ps") %>% filter(method != "isomap")

sup_fig1_plots <- list_plots(sup_fig1_dat)
sup_fig1 <- nice_plts(sup_fig1_plots) 
for (i in seq_along(sup_fig1)) {
  sup_fig1[[i]] <- sup_fig1[[i]] + emb_labs
}

plot_grid(plotlist = sup_fig1, nrow = 6, ncol = 3)
# ----

# Sub Fig 2 -----
sup_fig2_dat <- get_optvals(nl_dats_l2, space = "fs") %>% filter(method != "isomap")

sup_fig2_plts <- list_plots(sup_fig2_dat)
sup_fig2_plt <- nice_plts(sup_fig2_plts)
for (i in seq_along(sup_fig2_plt)) {
  sup_fig2_plt[[i]] <- sup_fig2_plt[[i]] + emb_labs
}

plot_grid(plotlist = sup_fig2_plt, nrow = 3, ncol = 5)

# Sub Fig 3 -----
sup_fig3_dat <- get_optvals(nl_geo_dats_l2, space = "fs") %>% filter(method != "isomap")

sup_fig3_plts <- list_plots(sup_fig3_dat)
sup_fig3_plt <- nice_plts(sup_fig3_plts)
for (i in seq_along(sup_fig3_plt)) {
  sup_fig3_plt[[i]] <- sup_fig3_plt[[i]] + emb_labs
}

plot_grid(plotlist = sup_fig3_plt, nrow = 5, ncol = 3)
# ----

