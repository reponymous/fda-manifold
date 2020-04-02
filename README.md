# fda-manifold

For reproduction: Download the repo as a RStudio Project

To reproduce the Tables and Figures in the paper use file: results.R

To reproduce the Experiment from scratch:

- data_generating_processes.R: computes the the functional data sets (setting a1-l - a3-tp)
- generate_dists.R: computes distances matrices based on the data (This is takes some time!)
- optimization_geo.R and optimization_nongeo.R tune the embeddings (This takes even more time: We split the computations to compute on two independent machine with 30 kernel each, i.e. in total 60 kernel. Computations times roughly 6-8 hours if computed with 30 kernels for one file. 
