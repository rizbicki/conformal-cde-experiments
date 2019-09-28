source("requirements.R")
source("base_functions.R")

which_plot=c("cd_split_local",
             "reg_split",
             "dist_split",
             "cd_split_global",
             "reg_split_w",
             "quantile_split")

data_plot <- read_all_rds("rds/gamma/")
plot_perfomance_n_paper(data_plot,title = "Asymmetric",folder = "gamma/",which_plot=which_plot)

data_plot <- read_all_rds("rds/gamma_v2/")
plot_perfomance_n_paper(data_plot,title = "Asymmetric",folder = "gamma2/",which_plot=which_plot)

data_plot <- read_all_rds("rds/gaussian_het/")
plot_perfomance_n_paper(data_plot,title = "Heteroscedastic",folder = "gaussian_het/",which_plot=which_plot)

data_plot <- read_all_rds("rds/gaussian_hom/")
plot_perfomance_n_paper(data_plot,title = "Homoscedastic",folder = "gaussian_hom/",which_plot=which_plot)

data_plot <- read_all_rds("rds/bimodal/")
plot_perfomance_n_paper(data_plot,title = "Bimodal",folder = "bimodal/",which_plot=which_plot)

data_plot <- read_all_rds("rds/quantile/")
plot_perfomance_n_paper(data_plot,title = "Quantile",folder = "quantile/",which_plot=which_plot)
