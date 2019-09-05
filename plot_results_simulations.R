source("requirements.R")
source("base_functions.R")

data_plot <- read_all_rds("rds/gamma/")
plot_perfomance_n(data_plot)

data_plot <- read_all_rds("rds/gaussian_het/")
plot_perfomance_n(data_plot)

data_plot <- read_all_rds("rds/gaussian_hom/")
plot_perfomance_n(data_plot)

data_plot <- read_all_rds("rds/bimodal/")
plot_perfomance_n(data_plot)
