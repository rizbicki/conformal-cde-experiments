source("requirements.R")
source("base_functions.R")

folder <- "1d_examples_bimodal/"
data_plot <- read_all_rds_1d(paste0("rds/",folder),
                             remove="regions")
data_plot_regions <- read_all_rds_region_1d(paste0("rds/",folder))

n_test <- 250
x_grid <- seq(-1.5,1.5,length.out = n_test)
data_plot$x_grid <- x_grid
data_plot$coverage_plus <-data_plot$coverage+2*data_plot$error
data_plot$coverage_minus <-data_plot$coverage-2*data_plot$error
data_plot$method <- revalue(data_plot$method, c("cd_split_local"="CD-split", 
                                                "reg_split"="Reg-split",
                                                "dist_split"="Dist-split"))

plot_conditional_coverage_1d(data_plot,folder)


plot_predictive_regions_1d(data_plot_regions,folder)

