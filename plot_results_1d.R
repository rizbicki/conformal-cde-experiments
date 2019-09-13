source("requirements.R")
source("base_functions.R")

data_plot <- read_all_rds_1d("rds/1d_examples_bimodal/")
n_test <- 500
x_grid <- seq(-1.5,1.5,length.out = n_test)
data_plot$x_grid <- x_grid
data_plot$coverage_plus <-data_plot$coverage+data_plot$error
data_plot$coverage_minus <-data_plot$coverage-data_plot$error
#plot(x_grid,data_plot %>%
#       filter(method=="cd_split_local") %>%
#       select(coverage) %>% pull())

ggplot()+
  geom_line(data = data_plot,aes(x=x_grid,y=coverage,color=method))+
theme_bw()+coord_cartesian(ylim=c(0,1))


ggplot(data = data_plot %>% filter(method=="reg_split"))+
  geom_line(aes(x=x_grid,y=coverage,color=method))+
  geom_ribbon(aes(x=x_grid,ymin=coverage_minus,ymax=coverage_plus,fill=method))+
  theme_bw()+coord_cartesian(ylim=c(0,1))
