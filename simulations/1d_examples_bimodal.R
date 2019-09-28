source("../requirements.R")
source("../base_functions.R")

folder <- "../rds/1d_examples_bimodal/"
dir.create(folder, showWarnings = FALSE)

# if x is given, only generate response again
generate_bimodal <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-1.5,1.5),n,d)
  }
  f=(x[,1]-1)^2*(x[,1]+1)
  g=rep(0,n)
  g[x[,1]> -0.5]=2*sqrt(x[x[,1]> -0.5,1]+0.5)
  s=1/4+abs(x[,1])
  # response
  y=ifelse(runif(n)>0.5,f-g,f+g)+rnorm(n,0,sqrt(s))
  return(list(x=x,y=y))
}

n_fits <- 1 # total numer of I1 datasets
n_repetitions <- 500 # total numer of I2 datasets
n_train <- 5000
n_test <- 250 # to check coverage
d <- 1
k <- 100
percent_train <- 0.7
alpha <- 0.1

generate_data <- function(n,x=NULL) {generate_bimodal(n=n,d=d,x=x)}

data_test_aux <- generate_data(n_test,
                               x=matrix(seq(-1.5,1.5,length.out = n_test))) # used to fix x test
plot(data_test_aux$x[,1],data_test_aux$y)

cd_split_global <- list()
cd_split_local <- list()
dist_split <- list()
reg_split <- list()
reg_split_w <- list()

rep <- 1
bands_global <- list()
bands_local <- list()
bands_dist <- list()
bands_reg <- list()
bands_reg_w <- list()
for(n_fits_index in 1:n_fits)
{
  cat(".")
  
  data_I1 <- generate_data(n=n_train)
  
  which_train <- sample(1:length(data_I1$y),length(data_I1$y)*percent_train)
  cde_fit <- fit_density_knn(xTrain=data_I1$x[which_train,,drop=FALSE],
                             yTrain = data_I1$y[which_train,drop=FALSE],
                             xValidation=data_I1$x[-which_train,,drop=FALSE],
                             yValidation = data_I1$y[-which_train,drop=FALSE],
                             nCores=8)
  
  regression_fit <- fit_regression_knn(xTrain=data_I1$x[which_train,,drop=FALSE],
                                       yTrain = data_I1$y[which_train,drop=FALSE],
                                       xValidation=data_I1$x[-which_train,,drop=FALSE],
                                       yValidation = data_I1$y[-which_train,drop=FALSE])
  regression_fit_mean_error <- fit_regression_mean_error_knn(xTrain=data_I1$x[which_train,,drop=FALSE],
                                                             yTrain = data_I1$y[which_train,drop=FALSE],
                                                             xValidation=data_I1$x[-which_train,,drop=FALSE],
                                                             yValidation = data_I1$y[-which_train,drop=FALSE])
  
  
  for(ll in 1:n_repetitions)
  {
    print(rep)
    data_I2 <- generate_data(n=n_train)
    pred_I2 <- predict(cde_fit,data_I2$x)
    t_grid <- seq(0,max(pred_I2$CDE),length.out = 500)
    
    # CD-split global
    fit_cd_split_global <- cd_split_prediction_bands(cde_fit,
                                                     xTrain=data_I2$x,yTrain = data_I2$y,
                                                     k=nrow(data_I2$x),
                                                     xTest=data_test_aux$x,
                                                     t_grid=t_grid,
                                                     alpha=alpha)
    if(rep==1)
    {
      data_plot_original <- data.frame(x=data_I2$x,y = data_I2$y)
      pred_bands <- matrix(unlist(sapply(fit_cd_split_global$prediction_bands,function(x)fit_cd_split_global$y_grid[x])),
                           ncol=1,
                           byrow = TRUE)
      x <- matrix(unlist(sapply(1:length(fit_cd_split_global$prediction_bands),function(ii){
        rep(data_test_aux$x[ii],length.out=sum(fit_cd_split_global$prediction_bands[[ii]]))
      })),ncol=1,byrow = TRUE)
      data_plot_bands <- data.frame(x=x,
                                    y = pred_bands)
      saveRDS(list(data_plot_original=data_plot_original,
                   data_plot_bands=data_plot_bands),
              paste0(folder,"cd_global_regions.RDS"))
      
      #ggplot()+
      #  geom_point(data = data_plot_original,aes(x=x,y=y))+
      #  geom_point(data=data_plot_bands,aes(x=x,y=y),alpha=0.1,color="blue")
      
    }
    
    #CD-split local
    fit_cd_split_local <- cd_split_prediction_bands(cde_fit,
                                                    xTrain=data_I2$x,
                                                    yTrain = data_I2$y,
                                                    k=k,
                                                    xTest=data_test_aux$x,
                                                    t_grid=t_grid,
                                                    alpha=alpha)
    if(rep==1)
    {
      data_plot_original <- data.frame(x=data_I2$x,y = data_I2$y)
      pred_bands <- matrix(unlist(sapply(fit_cd_split_local$prediction_bands,function(x)fit_cd_split_local$y_grid[x])),
                           ncol=1,
                           byrow = TRUE)
      x <- matrix(unlist(sapply(1:length(fit_cd_split_local$prediction_bands),function(ii){
        rep(data_test_aux$x[ii],length.out=sum(fit_cd_split_local$prediction_bands[[ii]]))
      })),ncol=1,byrow = TRUE)
      data_plot_bands <- data.frame(x=x,
                                    y = pred_bands)
      saveRDS(list(data_plot_original=data_plot_original,
                   data_plot_bands=data_plot_bands),
              paste0(folder,"cd_local_regions.RDS"))
      
      #ggplot()+
      #  geom_point(data = data_plot_original,aes(x=x,y=y))+
      #  geom_point(data=data_plot_bands,aes(x=x,y=y),alpha=0.1,color="blue")
      
    }
    
    # Dist-split
    fit_dist_split <- dist_split_prediction_bands(cde_fit,
                                                  xTrain=data_I2$x,
                                                  yTrain = data_I2$y,
                                                  xTest=data_test_aux$x,
                                                  alpha=alpha)
    if(rep==1)
    {
      data_plot_original <- data.frame(x=data_I2$x,y = data_I2$y)
      pred_bands <- matrix(unlist(sapply(fit_dist_split$prediction_bands,function(x)fit_dist_split$y_grid[x])),
                           ncol=1,
                           byrow = TRUE)
      x <- matrix(unlist(sapply(1:length(fit_dist_split$prediction_bands),function(ii){
        rep(data_test_aux$x[ii],length.out=sum(fit_dist_split$prediction_bands[[ii]]))
      })),ncol=1,byrow = TRUE)
      data_plot_bands <- data.frame(x=x,
                                    y = pred_bands)
      saveRDS(list(data_plot_original=data_plot_original,
                   data_plot_bands=data_plot_bands),
              paste0(folder,"dist_split_regions.RDS"))
      
      #ggplot()+
      #  geom_point(data = data_plot_original,aes(x=x,y=y))+
      #  geom_point(data=data_plot_bands,aes(x=x,y=y),alpha=0.1,color="blue")
      
    }
    
    
    # Reg-split
    fit_reg_split <- reg_split_prediction_bands(regression_fit,
                                                xTrain=data_I2$x,yTrain = data_I2$y,
                                                xTest=data_test_aux$x,
                                                alpha=alpha,
                                                y_grid = pred_I2$z)
    if(rep==1)
    {
      data_plot_original <- data.frame(x=data_I2$x,y = data_I2$y)
      pred_bands <- matrix(unlist(sapply(fit_reg_split$prediction_bands,function(x)fit_reg_split$y_grid[x])),
                           ncol=1,
                           byrow = TRUE)
      x <- matrix(unlist(sapply(1:length(fit_reg_split$prediction_bands),function(ii){
        rep(data_test_aux$x[ii],length.out=sum(fit_reg_split$prediction_bands[[ii]]))
      })),ncol=1,byrow = TRUE)
      data_plot_bands <- data.frame(x=x,
                                    y = pred_bands)
      saveRDS(list(data_plot_original=data_plot_original,
                   data_plot_bands=data_plot_bands),
              paste0(folder,"reg_split_regions.RDS"))
      
      #ggplot()+
      #  geom_point(data = data_plot_original,aes(x=x,y=y))+
      #  geom_point(data=data_plot_bands,aes(x=x,y=y),alpha=0.1,color="blue")
      
    }
    
    # Reg-split weighted
    fit_reg_weighted_split <- reg_weighted_split_prediction_bands(regression_fit_mean_error$fit_mean,
                                                                  regression_fit_mean_error$fit_error,
                                                                  xTrain=data_I2$x,yTrain = data_I2$y,
                                                                  xTest=data_test_aux$x,
                                                                  alpha=alpha,
                                                                  y_grid = pred_I2$z)
    if(rep==1)
    {
      data_plot_original <- data.frame(x=data_I2$x,y = data_I2$y)
      pred_bands <- matrix(unlist(sapply(fit_reg_weighted_split$prediction_bands,function(x)fit_reg_weighted_split$y_grid[x])),
                           ncol=1,
                           byrow = TRUE)
      x <- matrix(unlist(sapply(1:length(fit_reg_weighted_split$prediction_bands),function(ii){
        rep(data_test_aux$x[ii],length.out=sum(fit_reg_weighted_split$prediction_bands[[ii]]))
      })),ncol=1,byrow = TRUE)
      data_plot_bands <- data.frame(x=x,
                                    y = pred_bands)
      saveRDS(list(data_plot_original=data_plot_original,
                   data_plot_bands=data_plot_bands),
              paste0(folder,"reg_split_w_regions.RDS"))
      
      #ggplot()+
      #  geom_point(data = data_plot_original,aes(x=x,y=y))+
      #  geom_point(data=data_plot_bands,aes(x=x,y=y),alpha=0.1,color="blue")
      
    }
    
    data_test <- generate_data(n=n_test,x=data_test_aux$x)
    
    
    # CD-split global
    bands_global[[rep]] <- cd_split_prediction_bands_evalY(fit_cd_split_global,
                                                           yTest=data_test$y)
    
    # CD-split local
    bands_local[[rep]] <- cd_split_prediction_bands_evalY(fit_cd_split_local,
                                                          yTest=data_test$y)
    
    # Dist-split
    bands_dist[[rep]] <- dist_split_prediction_bands_evalY(fit_dist_split,
                                                           yTest=data_test$y)
    
    # reg-split
    bands_reg[[rep]] <- reg_split_prediction_bands_evalY(fit_reg_split,
                                                         yTest=data_test$y)
    
    
    # reg-split weighted
    bands_reg_w[[rep]] <- reg_weighted_split_prediction_bands_evalY(fit_reg_weighted_split,
                                                                    yTest=data_test$y)
    rep <- rep+1
    gc()
  }
}
cd_split_global[[1]] <- eval_prediction_bands_conditional_coverage(xTest=data_test$x,
                                                                   bands_global,
                                                                   alpha=alpha)
saveRDS(cd_split_global,file = paste0(folder,"cd_split_global.RDS"))

cd_split_local[[1]] <- eval_prediction_bands_conditional_coverage(xTest=data_test$x,
                                                                  bands_local,
                                                                  alpha=alpha)
saveRDS(cd_split_local,file = paste0(folder,"cd_split_local.RDS"))

dist_split[[1]] <- eval_prediction_bands_conditional_coverage(xTest=data_test$x,
                                                              bands_dist,
                                                              alpha=alpha)
saveRDS(dist_split,file = paste0(folder,"dist_split.RDS"))

reg_split[[1]] <- eval_prediction_bands_conditional_coverage(xTest=data_test$x,
                                                             bands_reg,
                                                             alpha=alpha)
saveRDS(reg_split,file = paste0(folder,"reg_split.RDS"))

reg_split_w[[1]] <- eval_prediction_bands_conditional_coverage(xTest=data_test$x,
                                                               bands_reg_w,
                                                               alpha=alpha)
saveRDS(reg_split_w,file = paste0(folder,"reg_split_w.RDS"))




