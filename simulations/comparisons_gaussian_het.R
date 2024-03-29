source("../requirements.R")
source("../base_functions.R")

folder <- "../rds/gaussian_het/"
dir.create(folder, showWarnings = FALSE)

# if x is given, only generate response again
generate_het_gaussian <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-5,5),n,d)
  }
  # response
  y=rnorm(nrow(x),x[,1],abs(x[,1]))
  return(list(x=x,y=y))
}

n_fits <- 15 # total numer of I1 datasets
n_repetitions <- 250 # total numer of I2 datasets
n_each_set_grid <- c(200,500,1000,2500,5000) # size of I1 and I2
n_test <- 500 # to check coverage
d <- 20
k <- 100
percent_train <- 0.7
alpha <- 0.1

generate_data <- function(n,x=NULL) {generate_het_gaussian(n=n,d=d,x=x)}

data_test_aux <- generate_data(n=n_test) # used to fix x test
cd_split_global <- list()
cd_split_local <- list()
dist_split <- list()
quantile_split <- list()
reg_split <- list()
reg_split_w <- list()
for(n_each_index in 1:length(n_each_set_grid))
{
  print(n_each_index/length(n_each_set_grid))
  rep <- 1
  bands_global <- list()
  bands_local <- list()
  bands_dist <- list()
  bands_quantile <- list()
  bands_reg <- list()
  bands_reg_w <- list()
  for(n_fits_index in 1:n_fits)
  {
    cat(".")
    
    data_I1 <- generate_data(n=n_each_set_grid[n_each_index])
    
    which_train <- sample(1:length(data_I1$y),length(data_I1$y)*percent_train)
    cde_fit <- fit_density_forest(xTrain=data_I1$x[which_train,,drop=FALSE],
                                  yTrain = data_I1$y[which_train,drop=FALSE],
                                  xValidation=data_I1$x[-which_train,,drop=FALSE],
                                  yValidation = data_I1$y[-which_train,drop=FALSE])
    
    quantile_fit <- fit_quantile_forest(xTrain=data_I1$x[which_train,,drop=FALSE],
                                        yTrain = data_I1$y[which_train,drop=FALSE],
                                        xValidation=data_I1$x[-which_train,,drop=FALSE],
                                        yValidation = data_I1$y[-which_train,drop=FALSE])
    
    regression_fit <- fit_regression_forest(xTrain=data_I1$x[which_train,,drop=FALSE],
                                            yTrain = data_I1$y[which_train,drop=FALSE],
                                            xValidation=data_I1$x[-which_train,,drop=FALSE],
                                            yValidation = data_I1$y[-which_train,drop=FALSE])
    regression_fit_mean_error <- fit_regression_mean_error_forest(xTrain=data_I1$x[which_train,,drop=FALSE],
                                                                  yTrain = data_I1$y[which_train,drop=FALSE],
                                                                  xValidation=data_I1$x[-which_train,,drop=FALSE],
                                                                  yValidation = data_I1$y[-which_train,drop=FALSE])
    
    
    
    for(ll in 1:n_repetitions)
    {
      data_I2 <- generate_data(n=n_each_set_grid[n_each_index])
      pred_I2 <- predict(cde_fit,data_I2$x)
      t_grid <- seq(0,max(pred_I2$CDE),length.out = 500)
      
      # CD-split global
      fit_cd_split_global <- cd_split_prediction_bands(cde_fit,
                                                       xTrain=data_I2$x,yTrain = data_I2$y,
                                                       k=nrow(data_I2$x),
                                                       xTest=data_test_aux$x,
                                                       t_grid=t_grid,
                                                       alpha=alpha)
      # CD-split local
      fit_cd_split_local <- cd_split_prediction_bands(cde_fit,
                                                      xTrain=data_I2$x,
                                                      yTrain = data_I2$y,
                                                      k=k,
                                                      xTest=data_test_aux$x,
                                                      t_grid=t_grid,
                                                      alpha=alpha)
      
      # Dist-split 
      fit_dist_split <- dist_split_prediction_bands(cde_fit,
                                                    xTrain=data_I2$x,yTrain = data_I2$y,
                                                    xTest=data_test_aux$x,
                                                    alpha=alpha)
      
      # Quantile-split 
      fit_quantile_split <- quantile_split_prediction_bands(quantile_fit,
                                                            xTrain=data_I2$x,
                                                            yTrain = data_I2$y,
                                                            xTest=data_test_aux$x,
                                                            alpha=alpha)
      
      
      data_test <- generate_data(n=n_test,x=data_test_aux$x)
      
      # Reg-split
      fit_reg_split <- reg_split_prediction_bands(regression_fit,
                                                  xTrain=data_I2$x,yTrain = data_I2$y,
                                                  xTest=data_test_aux$x,
                                                  alpha=alpha,
                                                  y_grid = pred_I2$z)
      
      # Reg-split
      fit_reg_weighted_split <- reg_weighted_split_prediction_bands(regression_fit_mean_error$fit_mean,
                                                                    regression_fit_mean_error$fit_error,
                                                                    xTrain=data_I2$x,yTrain = data_I2$y,
                                                                    xTest=data_test_aux$x,
                                                                    alpha=alpha,
                                                                    y_grid = pred_I2$z)
      
      #print(rep)
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
      
      # Quantile-split 
      bands_quantile[[rep]] <- quantile_split_prediction_bands_evalY(fit_quantile_split,
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
  cd_split_global[[n_each_index]] <- eval_prediction_bands(xTest=data_test$x,
                                                           bands_global,
                                                           alpha=alpha)
  cd_split_global[[n_each_index]]$n <- n_each_set_grid[n_each_index]
  saveRDS(cd_split_global,file = paste0(folder,"cd_split_global.RDS"))
  rm(bands_global)
  gc()
  
  cd_split_local[[n_each_index]] <- eval_prediction_bands(xTest=data_test$x,
                                                          bands_local,
                                                          alpha=alpha)
  cd_split_local[[n_each_index]]$n <- n_each_set_grid[n_each_index]
  saveRDS(cd_split_local,file = paste0(folder,"cd_split_local.RDS"))
  rm(bands_local)
  gc()
  
  dist_split[[n_each_index]] <- eval_prediction_bands(xTest=data_test$x,
                                                      bands_dist,
                                                      alpha=alpha)
  dist_split[[n_each_index]]$n <- n_each_set_grid[n_each_index]
  saveRDS(dist_split,file = paste0(folder,"dist_split.RDS"))
  rm(bands_dist)
  gc()
  
  quantile_split[[n_each_index]] <- eval_prediction_bands(xTest=data_test$x,
                                                          bands_quantile,
                                                          alpha=alpha)
  quantile_split[[n_each_index]]$n <- n_each_set_grid[n_each_index]
  saveRDS(quantile_split,file = paste0(folder,"quantile_split.RDS"))
  
  reg_split[[n_each_index]] <- eval_prediction_bands(xTest=data_test$x,
                                                     bands_reg,
                                                     alpha=alpha)
  reg_split[[n_each_index]]$n <- n_each_set_grid[n_each_index]
  saveRDS(reg_split,file = paste0(folder,"reg_split.RDS"))
  rm(bands_reg)
  gc()
  
  reg_split_w[[n_each_index]] <- eval_prediction_bands(xTest=data_test$x,
                                                       bands_reg_w,
                                                       alpha=alpha)
  reg_split_w[[n_each_index]]$n <- n_each_set_grid[n_each_index]
  saveRDS(reg_split_w,file = paste0(folder,"reg_split_w.RDS"))
  rm(bands_reg_w)
  gc()
  
}
