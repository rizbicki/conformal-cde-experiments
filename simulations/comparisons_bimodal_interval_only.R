source("../requirements.R")
source("../base_functions.R")

folder <- "../rds/bimodal/"
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

n_fits <- 20 # total numer of I1 datasets
n_repetitions <- 100 # total numer of I2 datasets
n_each_set_grid <- c(200,500,1000,2500,5000) # size of I1 and I2
n_test <- 500 # to check coverage
d <- 20
k <- 100
percent_train <- 0.7
alpha <- 0.1

generate_data <- function(n,x=NULL) {generate_bimodal(n=n,d=d,x=x)}

data_test_aux <- generate_data(n=n_test) # used to fix x test

interval_split <- list()
for(n_each_index in 1:length(n_each_set_grid))
{
  print(n_each_index/length(n_each_set_grid))
  rep <- 1
  bands_interval <- list()
  for(n_fits_index in 1:n_fits)
  {
    cat(".")
    
    data_I1 <- generate_data(n=n_each_set_grid[n_each_index])
    
    which_train <- sample(1:length(data_I1$y),length(data_I1$y)*percent_train)
    cde_fit <- fit_density_forest(xTrain=data_I1$x[which_train,,drop=FALSE],
                                  yTrain = data_I1$y[which_train,drop=FALSE],
                                  xValidation=data_I1$x[-which_train,,drop=FALSE],
                                  yValidation = data_I1$y[-which_train,drop=FALSE])
    for(ll in 1:n_repetitions)
    {
      data_I2 <- generate_data(n=n_each_set_grid[n_each_index])
      pred_I2 <- predict(cde_fit,data_I2$x)
      t_grid <- seq(0,max(pred_I2$CDE),length.out = 500)
      
      
      
      # Interval-split 
      fit_interval_split <- cd_split_intervals_prediction_bands(cde_fit,
                                                                xTrain=data_I2$x,
                                                                yTrain = data_I2$y,
                                                                xTest=data_test_aux$x,
                                                                k=k,
                                                                t_grid=t_grid,
                                                                alpha=alpha)
      
      
      data_test <- generate_data(n=n_test,x=data_test_aux$x)
      
      
      # Interval-split 
      bands_interval[[rep]] <- cd_split_interval_prediction_bands_evalY(fit_interval_split,
                                                                        yTest=data_test$y)
      
      
      rep <- rep+1
      gc()
    }
  }
  
  interval_split[[n_each_index]] <- eval_prediction_bands(xTest=data_test$x,
                                                          bands_interval,
                                                          alpha=alpha)
  interval_split[[n_each_index]]$n <- n_each_set_grid[n_each_index]
  saveRDS(interval_split,file = paste0(folder,"interval_split.RDS"))
  rm(bands_interval)
  
}

