fit_density_forest <- function(xTrain,yTrain,xValidation,yValidation,nCores=6)
{
  fit=fitFlexCoDE(xTrain=xTrain,zTrain=yTrain,
                  xValidation=xValidation,zValidation=yValidation,nIMax = 35,
                  regressionFunction = regressionFunction.Forest,
                  regressionFunction.extra = list(nCores=nCores))
  return(fit)
}

fit_regression_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit <- randomForest(x=rbind(xTrain,xValidation),
                      y=c(yTrain,yValidation))
  #class(fit) <- "forest"
  return(fit)
}

fit_regression_mean_error_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit_mean <- randomForest(x=xTrain,
                           y=yTrain)
  pred_val<- predict(fit_mean, xValidation)
  fit_error <- randomForest(x=xValidation,
                            y=abs(yValidation-pred_val))
  fit <- list(fit_mean=fit_mean,
              fit_error=fit_error)
  class(fit) <- "forest_weighted"
  return(fit)
}

fit_quantile_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit <- quantregForest(x=rbind(xTrain,xValidation),
                        y=c(yTrain,yValidation),nthreads = 7)
  return(fit)
}


fit_density_knn <- function(xTrain,yTrain,xValidation,yValidation,nCores=6)
{
  fit=fitFlexCoDE(xTrain=xTrain,zTrain=yTrain,
                  xValidation=xValidation,zValidation=yValidation,nIMax = 20,
                  regressionFunction = regressionFunction.NN,
                  regressionFunction.extra = list(nCores=nCores))
  return(fit)
}

fit_regression_knn <- function(xTrain,yTrain,xValidation,yValidation)
{
  k_grid <- round(seq(1,1000,length.out = 100))
  error <- rep(NA,length(k_grid))
  for(ii in 1:length(k_grid))
  {
    pred <- FNN::knn.reg(train=rbind(xTrain,xValidation),
                         y=c(yTrain,yValidation),k=k_grid[ii])$pred
    error[ii] <- mean((pred-c(yTrain,yValidation))^2)
  }
  fit <- list(x=rbind(xTrain,xValidation),
              y=c(yTrain,yValidation),
              k=k_grid[which.min(error)])
  class(fit) <- "KNN"
  return(fit)
}

predict.KNN <- function(fit,xTest)
{
  pred <- FNN::knn.reg(train=fit$x,
                       y=fit$y,k=fit$k,test=xTest)$pred
  return(pred)
}

fit_regression_mean_error_knn <- function(xTrain,yTrain,xValidation,yValidation)
{
  k_grid <- round(seq(1,1000,length.out = 100))
  error <- rep(NA,length(k_grid))
  for(ii in 1:length(k_grid))
  {
    pred <- FNN::knn.reg(train=xTrain,
                         y=yTrain,k=k_grid[ii])$pred
    error[ii] <- mean((pred-yTrain)^2)
  }
  best_k_mean <-k_grid[which.min(error)]
  
  pred_val<- FNN::knn.reg(train=xTrain,
                          y=yTrain,k=best_k_mean,test = xValidation)$pred
  
  for(ii in 1:length(k_grid))
  {
    pred <- FNN::knn.reg(train=xValidation,
                         y=abs(yValidation-pred_val),k=k_grid[ii])$pred
    error[ii] <- mean((pred-abs(yValidation-pred_val))^2)
  }
  best_k_error <-k_grid[which.min(error)]
  fit_mean <- list(xTrain=xTrain,
                   yTrain=yTrain,
                   k=best_k_mean)
  class(fit_mean) <- "KNN"
  
  fit_error <- list(xTrain=xValidation,
                    yTrain=abs(yValidation-pred_val),
                    k=best_k_error)
  class(fit_error) <- "KNN"
  
  fit <- list(fit_mean=fit_mean,fit_error=fit_error)
  class(fit) <- "KNN_weighted"
  return(fit)
}

# finds k nearest neighbors in xTrain of each xTest
which_neighbors <-  function(xTrain,xTest,k){
  return(FNN::get.knnx(data=xTrain,query = xTest,k=k)$nn.index)
}

# computes the profile of a density (OLD VERSION)
# profile_density <- function(t_grid,y_grid,cde_estimate)
# {
#   g <- rep(NA,length(t_grid))
#   for(i in seq_along(t_grid))
#   {
#     cde_estimate <- cde_estimate[cde_estimate>t_grid[i]]
#     g[i] <- sum(cde_estimate)
#   }
#   return(g*(y_grid[2]-y_grid[1]))
# }


# computes the profile of a density
profile_density <- function(t_grid,y_grid,cde_estimate)
{
  v2 <- cde_estimate[order(cde_estimate)]
  v2s <- rev(cumsum(rev(v2)))*(y_grid[2]-y_grid[1])
  v2s <- v2s[findInterval(t_grid, v2) + 1]
  v2s[which(is.na(v2s))] <- 0
  return(v2s)
}

profile_density_m <- function(t_grid,y_grid,cde_estimates)
{
  g <- matrix(NA,nrow(cde_estimates),length(t_grid))
  for(i in seq_along(t_grid))
  {
    g[,i] <- rowSums(cde_estimates*(cde_estimates>t_grid[i]))
  }
  return(g*(y_grid[2]-y_grid[1]))
}

# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
# t_grid is a grid of values for the f(.|x)
cd_split_prediction_bands <- function(cde_fit,
                                      xTrain,yTrain,
                                      k=nrow(xTrain),
                                      xTest,
                                      t_grid,
                                      alpha=0.1)
{
  pred_test <- predict(cde_fit,xTest)
  pred_train <- predict(cde_fit,xTrain)
  # observed densities:
  which_select <- cbind(1:length(yTrain),
                        which_neighbors(as.matrix(pred_train$z),
                                        as.matrix(yTrain),k=1))
  conformity_score_train <- pred_train$CDE[which_select]
  
  #prediction_bands <- list()
  prediction_bands_which_belong <- list()
  if(k!=nrow(xTrain))
  {
    ths <- rep(NA,nrow(xTest))
    g_train <- matrix(NA,nrow(xTrain),length(t_grid))
    for(ii in 1:nrow(xTrain))
    {
      g_train[ii,] <- profile_density(t_grid,pred_train$z,
                                      pred_train$CDE[ii,])
    }
    
    #aux <- profile_density_m(t_grid,pred_train$z,
    #                           pred_train$CDE)
    
    g_test <- matrix(NA,nrow(xTest),length(t_grid))
    for(ii in 1:nrow(xTest))
    {
      g_test[ii,] <- profile_density(t_grid,pred_test$z,
                                     pred_test$CDE[ii,])
    }
    neighbors <- which_neighbors(g_train,g_test,k=k)
    for(ii in 1:nrow(xTest))
    {
      ths[ii] <- quantile(conformity_score_train[neighbors[ii,]],probs=alpha)
      prediction_bands_which_belong[[ii]] <- pred_test$CDE[ii,]>=ths[ii]
    }
    
  } else {
    ths <- quantile(conformity_score_train,probs=alpha)
    for(ii in 1:nrow(xTest))
    {
      #prediction_bands[[ii]] <- pred_test$z[pred_test$CDE[ii,]>=ths]
      prediction_bands_which_belong[[ii]] <- pred_test$CDE[ii,]>=ths
    }
  }
  
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=pred_test$z,ths=ths,pred_test=pred_test))
  
  
}

# Returns fit_cd_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
cd_split_prediction_bands_evalY <- function(fit_cd_split,yTest)
{
  which_select <- cbind(1:length(yTest),
                        which_neighbors(as.matrix(fit_cd_split$pred_test$z),
                                        as.matrix(yTest),k=1))
  conformity_score_test <- fit_cd_split$pred_test$CDE[which_select]
  yTest_covered <- conformity_score_test >= fit_cd_split$ths
  fit_cd_split$yTest_covered <- yTest_covered
  fit_cd_split$pred_test <- NULL
  fit_cd_split$prediction_bands_size <- sapply(fit_cd_split$prediction_bands,function(x)mean(x))*diff(range(fit_cd_split$y_grid))
  fit_cd_split$prediction_bands <- NULL
  return(fit_cd_split)
}

lower_upper_limits <- function(density,z_grid,t)
{
  first_interval <- density>=t
  rle_x <- rle(as.numeric(first_interval))
  end = cumsum(rle_x$lengths)
  start = c(1, lag(end)[-1] + 1)
  lower <- z_grid[start[which(rle_x$values==1)[1]]]
  upper <- z_grid[end[which(rle_x$values==1)[1]]]
  return(c(lower,upper))
}

# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
# t_grid is a grid of values for the f(.|x)
cd_split_intervals_prediction_bands <- function(cde_fit,
                                                xTrain,yTrain,
                                                k=nrow(xTrain),
                                                xTest,
                                                t_grid,
                                                alpha=0.1)
{
  pred_test <- predict(cde_fit,xTest)
  pred_train <- predict(cde_fit,xTrain)
  # observed densities:
  which_train_belong_to_each_interval <- matrix(FALSE,nrow(xTrain),length(t_grid))
  for(ii in 1:nrow(xTrain))
  {
    for(tt in seq_along(t_grid))
    {
      limits <- lower_upper_limits(pred_train$CDE[ii,],pred_train$z,t_grid[tt])
      if(is.na(limits[1])|is.na(limits[2]))
      {
        break;
      }
      which_train_belong_to_each_interval[ii,tt] <- (limits[1]<=yTrain[ii])&(yTrain[ii]<=limits[2])  
      if(!which_train_belong_to_each_interval[ii,tt])
        break;
    }
  }
  prediction_bands_limits <- matrix(NA,nrow(xTest),2)
  if(k!=nrow(xTrain))
  {
    ths <- rep(NA,nrow(xTest))
    g_train <- matrix(NA,nrow(xTrain),length(t_grid))
    for(ii in 1:nrow(xTrain))
    {
      g_train[ii,] <- profile_density(t_grid,pred_train$z,
                                      pred_train$CDE[ii,])
    }
    
    g_test <- matrix(NA,nrow(xTest),length(t_grid))
    for(ii in 1:nrow(xTest))
    {
      g_test[ii,] <- profile_density(t_grid,pred_test$z,
                                     pred_test$CDE[ii,])
    }
    neighbors <- which_neighbors(g_train,g_test,k=k)
    
    
    for(ii in 1:nrow(xTest))
    {
      col <- which.min(abs(1-alpha-apply(which_train_belong_to_each_interval[neighbors[ii,],],2,mean)))
      ths[ii] <- t_grid[col]
      limits <- lower_upper_limits(density = pred_test$CDE[ii,],
                                   z_grid = pred_train$z,
                                   t=ths[ii])
      prediction_bands_limits[ii,] <- limits
      if(is.na(limits[1])|is.na(limits[2]))
      {
        prediction_bands_limits[ii,] <- c(0,0)
      }
    }
    
  } else {
    col <- which.min(abs(1-alpha-apply(which_train_belong_to_each_interval,2,mean)))
    ths <- t_grid[col]
    for(ii in 1:nrow(xTest))
    {
      limits <- lower_upper_limits(pred_test$CDE[ii,],pred_train$z,ths)
      prediction_bands_limits[ii,] <- limits
      if(is.na(limits[1])|is.na(limits[2]))
      {
        prediction_bands_limits[ii,] <- c(0,0)
      }
    }
  }
  
  return(list(prediction_bands_limits=prediction_bands_limits,
              y_grid=pred_test$z,ths=ths,pred_test=pred_test))
  
  
}


# Returns fit_cd_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
cd_split_interval_prediction_bands_evalY <- function(fit_cd_interval_split,
                                                     yTest)
{
  yTest_covered <- (fit_cd_interval_split$prediction_bands_limits[,1]<=yTest)&(yTest<=fit_cd_interval_split$prediction_bands_limits[,2])
  fit_cd_interval_split$yTest_covered <- yTest_covered
  fit_cd_interval_split$pred_test <- NULL
  fit_cd_interval_split$prediction_bands_size <- fit_cd_interval_split$prediction_bands_limits[,2]-fit_cd_interval_split$prediction_bands_limits[,1]
  fit_cd_interval_split$prediction_bands <- NULL
  return(fit_cd_interval_split)
}


cum_dist <- function(y_grid,cde_estimates,y_values)
{
  which_closest <- FNN::get.knnx(data = y_grid,
                                 query=y_values,k=1)$nn.index
  apply(as.matrix(1:nrow(cde_estimates)),1,function(xx){
    return(sum(cde_estimates[xx,1:which_closest[xx]])*diff(y_grid)[1])
  })
}

# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
dist_split_prediction_bands <- function(cde_fit,
                                        xTrain,yTrain,
                                        xTest,
                                        alpha=0.1,
                                        yTest=NULL, 
                                        median=FALSE)
{
  pred_test <- predict(cde_fit,xTest)
  pred_train <- predict(cde_fit,xTrain)
  
  cum_dist_evaluated_train <- cum_dist(pred_train$z,pred_train$CDE,yTrain)
  if(median==FALSE)
  {
    ths <-  quantile(cum_dist_evaluated_train,
                     probs = c(alpha/2,1-alpha/2))
  } else {
    ths <- quantile(abs(cum_dist_evaluated_train-0.5),probs = 1-alpha)
  }
  prediction_bands_which_belong <- list()
  FTest <- matrix(NA,nrow(xTest),length(pred_train$z))
  for (ii in 1:nrow(xTest)){
    FTest[ii,] <- cumsum(pred_test$CDE[ii,])*diff(pred_train$z)[1]
    if(median==FALSE)
    {
      prediction_bands_which_belong[[ii]] <- FTest[ii,]>=ths[1]&FTest[ii,]<=ths[2]  
    } else {
      prediction_bands_which_belong[[ii]] <- abs(FTest[ii,]-0.5)<=ths
    }
    
    
  }
  
  
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=pred_train$z,FTest=FTest,ths=ths,median=median))
  
}


# Returns fit_dist_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
dist_split_prediction_bands_evalY <- function(fit_dist_split,
                                              yTest)
{
  yTest_covered <- rep(NA,length(yTest))
  ths <- fit_dist_split$ths
  for (ii in 1:length(yTest)){
    which_closest <- which.min(abs(fit_dist_split$y_grid-yTest[ii]))
    if(fit_dist_split$median==FALSE)
    {
      yTest_covered[ii] <- fit_dist_split$FTest[ii,which_closest]>=ths[1]&fit_dist_split$FTest[ii,which_closest]<=ths[2]  
    } else {
      yTest_covered[ii] <- abs(fit_dist_split$FTest[ii,which_closest]-0.5)<=ths  
    }
  }
  fit_dist_split$yTest_covered <- yTest_covered
  fit_dist_split$prediction_bands_size <- sapply(fit_dist_split$prediction_bands,
                                                 function(x)
                                                   mean(x))*diff(range(fit_dist_split$y_grid))
  fit_dist_split$prediction_bands <- NULL
  fit_dist_split$FTest <- NULL
  return(fit_dist_split)
}


# Returns fit_dist_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
quantile_split_prediction_bands_evalY <- function(fit_quantile_split,
                                                  yTest)
{
  fit_quantile_split$yTest_covered <- fit_quantile_split$limits[,1]<=yTest & fit_quantile_split$limits[,2]>=yTest
  fit_quantile_split$prediction_bands_size <- fit_quantile_split$limits[,2]-fit_quantile_split$limits[,1]
  fit_quantile_split$limits <- NULL
  return(fit_quantile_split)
}


# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
quantile_split_prediction_bands <- function(quantile_fit,
                                            xTrain,yTrain,
                                            xTest,
                                            alpha=0.1,
                                            yTest=NULL)
{
  pred_test <- predict(quantile_fit,xTest,what=c(alpha/2,1-alpha/2))
  pred_train <- predict(quantile_fit,xTrain,what=c(alpha/2,1-alpha/2))
  E <- apply(cbind(pred_train[,1]-yTrain,yTrain-pred_train[,2]),1,max)
  Q <- quantile(E,probs = 1-alpha)
  
  limits <- cbind(pred_test[,1]-Q,pred_test[,2]+Q) 
  return(list(limits=limits))
  
}

# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
# t_grid is a grid of values for the f(.|x)
reg_split_prediction_bands <- function(reg_fit,
                                       xTrain,yTrain,
                                       xTest,
                                       alpha=0.1,
                                       y_grid)
{
  pred_test <- predict(reg_fit,xTest)
  pred_train <- predict(reg_fit,xTrain)
  # observed densities:
  conformity_score_train <- -(pred_train-yTrain)^2
  
  #prediction_bands <- list()
  prediction_bands_which_belong <- list()
  
  ths <- quantile(conformity_score_train,probs=alpha)
  for(ii in 1:nrow(xTest))
  {
    prediction_bands_which_belong[[ii]] <- -(y_grid-pred_test[ii])^2>=ths
  }
  
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=y_grid,pred_test=pred_test,ths=ths))
  
}


# Returns fit_reg_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
reg_split_prediction_bands_evalY <- function(fit_reg_split,
                                             yTest)
{
  conformity_score_test <- -(fit_reg_split$pred_test-yTest)^2
  yTest_covered <- conformity_score_test >= fit_reg_split$ths
  fit_reg_split$yTest_covered <- yTest_covered
  fit_reg_split$prediction_bands_size <- sapply(fit_reg_split$prediction_bands,
                                                function(x)
                                                  mean(x))*diff(range(fit_reg_split$y_grid))
  fit_reg_split$prediction_bands <- NULL
  fit_reg_split$pred_test <- NULL
  return(fit_reg_split)
}

# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
# t_grid is a grid of values for the f(.|x)
reg_weighted_split_prediction_bands <- function(reg_fit_mean,
                                                reg_fit_error,
                                                xTrain,yTrain,
                                                xTest,
                                                alpha=0.1,
                                                y_grid)
{
  pred_train_mean <- predict(reg_fit_mean,xTrain)
  pred_test_mean <- predict(reg_fit_mean,xTest)
  pred_train_error <- predict(reg_fit_error,xTrain)
  pred_test_error <- predict(reg_fit_error,xTest)
  
  # observed densities:
  conformity_score_train <- -abs(pred_train_mean-yTrain)/pred_train_error
  
  #prediction_bands <- list()
  prediction_bands_which_belong <- list()
  
  ths <- quantile(conformity_score_train,probs=alpha)
  for(ii in 1:nrow(xTest))
  {
    prediction_bands_which_belong[[ii]] <- -abs(y_grid-pred_test_mean[ii])/pred_test_error[ii] >=ths
  }
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=y_grid,pred_test_mean=pred_test_mean,
              pred_test_error=pred_test_error,
              ths=ths))
  
  
}

reg_weighted_split_prediction_bands_evalY <- function(fit_reg_weighted_split,
                                                      yTest)
{
  conformity_score_test <- -abs(fit_reg_weighted_split$pred_test_mean-yTest)/fit_reg_weighted_split$pred_test_error
  yTest_covered <- conformity_score_test >= fit_reg_weighted_split$ths
  fit_reg_weighted_split$yTest_covered <- yTest_covered
  fit_reg_weighted_split$prediction_bands_size <- sapply(fit_reg_weighted_split$prediction_bands,
                                                         function(x)
                                                           mean(x))*diff(range(fit_reg_weighted_split$y_grid))
  fit_reg_weighted_split$prediction_bands <- NULL
  fit_reg_weighted_split$pred_test_mean <- NULL
  fit_reg_weighted_split$pred_test_error <- NULL
  
  return(fit_reg_weighted_split)
}

# returns statistics for evaluating prediction bands
# xTest is the test sample that was used (x"s are always the same)
# bands  list of outputs of cd_split_prediction_bands
# 1-alpha is the nominal coverage
eval_prediction_bands <- function(xTest,
                                  bands,
                                  alpha)
{
  coverage <- lapply(bands, function (x) x[[c('yTest_covered')]])
  coverage <- do.call(rbind,coverage)
  coverage_mean <- colMeans(coverage)
  global_coverage <- mean(coverage_mean)
  
  B <- 200
  mean_absolute_deviation_coverage_vec <- rep(NA,B)
  for(b in 1:B)
  {
    which_sample_boot <- sample(1:nrow(coverage),
                                size = nrow(coverage),
                                replace=TRUE)
    coverage_mean_boot <- colMeans(coverage[which_sample_boot,])
    mean_absolute_deviation_coverage_vec[b] <- mean(abs(coverage_mean_boot-(1-alpha)))
  }
  mean_absolute_deviation_coverage <- mean(abs(coverage_mean-(1-alpha)))
  mean_absolute_deviation_coverage_se <- sd(mean_absolute_deviation_coverage_vec)
  
  size <- lapply(bands, function (x) x[[c('prediction_bands_size')]])
  size <- do.call(rbind,size)
  size_mean <- colMeans(size)
  average_size <- mean(size)
  average_size_vec <- rep(NA,B)
  for(b in 1:B)
  {
    which_sample_boot <- sample(1:nrow(size),
                                size = nrow(size),
                                replace=TRUE)
    average_size_vec[b] <- mean(colMeans(size[which_sample_boot,]))
  }
  average_size_se <- sd(average_size_vec)
  
  mean_absolute_deviation_size <- mean(abs(size-average_size))
  
  return(list(global_coverage=global_coverage,
              mean_absolute_deviation_coverage=mean_absolute_deviation_coverage,
              mean_absolute_deviation_coverage_se=mean_absolute_deviation_coverage_se,
              average_size=average_size,
              average_size_se=average_size_se,
              mean_absolute_deviation_size=mean_absolute_deviation_size))
}

# returns statistics for evaluating prediction bands
# xTest is the test sample that was used (x"s are always the same)
# bands  list of outputs of cd_split_prediction_bands
# 1-alpha is the nominal coverage
eval_prediction_bands_conditional_coverage <- function(xTest,
                                                       bands,
                                                       alpha)
{
  coverage <- lapply(bands, function (x) x[[c('yTest_covered')]])
  coverage <- do.call(rbind,coverage)
  coverage_mean <- colMeans(coverage)
  
  return(list(coverage_mean=coverage_mean,
              coverage_mean_se=sqrt((coverage_mean*(1-coverage_mean))/nrow(coverage))))
}

read_all_rds <- function(folder)
{
  
  data_plot <- paste0(folder,list.files(pattern = ".RDS",path = folder)) %>%
    map(readRDS)
  
  data_plot <- lapply(data_plot, function(x) {
    data <- matrix(NA,length(x),length(x[[1]]))
    colnames(data) <- names(x[[1]])
    for(ii in 1:length(x))
    {
      data[ii,] <- unlist(x[[ii]])
    }
    return(data)
  })
  names(data_plot) <- tools::file_path_sans_ext(list.files(pattern = ".RDS",path = folder))
  data_plot <- ldply(data_plot, data.frame)
  return(data_plot)
}

read_all_rds_1d <- function(folder,remove="")
{
  
  list_files <- list.files(pattern = ".RDS",path = folder)
  list_files <- list_files[-grep(remove,list_files)]
  data_plot <- paste0(folder,list_files) %>%
    map(readRDS)
  data_plot_error <- data_plot
  
  data_plot <- lapply(data_plot, function(x) {
    return(x[[1]]$coverage_mean)
  })
  names(data_plot) <- tools::file_path_sans_ext(list_files)
  data_plot <- ldply(data_plot, data.frame)
  
  data_plot_error <- lapply(data_plot_error, function(x) {
    return(x[[1]]$coverage_mean_se)
  })
  names(data_plot_error) <- tools::file_path_sans_ext(list_files)
  data_plot_error <- ldply(data_plot_error, data.frame)
  
  colnames(data_plot) <- c("method","coverage")
  data_plot$error <- data_plot_error[,2]
  return(data_plot)
}

read_all_rds_region_1d <- function(folder)
{
  list_files <- list.files(pattern = ".RDS",path = folder)
  list_files <- list_files[grep("regions",list_files)]
  lapply(paste0(folder,list_files),readRDS)
  d <- map(paste0(folder,list_files), readRDS)
  names(d) <- paste0(folder,list_files)
  names(d) <- tools::file_path_sans_ext(list_files)
  
  return(d)
}

plot_perfomance_n <- function(data_plot)
{
  
  g <- ggplot(data_plot) +
    geom_line(aes(x=n,y=global_coverage,color=.id,linetype=.id),size=2)+
    theme_minimal(base_size = 14)+ ylab("Global coverage")+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    expand_limits(y = c(0,1))+
    theme(legend.title = element_blank())
  print(g)
  
  g <- ggplot(data_plot) +
    geom_line(aes(x=n,y=mean_absolute_deviation_coverage,color=.id,linetype=.id),size=2)+
    theme_minimal(base_size = 14)+ ylab("Conditonal coverage absolute deviation")+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    expand_limits(y = 0)+
    theme(legend.title = element_blank())
  print(g)
  
  g <- ggplot(data_plot) +
    geom_line(aes(x=n,y=average_size,color=.id,linetype=.id),size=2)+
    theme_minimal(base_size = 14)+ ylab("Average size")+
    theme(legend.title = element_blank())
  print(g)
  
  g <- ggplot(data_plot) +
    geom_line(aes(x=n,y=mean_absolute_deviation_size,color=.id,linetype=.id),size=2)+
    theme_minimal(base_size = 14)+ ylab("Size absolute deviation")+
    theme(legend.title = element_blank())
  print(g)
}

plot_perfomance_n_paper <- function(data_plot,title,folder,which_plot=NULL)
{
  if(!is.null(which_plot))
  {
    data_plot <- data_plot %>% filter(.id %in% which_plot)
  }
  data_plot$.id <- revalue(data_plot$.id, c("cd_split_local"="Local CD-split", 
                                            "reg_split"="Reg-split",
                                            "dist_split"="Dist-split",
                                            "cd_split_global"="Global CD-split",
                                            "reg_split_w"="Local Reg-split",
                                            "quantile_split"="Quantile-split"))
  data_plot$mean_absolute_deviation_coverage_plus <-data_plot$mean_absolute_deviation_coverage+2*data_plot$mean_absolute_deviation_coverage_se
  data_plot$mean_absolute_deviation_coverage_minus <-data_plot$mean_absolute_deviation_coverage-2*data_plot$mean_absolute_deviation_coverage_se
  data_plot$average_size_plus <-data_plot$average_size+2*data_plot$average_size_se
  data_plot$average_size_minus <-data_plot$average_size-2*data_plot$average_size_se
  
  
  g <- ggplot(data_plot) +
    geom_line(aes(x=n,y=mean_absolute_deviation_coverage,color=.id,linetype=.id),size=2.2)+
    theme_bw()+ ylab("Conditonal coverage absolute deviation")+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    geom_ribbon(aes(x=n,ymin=mean_absolute_deviation_coverage_minus,
                    ymax=mean_absolute_deviation_coverage_plus,fill=.id),
                alpha=0.3)+
    expand_limits(y = 0)+
    theme(legend.title = element_blank())+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"),
          legend.title = element_blank(),
          legend.position = "top", 
          plot.title = element_text(size = 25, face = "bold"),
          legend.text=element_text(size=18),
          legend.key.width=unit(2,"cm"))+
    ggtitle(title)+
    guides(color=guide_legend(nrow=2))
  ggsave(paste0("figures/",folder,"/",title,"_conditional_coverage.png"),
         width = 8.6,height = 6)
  
  g <- ggplot(data_plot) +
    geom_line(aes(x=n,y=average_size,color=.id,linetype=.id),size=2.2)+
    geom_ribbon(aes(x=n,ymin=average_size_minus,
                    ymax=average_size_plus,fill=.id),
                alpha=0.3)+
    theme_bw()+ ylab("Average size")+
    theme(legend.title = element_blank())+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"),
          legend.title = element_blank(),
          legend.position = "top", 
          plot.title = element_text(size = 25, face = "bold"),
          legend.text=element_text(size=18),
          legend.key.width=unit(2,"cm"))+
    ggtitle(title)+
    guides(color=guide_legend(nrow=2))
  ggsave(paste0("figures/",folder,"/",title,"_size.png"),
         width = 8.6,height = 6)
  
}

plot_conditional_coverage_1d <- function(data_plot,folder)
{
  g <- ggplot(data = data_plot %>% filter(method%in%c("CD-split","Reg-split","Dist-split")))+
    geom_line(aes(x=x_grid,y=coverage,color=method),size=0.6)+
    geom_ribbon(aes(x=x_grid,ymin=coverage_minus,ymax=coverage_plus,fill=method),alpha=0.3)+
    theme_bw()+coord_cartesian(ylim=c(0,1))+
    ylab("Conditional coverage")+
    geom_abline(intercept=0.9,slope=0,size=1.2)+
    scale_y_continuous(labels = scales::percent)+
    xlab("x")+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"),
          legend.title = element_blank(),
          legend.position = c(0.01, 0.3), 
          legend.justification = c(0, 1),
          plot.title = element_text(size = 25, face = "bold"),
          legend.text=element_text(size=18))+
    ggtitle("Conditional coverage")
  ggsave(paste0("figures/",folder,"/1d_examples_bimodal_conditional_coverage.png"),
         plot=g,width = 10,height = 6)
}

plot_predictive_regions_1d <- function(data_plot_regions,folder)
{
  ggplot()+
    geom_point(data = data_plot_regions$cd_local_regions$data_plot_original,
               aes(x=x,y=y))+
    geom_point(data=data_plot_regions$cd_local_regions$data_plot_bands,
               aes(x=x,y=y),alpha=0.02,color="red")+
    theme_bw()+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"),
          legend.title = element_blank(),
          plot.title = element_text(size = 25, face = "bold"))+
    ggtitle("CD-split")
  ggsave(paste0("figures/",folder,"/1d_examples_bimodal_cd_split.png"),width = 10,height = 6)
  
  ggplot()+
    geom_point(data = data_plot_regions$dist_split_regions$data_plot_original,
               aes(x=x,y=y))+
    geom_point(data=data_plot_regions$dist_split_regions$data_plot_bands,
               aes(x=x,y=y),alpha=0.02,color="green")+
    theme_bw()+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"),
          legend.title = element_blank(),
          plot.title = element_text(size = 25, face = "bold"))+
    ggtitle("Dist-split")
  ggsave(paste0("figures/",folder,"/1d_examples_bimodal_dist_split.png"),width = 10,height = 6)
  
  
  ggplot()+
    geom_point(data = data_plot_regions$reg_split_regions$data_plot_original,
               aes(x=x,y=y))+
    geom_point(data=data_plot_regions$reg_split_regions$data_plot_bands,
               aes(x=x,y=y),alpha=0.02,color="blue")+
    theme_bw()+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"),
          legend.title = element_blank(),
          plot.title = element_text(size = 25, face = "bold"))+
    ggtitle("Reg-split")
  ggsave(paste0("figures/",folder,"/1d_examples_bimodal_reg_split.png"),width = 10,height = 6)
  
  
  ggplot()+
    geom_point(data = data_plot_regions$reg_split_w_regions$data_plot_original,
               aes(x=x,y=y))+
    geom_point(data=data_plot_regions$reg_split_w_regions$data_plot_bands,
               aes(x=x,y=y),alpha=0.02,color="blue")+
    theme_bw()+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20,face="bold"),
          legend.title = element_blank(),
          plot.title = element_text(size = 25, face = "bold"))+
    ggtitle("Reg-split weighted")
  ggsave(paste0("figures/",folder,"/1d_examples_bimodal_reg_split_w.png"),width = 10,height = 6)
  
  
}
