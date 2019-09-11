fit_density_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit=fitFlexCoDE(xTrain=xTrain,zTrain=yTrain,
                  xValidation=xValidation,zValidation=yValidation,nIMax = 20,
                  regressionFunction = regressionFunction.Forest,
                  regressionFunction.extra = list(nCores=6))
  return(fit)
}

fit_regression_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit <- randomForest(x=rbind(xTrain,xValidation),
                      y=c(yTrain,yValidation))
  return(fit)
}

fit_regression_mean_error_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit_mean <- randomForest(x=xTrain,
                           y=yTrain)
  pred_val<- predict(fit_mean, xValidation)
  fit_error <- randomForest(x=xValidation,
                            y=abs(yValidation-pred_val))
  return(list(fit_mean=fit_mean,
              fit_error=fit_error))
}

# finds k nearest neighbors in xTrain of each xTest
which_neighbors <-  function(xTrain,xTest,k){
  return(FNN::get.knnx(data=xTrain,query = xTest,k=k)$nn.index)
}

# computes the profile of a density
profile_density <- function(t_grid,y_grid,cde_estimate)
{
  g <- rep(NA,length(t_grid))
  for(i in seq_along(t_grid))
  {
    g[i] <- sum(cde_estimate[cde_estimate>t_grid[i]]*(y_grid[2]-y_grid[1]))
  }
  return(g)
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
    
    for(ii in 1:nrow(xTest))
    {
      g_test <- profile_density(t_grid,pred_test$z,
                                pred_test$CDE[ii,])
      neighbors <- which_neighbors(g_train,t(g_test),k=k)
      ths[ii] <- quantile(conformity_score_train[neighbors],probs=alpha)
      #prediction_bands[[ii]] <- pred_test$z[pred_test$CDE[ii,]>=ths]
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
# t_grid is a grid of values for the f(.|x)
dist_split_prediction_bands <- function(cde_fit,
                                        xTrain,yTrain,
                                        xTest,
                                        alpha=0.1,
                                        yTest=NULL)
{
  pred_test <- predict(cde_fit,xTest)
  pred_train <- predict(cde_fit,xTrain)
  
  ths = quantile(cum_dist(pred_train$z,pred_train$CDE,yTrain),
                 probs = c(alpha/2,1-alpha/2))
  
  prediction_bands_which_belong <- list()
  FTest <- matrix(NA,nrow(xTest),length(pred_train$z))
  for (ii in 1:nrow(xTest)){
    FTest[ii,] <- cumsum(pred_test$CDE[ii,])*diff(pred_train$z)[1]
    prediction_bands_which_belong[[ii]] <- FTest[ii,]>=ths[1]&FTest[ii,]<=ths[2]
    
  }
  
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=pred_train$z,FTest=FTest,ths=ths))  
  
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
    yTest_covered[ii] <- fit_dist_split$FTest[ii,which_closest]>=ths[1]&fit_dist_split$FTest[ii,which_closest]<=ths[2]
    
  }
  fit_dist_split$yTest_covered <- yTest_covered
  fit_dist_split$prediction_bands_size <- sapply(fit_dist_split$prediction_bands,
                                                 function(x)
                                                   mean(x))*diff(range(fit_dist_split$y_grid))
  fit_dist_split$prediction_bands <- NULL
  fit_dist_split$FTest <- NULL
  return(fit_dist_split)
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