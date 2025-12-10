##############
# GP_RF_tune #
##############

# function to fine tune the hyper parameter of a random forest model.

GP_RF_tune <- function(x, y, x_vs, ntree=c(100, 500, 1000),
                       mtry=NULL, 
                       nodesize=c(1, 5, 10),
                       mse=TRUE,
                       rsq=TRUE,
                       type="regression",
                       k_cross = 5){
  
  if(is.null(mtry)){
    Np <- ncol(x)
    mtry <- round(c(0.3, 0.6, 0.9) * Np) 
  }
  
  rf <- tune.randomForest(x = x, y = y, ntree=ntree, mtry=mtry,
                          nodesize=nodesize, mse=mse, rsq=rsq,
                          type=type,
                          tunecontrol=tune.control(cross = k_cross))
  
  # calculate y_hat (prediction)
  y_hat <- predict(rf$best.model, x_vs, type="response", predict.all=FALSE,
                   proximity=FALSE)
  
  return(list(m = rf$best.model, y_hat = y_hat))
  
}