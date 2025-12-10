####################
# GP_random_forest #
####################

GP_random_forest <- function(x, y, x_vs, ntree=500, mtry=NULL, 
                             nodesize=1,
                             importance=TRUE,
                             importanceSD=TRUE,
                             mse=TRUE,
                             rsq=TRUE,
                             type="regression"){
  
  if(is.null(mtry)){
    mtry <- floor(ncol(x)/3)
  }
  
  rf <- randomForest(x = x, y = y, ntree=ntree, mtry=mtry, nodesize=1,
                     importance=importance, importanceSD=importanceSD,
                     mse=mse, rsq=rsq, type=type, predicted=TRUE)
  
  # calculate y_hat (prediction)
  y_hat <- predict(rf, x_vs, type="response", predict.all=FALSE,
                   proximity=FALSE)
  
  if(importance){
    Bj_imp <- importance(rf)
  } else {
    Bj_imp <- NULL
  }
  
  return(list(m = rf, y_hat = y_hat, Bj = Bj_imp))
  
}