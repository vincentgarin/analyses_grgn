######################################
# Prediction génomique - main effect #
######################################

# librairies
library(mppR)
library(rrBLUP)
library(sommer)
library(glmnet)
library(ggplot2)
library(BGLR)
library(randomForest)
library(e1071)

# ad-hoc functions
source("./functions/GP_random_forest.R")
source("./functions/GP_RF_tune.R")

# data (à modifier)
# load(file = "./data/GP_models_geno.RData")
# load(file = "./data/GP_models_pheno.RData")
# load(file = "./data/GP_models_map.RData")

# process pheno data
GID <- rownames(geno)
map <- map[, -3]

# calculate kinship
K <- A.mat(X = geno)

d <- data.frame(GID = rownames(pheno), y = pheno[, 2])
pheno <- pheno[, 2, drop = FALSE]

k <- 5
n_rep <- 2

# space to store the results
methods_vec <- c("GLUP", "Ridge", "LASSO", "Enet",
                 "BayesA", "BayesB", "BayesC", "RF")

n_meth <- length(methods_vec)
res_list <- vector(mode = "list", length = n_meth)

for(m in 1:n_meth){
  res_list[[m]] <- cor_res <- matrix(NA, nrow = k, ncol = n_rep)
}
names(res_list) <- methods_vec 

# CV loop ----

set.seed(45683)

for(j in 1:n_rep){
  
  print(paste0("r = ", j))
  # partition the data (TS, VS)
  VS <- matrix(sample(1:nrow(pheno)), ncol = k)
  
  for(i in 1:k){
    
    # prepare data
    d_i <- d # copy the original data
    d_res <- d
    
    # mask validation set data
    d_i$y[VS[, i]] <- NA
    
    # calculate prediction model
    # GBLUP ----
    m1 <- kin.blup(data = d_i, geno = "GID", pheno = "y", K = K, GAUSS = FALSE)
    
    # extract BLUP and add to datset checking for order
    d_res$y_pred1 <- m1$g[d$GID]
    
    print(paste0("GBLUP, k=", i))
    
    # Ridge regression ----
    X_TS <- geno[-VS[, i], ]
    y_TS <- d$y[-VS[, i]]
    
    # remove NA values
    NA_pos <- is.na(y_TS)
    X_TS <- X_TS[!NA_pos, ]
    y_TS <- y_TS[!NA_pos]
    
    m3_cv <- cv.glmnet(x = X_TS, y = y_TS, family = "gaussian",
                       alpha = 0, nfolds = 10)
    
    m3_fit <- glmnet(x = X_TS, y = y_TS, family = "gaussian",
                     alpha = 0, lambda = m3_cv$lambda.min)
    B_m3 <- m3_fit$beta
    
    # calculate BLUP (X * B) add to datset checking for order
    u3 <- as.matrix(geno %*% B_m3)
    d_res$y_pred3 <- u3[d$GID, 1]
    
    print(paste0("Ridge reg., k=", i))
    
    # LASSO ----
    
    m4_cv <- cv.glmnet(x = X_TS, y = y_TS, family = "gaussian",
                       alpha = 0.5, nfolds = 10)
    
    m4_fit <- glmnet(x = X_TS, y = y_TS, family = "gaussian",
                     alpha = 0.5, lambda = m4_cv$lambda.min)
    B_m4 <- m4_fit$beta
    
    # calculate BLUP (X * B) add to datset checking for order
    u4 <- as.matrix(geno %*% B_m4)
    d_res$y_pred4 <- u4[d$GID, 1]
    
    print(paste0("LASSO, k=", i))
    
    # Enet ----
    m5_cv <- cv.glmnet(x = X_TS, y = y_TS, family = "gaussian",
                       alpha = 1, nfolds = 10)
    
    m5_fit <- glmnet(x = X_TS, y = y_TS, family = "gaussian",
                     alpha = 1, lambda = m5_cv$lambda.min)
    B_m5 <- m5_fit$beta
    
    # calculate BLUP (X * B) add to datset checking for order
    u5 <- as.matrix(geno %*% B_m5)
    d_res$y_pred5 <- u3[d$GID, 1]
    
    print(paste0("Elastic-net, k=", i))
    
    # BayesA ----
    mBA <- BGLR(y=d_i$y,ETA=list( list(X=geno,model='BayesA')), 
                nIter=1000,burnIn=100,saveAt="./results",
                verbose = FALSE)
    
    # add y_hat
    d_res$y_pred6 <- mBA$yHat
    
    print(paste0("BayesA, k=", i))
    
    # BayesB ----
    mBB <- BGLR(y=d_i$y,ETA=list(list(X=geno,model='BayesB', probIn = 0.33)), 
                nIter=1000,burnIn=100,saveAt="./results",
                verbose = FALSE)
    
    # add y_hat
    d_res$y_pred7 <- mBB$yHat
    
    print(paste0("BayesB, k=", i))
    
    # BayesCpi ----
    mBC <- BGLR(y=d_i$y,ETA=list(list(X=geno,model='BayesC')), 
                nIter=1000,burnIn=100,saveAt="./results",
                verbose = FALSE)
    
    # add y_hat
    d_res$y_pred8 <- mBC$yHat
    
    print(paste0("BayesCpi, k=", i))
    
    # Random forest ----
    
    # option to make a GP_random_forest_tune function
    # via tune.randomForest
    
    # use the same data partition as Ridge, LASSO and Enet
    
    # mRF <- GP_random_forest(x = X_TS, y = y_TS, x_vs = geno)
    # fine_tuning of RF take some time (Lourenco, parallelize)
    mRF <- GP_RF_tune(x = X_TS, y = y_TS, x_vs = geno, mtry = c(1000, 3000),
                      ntree = c(50, 100), nodesize = 5, k_cross = 3)
    # add y_hat
    d_res$y_pred9 <- mRF$y_hat[d$GID]
    
    print(paste0("Random forest, k=", i))
    
    
    # process the results ----
    
    colnames(d_res) <- c("GID","y", methods_vec)
    
    # calculate the correlation (obs, pred)
    d_vs <- d_res[VS[, i], ] # selection only VS data
    
    # cor inter method
    round(cor(d_vs[, 2:ncol(d_vs)], use = "complete.obs"), 3)
    
    
    # print(plot(x = d_vs$y, y = d_vs$y_pred, xlab = "obs", ylab = "pred",
    #            main = paste("Rep", j, paste0("k=",i))))
    
    for(m in 1:n_meth){
      res_list[[m]][i, j] <- cor(x = d_vs$y, y = d_vs[, m + 2],
                                 use = "complete.obs")
    }
    
    print(paste0("k = ", i))
  }
}


# plot the results ----

d <- c()
for(i in 1:length(res_list)){
  d_i <- data.frame(model = names(res_list)[i], cor = c(res_list[[i]]))
  d <- rbind(d, d_i)
}

pl <- ggplot(d, aes(x = model, y = cor)) +
  geom_boxplot() + 
  labs(title = "GP model comp", x = "models", y = "correlation(y_obs, y_pred)")

pl