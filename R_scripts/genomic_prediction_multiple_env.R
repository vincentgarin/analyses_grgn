##########################################
# genomic prediction - multi-environment #
##########################################

# library ----
library(stringr)
library(sommer)
library(BGGE)
library(BGLR)
library(reshape2)
library(tidyr)
library(ggplot2)

# ad-hoc functions ----
source("functions/fct_GPGxE_model_computation.R")
source("functions/utils_BGGE_model.R")
source("functions/utils_BGLR_model.R")
source("functions/utils_som_model.R")
source("functions/utils_unsBLUP.R")
source("functions/utils_extract_BGGE_estimates.R")
source("functions/utils_extract_som_estimates.R")
source("functions/utils_nm_0_lead.R")
source("functions/utils_plot_cor_within_env.R")
source("functions/utils_cor_within_env.R")

# set seed ----
set.seed(68711)

# load data ----
load(file = "data/geno/geno_012.RData")
load(file = "output/pheno/within_env_BLUE_YIELD.RData")

# process data ----
pheno <- BLUE_df
pheno <- pheno %>% pivot_wider(names_from = "env", values_from = "BLUE") %>%
  as.data.frame()

rm(BLUE_df)
rownames(pheno) <- pheno$geno
pheno <- pheno[, -1]

# harmonize genotype list
geno_cm <- intersect(rownames(pheno), rownames(geno))

pheno <- pheno[geno_cm, ]
pheno <- as.matrix(pheno)
geno <- geno[geno_cm, ]

# kinship calculation
K <- A.mat(X = geno - 1)

# single models computation ----
m1 <- GPGxE_model_computation(model = "BGGE_MM", pheno = pheno, K = K)
m2 <- GPGxE_model_computation(model = "BGGE_MDe", pheno = pheno, K = K)
m3 <- GPGxE_model_computation(model = "BGLR_DIAG", pheno = pheno, K = K)
m4 <- GPGxE_model_computation(model = "BGLR_UN", pheno = pheno, K = K)
m5 <- GPGxE_model_computation(model = "BGLR_FA", pheno = pheno, K = K)

pred1 <- m1$y_hat
pred2 <- m2$y_hat
pred3 <- m3$y_hat
pred4 <- m4$y_hat
pred5 <- m5$y_hat

# test set definition:: Obs Geno Unobserved Env ----
scen <- "Obs. Geno Unobs. Env"
G <- rownames(pheno)
E <- colnames(pheno)

G_TR <- G_TS <- G
E_TS <- E[2]
E_TR <- E[!(E %in% E_TS)]

# copy of the phenotype data
pheno_s <- pheno

# mask the values to be predicted (test set)
pheno_s[, (E %in% E_TS)] <- NA

# remove the whole environment since completely missing
pheno_s <- pheno_s[, E_TR]

# model computation
m <- GPGxE_model_computation(model = "BGLR_FA", pheno = pheno_s, K = K)

# pred
pred <- m$g_hat
rownames(pred) <- rownames(pheno_s)
colnames(pred) <- colnames(pheno_s)
obs <- pheno[, E_TS]
pred <- rowMeans(pred, na.rm = TRUE)
d_obs_pred <- data.frame(obs = obs, pred = pred)

plot(x = d_obs_pred$pred, y = d_obs_pred$obs)
p <- ggplot(d_obs_pred, aes(x = pred, y = obs)) +
  geom_point() +                     # Add scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = scen) +   # Add title
  theme_minimal()
p
cor(x = d_obs_pred$pred, y = d_obs_pred$obs, use = "complete.obs")

# test set definition:: Unobserved Geno Observed Env ----
scen <- "Unobs. Geno Obs. Env"
G <- rownames(pheno)
E <- colnames(pheno)

geno_sel_ind <- sample(1:length(G), size = round(length(G)/5))
G_TS <- G[geno_sel_ind]
G_TR <- G[!(G %in% G_TS)]
E_TS <- E_TR <- E

# copy of the phenotype data
pheno_s <- pheno

# mask the values to be predicted (test set)
pheno_s[(G %in% G_TS), ] <- NA

# model computation
m <- GPGxE_model_computation(model = "BGLR_FA", pheno = pheno_s, K = K)

# pred
pred <- m$g_hat
rownames(pred) <- rownames(pheno_s)
colnames(pred) <- colnames(pheno_s)

pred <- pred[G_TS, ]
pred$Genotype <- rownames(pred)
pred <- pred %>%
  pivot_longer(
    cols = -Genotype,          # All columns except 'geno' will be pivoted
    names_to = "Env",      # New column for environment names
    values_to = "pred"     # New column for BLUE values
  )

# obs
obs <- data.frame(pheno[G_TS, ])
obs$Genotype <- rownames(obs)
obs <- obs %>%
  pivot_longer(
    cols = -Genotype,          # All columns except 'geno' will be pivoted
    names_to = "Env",      # New column for environment names
    values_to = "obs"     # New column for BLUE values
  )

d_obs_pred <- merge(x = obs, y = pred, by = c("Genotype", "Env"))
colnames(d_obs_pred)[3] <- "trait" 
cor_wth_env <- cor_within_env(d_obs_pred)
print(cor_wth_env)

plot_cor_within_env(d = d_obs_pred, title = scen)

# test set definition:: Unobserved Geno Unobserved Env ----
G <- rownames(pheno)
E <- colnames(pheno)

geno_sel_ind <- sample(1:length(G), size = round(length(G)/3))
G_TS <- G[geno_sel_ind]
G_TR <- G[!(G %in% G_TS)]
E_TS <- E[5]
E_TR <- E[!(E %in% E_TS)]

# copy of the phenotype data
pheno_s <- pheno

# mask the values to be predicted (test set)
pheno_s[(G %in% G_TS), (E %in% E_TS)] <- NA

# model computation
m <- GPGxE_model_computation(model = "BGLR_FA", pheno = pheno_s, K = K)

# pred
pred <- m$g_hat
rownames(pred) <- rownames(pheno_s)
colnames(pred) <- colnames(pheno_s)

pred <- pred[G_TS, E_TS]

# obs
obs <- pheno[(G %in% G_TS), E_TS]

d_obs_pred <- data.frame(trait = obs, pred = pred, Env = E_TS)
cor_wth_env <- cor_within_env(d_obs_pred)
print(cor_wth_env)

plot_cor_within_env(d = d_obs_pred, title = scen)
