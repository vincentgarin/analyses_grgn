################
# GxE analysis #
################

# libraries ----
library(dplyr)
library(readxl)
library(lme4)
library(nlme)
library(SpATS)
library(metan)
library(ggplot2)
library(statgenGxE)
library(statgenSTA)

# load data ----

load(file = "data/pheno.RData")

# selected trait ----
trait_sel <- "PANW"
# within environment BLUE computation (SpATS models) ----

env_id <- unique(pheno$env)

d_env_i <- pheno[pheno$env == env_id[1], ]

m <- SpATS(response = "PANW", genotype = "geno",
           genotype.as.random = TRUE,
           spatial = ~SAP(col, row, nseg = c(20,20)),
           fixed = NULL, random = '~ rep + rep:block + row_f + col_f',
           data = d_env_i,
           control = list(maxit = 50, tolerance = 1e-06, monitoring = 1))
plot(m)
h2 <- getHeritability(m)

# BLUP
pred <- predict(m, which = 'geno')
BLUP <- pred[, c("geno", "predicted.values")]
colnames(BLUP)[2] <- "BLUP"

# BLUE
m <- SpATS(response = trait_sel, genotype = "geno",
           genotype.as.random = FALSE,
           spatial = ~SAP(col, row, nseg = c(20,20)),
           fixed = NULL, random = '~ rep + rep:block + row_f + col_f',
           data = d_env_i,
           control = list(maxit = 50, tolerance = 1e-06, monitoring = 1))
plot(m)

pred <- predict(m, which = 'geno')
BLUE <- pred[, c("geno", "predicted.values")]
colnames(BLUP)[2] <- "BLUE"

# within environment BLUE computation: iteration over all env ----

# Get unique environment IDs
env_ids <- unique(pheno$env)

# Initialize a vector to store heritability values
h2_vector <- numeric(length(env_ids))

# Initialize data.frames to store BLUP and BLUE values
BLUP_df <- data.frame(geno = character(), BLUP = numeric())
BLUE_df <- data.frame(geno = character(), BLUE = numeric())

# Loop over each environment ID
for (i in seq_along(env_ids)) {
  env_id <- env_ids[i]
  
  # Subset data for the current environment
  d_env_i <- pheno[pheno$env == env_id, ]
  
  if(!all(is.na(d_env_i$YIELD))){
    
    # Fit SpATS model for BLUP
    m_blup <- SpATS(
      response = trait_sel,
      genotype = "geno",
      genotype.as.random = TRUE,
      spatial = ~SAP(col, row, nseg = c(20, 20)),
      fixed = NULL,
      random = ~ rep + rep:block + row_f + col_f,
      data = d_env_i,
      control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)
    )
    
    # Plot the model (optional)
    plot(m_blup)
    
    # Get heritability
    h2_vector[i] <- getHeritability(m_blup)
    
    # Predict BLUP
    pred_blup <- predict(m_blup, which = 'geno')
    BLUP_env <- pred_blup[, c("geno", "predicted.values")]
    colnames(BLUP_env)[2] <- "BLUP"
    BLUP_env$env <- env_id  # Add environment ID for reference
    
    # Append BLUP to the BLUP data.frame
    BLUP_df <- rbind(BLUP_df, BLUP_env)
    
    # Fit SpATS model for BLUE
    m_blue <- SpATS(
      response = trait_sel,
      genotype = "geno",
      genotype.as.random = FALSE,
      spatial = ~SAP(col, row, nseg = c(20, 20)),
      fixed = NULL,
      random = ~ rep + rep:block + row_f + col_f,
      data = d_env_i,
      control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)
    )
    
    # Plot the model (optional)
    plot(m_blue)
    
    # Predict BLUE
    pred_blue <- predict(m_blue, which = 'geno')
    BLUE_env <- pred_blue[, c("geno", "predicted.values")]
    colnames(BLUE_env)[2] <- "BLUE"
    BLUE_env$env <- env_id  # Add environment ID for reference
    
    # Append BLUE to the BLUE data.frame
    BLUE_df <- rbind(BLUE_df, BLUE_env)
    
    print(paste0("Env:", as.character(env_id)))
    
  }
  
}

res_h2 <- data.frame(env = env_ids, h2 = h2_vector)
res_h2

# save results
save(BLUE_df, file = file.path("output/pheno/",
                               paste0("within_env_BLUE_", trait_sel, ".RData")))

# save results
save(BLUP_df, file = file.path("output/pheno/",
                               paste0("within_env_BLUP_", trait_sel, ".RData")))


AIC_vec <- rep(NA, 4)
BIC_vec <- rep(NA, 4)

# GxE models: Compound symmetry ----
m_cs <- lme(
  BLUE ~ env,  # Fixed effect of environment
  random = ~ 1 | geno,
  correlation = corCompSymm(form = ~ 1 | geno),
  data = BLUE_df
)

summary(m_cs)

# variance components
VCOV_est <- getVarCov(m_cs)
V_g <- VCOV_est[1, 1]
V_e <- (m_cs$sigma)^2


AIC_vec[1] <- AIC(m_cs)
BIC_vec[1] <- BIC(m_cs)

# GxE models: Diagonal ----
m_diag <- gls(
  BLUE ~ env,  # Fixed effect of environment
  weights = varIdent(form = ~ 1 | env),
  data = BLUE_df
)

summary(m_diag)

# variance components
V_ej <- (c(1.0000000, coef(m_diag$modelStruct$varStruct, unconstrained=F))*m_diag$sigma)^2
env_id <- as.character(unique(BLUE_df$env))
names(V_ej)[1] <- env_id[!(env_id %in% names(V_ej))]

AIC_vec[2] <- AIC(m_diag)
BIC_vec[2] <- BIC(m_diag)

# GxE models: CS + Diag ----
m_cs_diag <- lme(
  BLUE ~ env,  # Fixed effect of environment
  random = ~ 1 | geno,
  correlation = corCompSymm(form = ~ 1 | geno),
  weights = varIdent(form = ~ 1 | env),
  data = BLUE_df
)

summary(m_cs_diag)


# variance components
VCOV_est <- getVarCov(m_cs_diag)
V_g <- VCOV_est[1, 1]

V_ej <- (c(1.0000000, coef(m_cs_diag$modelStruct$varStruct, unconstrained=F))*m_cs_diag$sigma)^2
env_id <- as.character(unique(BLUE_df$env))
names(V_ej)[1] <- env_id[!(env_id %in% names(V_ej))]


AIC_vec[3] <- AIC(m_cs_diag)
BIC_vec[3] <- BIC(m_cs_diag)

# GxE models: unstructured ----
m_un <- lme(
  BLUE ~ env,  # Fixed effect of environment
  random = list(geno = pdSymm(form = ~ -1 + env)),
  data = BLUE_df
)

summary(m_un)

# variance components
VCOV_est <- getVarCov(m_un)
n_env <- length(unique(m_un$data$env))
V_g <- VCOV_est[1:n_env, 1:n_env]
V_e <- (m_un$sigma)^2

AIC_vec[4] <- AIC(m_un)
BIC_vec[4] <- BIC(m_un)


# GxE models: selection ----
model_id <- c("cs", "diag", "cs + diag", "un")
res_AIC <- data.frame(model = model_id, AIC = AIC_vec)

model_id <- c("cs", "diag", "cs + diag", "un")
res_BIC <- data.frame(model = model_id, BIC = BIC_vec)

# data processing for GGE, AMMI and Finlay Wilkinson models ----
TD <- createTD(data = pheno, genotype = "geno", trial = "env",
               loc = "location",
               repId = "rep", 
               trLat = "lat", 
               trLong = "lon",
               subBlock = "block",
               rowCoord = "row",
               colCoord = "col")

# specify the design
meta <- getMeta(TD = TD)
meta$trDesign <- rep("res.ibd", length(TD))
TD <- setMeta(TD = TD, meta = meta)

m_ST <- fitTD(TD = TD, traits = trait_sel, engine = "SpATS")

plot(m_ST, 
     plotType = "base", 
     what = "random")

# prediction intra-env
BLUE <- extractSTA(STA = m_ST, what = "BLUEs")
BLUP <- extractSTA(STA = m_ST, what = "BLUPs")
colnames(BLUE)[3] <- colnames(BLUP)[3] <- "trait"

TDGxE <- STAtoTD(STA = m_ST,
                 what = c("BLUEs", "seBLUEs",
                          "BLUPs", "seBLUPs"))

# GxE models: Finlay-Wilkinson ----
FW <- gxeFw(TD = TDGxE, trait = paste0("BLUEs_", trait_sel))
summary(FW)
d_geno <- FW$estimates
d_env <- FW$envEffs

# Create a prediction grid: each genotype × each environment value
plot_df <- d_geno %>%
  tidyr::expand_grid(EnvEff = d_env$EnvEff) %>%
  mutate(
    Pred = GenMean + Sens * EnvEff   # Finlay–Wilkinson regression line
  )

# plot
p <- ggplot(plot_df, aes(x = EnvEff, y = Pred, color = Genotype)) +
  geom_line(size = 0.5) +
  # optional: add the actual environmental mean point locations
  geom_point(data = d_env, aes(x = EnvEff, y = EnvEff), 
             inherit.aes = FALSE, shape = 4, size = 3) +
  labs(
    x = "Environmental mean",
    y = "Predicted genotype performance",
    title = "Finlay–Wilkinson Regression"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

p

# plot(FW, plotType = "line")
plot(FW, plotType = "scatter")
# GxE models: GGE ----
gge_model <- gge(BLUE, trial, genotype, trait)
plot(gge_model)

# GxE models: AMMI ----

AMMI <- gxeAmmi(TD = TDGxE, trait = paste0("BLUEs_", trait_sel))
summary(AMMI)
plot(AMMI, plotType = "AMMI2", scale = 0.5)
