#################
# Data analysis #
#################

# libraries ----
library(dplyr)
library(ggplot2)
library(lme4)
library(readxl)
library(desplot)

# load data ----
trial_files <- list.files(path = "data",
                          pattern = "Trait_Data.xlsx$", full.names = FALSE)

trial_id <- substr(x = trial_files, start = 1, stop = 4)

pheno <- c()

for(i in 1:length(trial_files)){
  trial_i <- read_xlsx(path = file.path("data", trial_files[i]))
  trial_i$Env <- trial_id[i]
  
  # possibility to add location and year variable
  
  pheno <- bind_rows(pheno, trial_i)
  
}

# process data ----

# reposition columns
pheno <- pheno %>% relocate(Env, .before = EXPT_DESIGN)

# select needed columns
pheno <- pheno %>% select(Env, LOCATION_NAME, SITE_LAT, SITE_LONG, EXPT_DESIGN,
                          PLOT_NO, REP_NO, BLOCK_NO, ROW, COL, FIELDMAP.RANGE,
                          FIELDMAP.COLUMN, DESIGNATION,
                          FLfL_C_day, Flo_C_day, PcleLng_M_cm,
                          PcleWid_M_cm, PcleDMYld_C_gPlot, GHvYld_C_kgha)

# rename columns (for interpretability)
colnames(pheno) <- c("env", "location", "lat", "lon", "design",
                     "plot_no", "rep", "block", "row", "col", "field_range",
                     "field_column", "geno", "FLAG", "FLOWER", "PANL", "PANW",
                     "PANDMY", "YIELD")


# convertir les variables
pheno$geno <- as.factor(pheno$geno)
pheno$env <- as.factor(pheno$env)
pheno$rep <- as.factor(pheno$rep)
pheno$block <- as.factor(pheno$block)
pheno$row_f <- as.factor(pheno$row)
pheno$col_f <- as.factor(pheno$col)

# Inspection générale des essais ----
desplot(data = pheno, form = FLOWER ~ row * col | env)

# Inspecter les variables
hist(pheno$FLOWER)

# Charger la bibliothèque ggplot2
ggplot(pheno, aes(x = FLOWER)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_wrap(~ env, scales = "free_y") +
  labs(title = "Histogramme de Floraison par environnement",
       x = "Floraison",
       y = "Fréquence") +
  theme_minimal()


# phenotypic model computation ----
m <- lmer(formula = FLOWER ~ env + (1|env:rep) + (1|env:rep:block) +
            (1|env:row) + (1|env:col) + (1|geno) + (1|geno:env),
          REML = TRUE, data = pheno)

summary(m)

# heritability computation ----
remat <- summary(m)
Verr <- remat$sigma^2
V_var <- remat$varcor
Vg <- V_var$geno[1, 1]
Vgxe <- V_var$`geno:env`[1, 1]
n_env = 6
n_rep = 2
h2 <- (Vg)/(Vg + (Vgxe/n_env) + (Verr/(n_rep * n_env)))
h2

# Adjusted means: BLUP ----
BLUP_geno <- ranef(m, whichel = "geno")
BLUP_geno <- BLUP_geno$geno



BLUP <- data.frame(geno = rownames(BLUP_geno),
                            BLUP = BLUP_geno[, 1])

# Adjusted means: BLUE ----
m <- lmer(formula = FLOWER ~ env + (1|env:rep) + (1|env:rep:block) +
            (1|env:row) + (1|env:col) + geno + (1|geno:env),
          REML = TRUE, data = pheno)

m_sum <- summary(m)
coeff <- m_sum$coefficients

Int <- coeff[rownames(coeff) == '(Intercept)', 1]
BLUE <- coeff[grepl(pattern = "geno", x = rownames(coeff)), ]
geno_id <- gsub(pattern = "geno", replacement = "", x = rownames(BLUE))

BLUE <- data.frame(geno = geno_id, BLUE = BLUE[, 1] + Int)

# add the reference value
d_used <- m@frame
geno_ref <- levels(d_used$geno)[1]
stopifnot(!geno_ref %in% BLUE$geno)
BLUE <- rbind(data.frame(geno = geno_ref, BLUE = Int), BLUE)

d_adj_means <- merge(BLUP, BLUE, by = "geno")

# comparison BLUE et BLUP
plot(x = d_adj_means$BLUP, y = d_adj_means$BLUE)
cor(x = d_adj_means$BLUP, y = d_adj_means$BLUE)

