##########################
# process phenotype data #
##########################

# libraries ----
library(dplyr)
library(readxl)

# process data ----

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

# save data ----
save(pheno, file = "data/pheno.RData")
