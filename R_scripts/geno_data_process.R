#Load packages----

library(tidyverse)
library(kableExtra)
#devtools::install_github("wolfemd/genomicMateSelectR", ref = 'master') 
library(genomicMateSelectR)
library(AGHmatrix)
library(ComplexHeatmap)

geno <- read.table("data/geno.hmp.txt", header = T, check.names = F)#import geno data

#recode geno (0,1,2)----

geno <- geno |>
  mutate(across(4:ncol(geno), ~{
    case_when(
      . == paste(reference, reference, sep = "") ~ 2,
      . == paste(recess, recess, sep = "") ~ 0,
      TRUE ~ 1
    )
  })) |>
  select(-c(reference, recess)) |> 
  column_to_rownames(var ="rs") |> 
  t()

geno[1:5,1:5]

#save data
save(geno, file="output/geno/geno_012.RData")
#Stop here



#No RILs & markers (Table)

geno|>
  dim() |>
  t() |>
  kbl(escape = F, align = 'c',
      col.names = c("Numbers of RIL", "Number of markers")) |>
  kable_classic("hover", full_width = F, position = "center", fixed_thead = T)

geno_data <- cbind(ID = rownames(geno), geno)

write.table(geno_data, "",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

#Data filtering-----

#mac@MacBookThierry Downloads % /Users/thierry/PopLDdecay/bin/PopLDdecay -InVCF progenies.vcf -MaxDist 10000 -MAF 0.05 -Miss 0.2 -OutStat progenies.stat
#mac@192 Downloads % vcftools --vcf progenies.vcf --maf 0.01 --non-ref-ac 1 --max-missing 0.25 --recode --out finalez_filtered

geno_data<- maf_filter(geno, thresh = 0.01) #MAF 1%


#We will use the AGHmatrix package (Amadeu et al., 2016) to build the G matrix
G_matrix = Gmatrix(geno, method = "VanRaden", ploidy = 2, missingValue = NA)#

dim(G_matrix)

Heatmap(G_matrix, show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(title = "Res"))#We can represent this matrix using a heatmap. 





