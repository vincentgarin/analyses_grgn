########################
# process genetic data #
########################

library(vcfR)

vcf <- read.vcfR("data/geno/Dart_ceraas.vcf")

# Extract genotypes (CORRECT object)
gt <- extract.gt(vcf, element = "GT")

# define parents
parents <- colnames(gt)[1:2]

parents_gt <- gt[, parents]

# Parent polymorphism filter
keep_parents <- apply(parents_gt, 1, function(x) {
  !any(is.na(x)) && x[1] != x[2]
})

vcf_sub <- vcf[keep_parents, ]

# Re-extract GT after filtering
gt_sub <- extract.gt(vcf_sub, element = "GT")

# Missing data per SNP
miss_snp <- rowMeans(is.na(gt_sub) | gt_sub == "./.")

# Minor allele frequency
maf_vals <- maf(vcf_sub)

# Combine filters
keep_qc <- (maf_vals[, 4] > 0.01) & (miss_snp < 0.25)

vcf_sub <- vcf_sub[keep_qc, ]

# transform into 012 format

gt <- extract.gt(vcf_sub, element = "GT")
maf_vals <- maf(vcf_sub)

gt_num <- matrix(NA_integer_, nrow = nrow(gt), ncol = ncol(gt))
rownames(gt_num) <- rownames(gt)
colnames(gt_num) <- colnames(gt)

gt_num[gt == "0/0" | gt == "0|0"] <- 0
gt_num[gt %in% c("0/1","1/0","0|1","1|0")] <- 1
gt_num[gt == "1/1" | gt == "1|1"] <- 2
flip <- maf_vals[, 4] > 0.5
gt_minor <- gt_num
gt_minor[flip, ] <- 2 - gt_num[flip, ]

summary(as.vector(gt_minor))
dim(gt_minor)

# arrange and remove the parents
geno_012 <- t(gt_minor)
geno_012 <- geno_012[-c(1:2), ]

# input missing values with a random value
for(i in 1:ncol(geno_012)){
  m_i <- geno_012[, i]
  if(sum(is.na(m_i)) > 0){
    m_i[is.na(m_i)] <- sample(x = c(0, 1, 2), size = sum(is.na(m_i)), replace = TRUE) 
    geno_012[, i] <- m_i
    }
  
}

sum(is.na(c(geno_012)))

# modify geno names
geno_names <- rownames(geno_012)
geno_names <- sub("_.*", "", geno_names)
rownames(geno_012) <- geno_names

geno <- geno_012
# save the data
save(geno, file = "data/geno/geno_012.RData")
