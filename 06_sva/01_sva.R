# 06 Surrogate Variable Analysis (SVA)

# Reference: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

setwd("~/bio/code/mpip/dex-stim-human-dna-methyl-qc/")

library(sva)
library(lumi) # convert beta values to M-values 

input.parameters.fn <- "input_parameters.csv"
input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# 1. Setting up the data from an ExpressionSet
x     <- load(pd_clean.fn)
pheno <- get(x)

x        <- load(beta.combat.expr.set.fn) #Betas_combated_ExprSet
expr.set <- get(x)
rm(x)

# beta.mtrx.fn <- "~/bio/datasets/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
# beta.mtrx    <- readRDS(beta.mtrx.fn)
# pheno  <- pheno[pheno$Sample_Name %in% colnames(beta.mtrx),]

# pheno <- fread("~/bio/datasets/pheno/pheno_full_for_kimono.csv", na.strings = "#N/A")
# pheno <- pheno[!is.na(DNAm_ID)]

#2. Convert beta values to M-values

m.mtrx <- beta2m(beta.mtrx)

# 3. Create the full and null model matrices
mod  <- model.matrix(~ Sample_Group, data = pheno)
mod0 <- model.matrix(~ 1, data = pheno)

# 4. Applying the 'sva' to estimate surrogate variables
n.sv  <- num.sv(m.mtrx, mod, method="leek")
svobj <- sva(m.mtrx, mod, mod0, n.sv = n.sv)



