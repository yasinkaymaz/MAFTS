#dev script.
library(tidyverse)
#tail -n+2 ~/Documents/MedGenetics/BioinfoAnalysis/gatk/DNA001r1_snpEff_genes.txt |sed 's/#//g' >  ./DNA001r1_snpEff_genes.txt
snpeff <- read.delim("data/DNA001r1_snpEff_genes.txt")
snpeff
snpeff[snpeff$variants_effect_missense_variant > 0,]
names(snpeff)
snpeff[, c("X.GeneName", "TranscriptId", "variants_effect_missense_variant", "variants_effect_synonymous_variant")]


snpeff[which(snpeff_data$GeneName == "BRAF"), c(rep(TRUE, 4), colSums(snpeff[which(snpeff_data$GeneName == "BRAF"),-c(1:4)])!=0)]
c(rep(TRUE, 4), colSums(snpeff[which(snpeff_data$GeneName == "BRAF"),-c(1:4)])!=0)

sampleSheet <- readxl::read_xlsx("data/sampleSheet.xlsx",col_names = T)
sampleSheet %>% filter(SampleID == "DNA004r2")
