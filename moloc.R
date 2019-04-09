## load moloc and read in datasets
library("moloc")
library("mvtnorm")
library("data.table")
library("dplyr")
library("TwoSampleMR")

## dataframes for each trait analysed
GWAS <- as.data.frame(data_genome[1])
eQTL <- as.data.frame(data_genome[2])
mQTL <- as.data.frame(data_genome[3])

## read in reg file
reg_file <- read.table('moloc_reg_file.txt', header=T)

## Create bed file with combination of ProbeIDs
library(GenomicRanges)
methyl_table <- DT[, list(CHR= unique(CHR), START = min(POS), STOP = max(POS)), by = ProbeID]
bed1.gr <- GRanges(seqnames = bed$CHR,IRanges(start = bed$START, end= bed$STOP))
bed2.gr <- GRanges(seqnames = methyl_table$CHR,IRanges(start =methyl_table$START, end= methyl_table$STOP))
my.overlap <- findOverlaps(query = bed2.gr, subject = bed1.gr)
bed = cbind(bed[my.overlap@subjectHits,], methyl_table[my.overlap@queryHits,])
bed = bed[,c(1,6,7,8,9,10)]
## Setting the prior_var to "default" uses coloc defaults
coloc.genome(listData = data_genome, bed, cores=1, have_alleles=TRUE, bychrpos=TRUE, prior_var=NULL, priors=c(1e-04, 1e-06, 1e-07), min_nsnps = 50, write=FALSE, outfolder = "test")
