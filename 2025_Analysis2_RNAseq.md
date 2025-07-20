# Kallisto and DeSeq2
```
## RNA-seq analysis with DESeq2
## Adapted from Stephen Turner, @genetics_blog
library(DESeq2)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("gplots")
# BiocManager::install("DESeq2")
library(apeglm)
# RNA-seq data from PRJNA315516
# https://www.ncbi.nlm.nih.gov//bioproject/PRJNA315516.
# 3 control samples ("ctl"), 3 samples grown under drought condition ("dro")
library(writexl)
# Import & pre-process ----------------------------------------------------

# Here https://support.bioconductor.org/p/105964/#105969
# it says to use the counts.matrix file (and not the TMM.EXPR.matrix file)

# Import data from featureCounts
## Previously ran at command line something like this:
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_dmrt1_knockouts/2023_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_Analysis2")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_dmrt1_knockouts/2023_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_Analysis2"
list.files(dir)


# MF wildtype ccdc ----
# counts <- read.table("MF_NEW.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
counts <- read.table("MF_NEW.isoform.counts.matrix", header=T, row.names = 1)
colnames(counts)
new_counts <- counts[,c(1:2,9:14)]
colnames(new_counts)
coldata <- read.table("MF_sample_genotype.txt", header=T, row.names = 1)
new_coldata <- coldata[c(1:2,9:14),]; new_coldata
new_coldata$batch <- as.factor(new_coldata$batch)
new_coldata$sex <- as.factor(new_coldata$sex)

dds <- DESeqDataSetFromMatrix(countData = round(new_counts),
                              colData = new_coldata,
                              design= ~ sex)
# save the unfiltered log2FoldChange to a dataframe 
dds <- DESeq(dds)
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
res <- results(dds)
MF_ccdc_unfiltered <- res;MF_ccdc_unfiltered
# Only sex related
sex_related_MF_ccdc <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_MF_ccdc, file="Sex_related_MF_ccdc_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
dim(dds)
#[1] 32161     8
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:14,];p
write.csv(p, file="MF_Kallisto_ccdc_DE_DeSeq2.csv", row.names = T)
# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./MF_Kallisto_ccdc_rsquare.xls")


# MF wildtype dmrt1L ----
colnames(counts)
new_counts <- counts[,c(3:5,15:19)]
colnames(new_counts)
coldata <- read.table("MF_sample_genotype.txt", header=T, row.names = 1)
new_coldata <- coldata[c(3:5,15:19),]; new_coldata
new_coldata$batch <- as.factor(new_coldata$batch)
new_coldata$sex <- as.factor(new_coldata$sex)
dds <- DESeqDataSetFromMatrix(countData = round(new_counts),
                              colData = new_coldata,
                              design= ~ sex)
# save the unfiltered log2FoldChange to a dataframe 
dds <- DESeq(dds)
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
res <- results(dds)
MF_dmrt1L_unfiltered <- res;MF_dmrt1L_unfiltered
# Only sex related
sex_related_MF_dmrt1L <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_MF_dmrt1L, file="Sex_related_MF_dmrt1L_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
dim(dds)
#[1] 31961     8
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:73,];p
write.csv(p, file="MF_Kallisto_dmrt1L_DE_DeSeq2.csv", row.names = T)
# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./MF_Kallisto_dmrt1L_rsquare.xls")


# MF wildtype dmrt1S ---- 
colnames(counts)
new_counts <- counts[,c(6:8,20:28)]
colnames(new_counts)
coldata <- read.table("MF_NEW_sample_genotype.txt", header=T, row.names = 1)
new_coldata <- coldata[c(6:8,20:28),]; new_coldata
new_coldata$batch <- as.factor(new_coldata$batch)
new_coldata$sex <- as.factor(new_coldata$sex)
colnames(new_counts)<-rownames(new_coldata)

#############
dds <- DESeqDataSetFromMatrix(countData = round(new_counts),
                              colData = new_coldata,
                              design= ~ batch + sex)
# save the unfiltered log2FoldChange to a dataframe 
dds <- DESeq(dds)
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
res <- results(dds)
MF_dmrt1S_unfiltered <- res;MF_dmrt1S_unfiltered
# Only sex related
sex_related_MF_dmrt1S <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_MF_dmrt1S, file="Sex_related_MF_dmrt1S_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
dim(dds)
#[1] 31881    12 # NEW dmrt1S
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F

# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:150,];p
write.csv(p, file="MF_Kallisto_dmrt1S_DE_DeSeq2.csv", row.names = T)
# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./MF_Kallisto_dmrt1S_rsquare.xls")

##################
###
###   DMRT1L knockout comparisons
###
##################


# wtko dmrt1L ----
# load count data (from Kalisto)
# read the count data from Kallisto that was combined from each sample
# into a single file into a dataframe called "dds"
# dmrt1L_counts <- read.table("dmrt1L.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
dmrt1L_counts <- read.table("dmrt1L.isoform.counts.matrix", header=T, row.names = 1)
dim(dmrt1L_counts)
coldata <- read.table("dmrt1L_sample_genotype.txt", header=T, row.names = 1)
coldata$genotype <- as.factor(coldata$genotype)
coldata$genotype <- relevel(coldata$genotype, ref="wt")
colnames(dmrt1L_counts) <- rownames(coldata)

# subset the males and females into separate dataframes
dmrt1L_counts_femalesonly <- dmrt1L_counts[,c(1,4:5,7:8,12:13,15:16)];colnames(dmrt1L_counts_femalesonly)
dmrt1L_counts_malesonly <- dmrt1L_counts[,c(2:3,6,9:11,14,17:20)];colnames(dmrt1L_counts_malesonly)
dmrt1L_counts_ko_MFonly <- dmrt1L_counts[,c(1:7,9:10,14:16)];colnames(dmrt1L_counts_ko_MFonly)

# subset the genotypes
genotypez_femsonly <- coldata[c(1,4:5,7:8,12:13,15:16),];genotypez_femsonly
genotypez_femsonly$genotype <- relevel(genotypez_femsonly$genotype, ref="wt")
genotypez_malesonly <- coldata[c(2:3,6,9:11,14,17:20),];genotypez_malesonly
genotypez_malesonly$genotype <- relevel(genotypez_malesonly$genotype, ref="wt")
genotypez_ko_MFonly <- coldata[c(1:7,9:10,14:16),];genotypez_ko_MFonly
genotypez_ko_MFonly$genotype <- relevel(genotypez_ko_MFonly$genotype, ref="wt")




# dmrt1L females ko vs wt ----
dds <- DESeqDataSetFromMatrix(countData = round(dmrt1L_counts_femalesonly),
                              colData = genotypez_femsonly,
                              design= ~ genotype)

# save the unfiltered log2FoldChange to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
wtkofemsonly_dmrt1L_unfiltered <- res;wtkofemsonly_dmrt1L_unfiltered
# Only sex related
sex_related_wtkofemsonly_dmrt1L <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_wtkofemsonly_dmrt1L, file="Sex_related_wtkofemsonly_dmrt1L_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# get normalized counts as detailed here: https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
normalized_countz <- counts(dds, normalized=T)
# Get sex related countz
sex_related_wtkofemsonly_dmrt1L_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(normalized_countz)), ]
write.csv(sex_related_wtkofemsonly_dmrt1L_normalized_countz, file="Sex_related_wtkofemsonly_dmrt1L_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)
#temp <- normalized_countz[grepl("gnl\\|XBXL10_1g8729\\|", rownames(normalized_countz)), ];temp

# Get germcell related countz
#germcell_related_wtkofemsonly_dmrt1L_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g35317\\||gnl\\|XBXL10_1g38173\\|", rownames(normalized_countz)), ]
#write.csv(germcell_related_wtkofemsonly_dmrt1L_normalized_countz, file="Germcell_related_wtkofemsonly_dmrt1L_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)



# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
# here we remove transcripts where the average count per sample is 1or more
#keep <- rowSums(counts(dds)) >=length(colnames(dds))
dds <- dds[keep,]
dim(dds)
# [1] 31866     9
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:953,];p
write.csv(p, file="wtkofemsonly_dmrt1L_DE_DeSeq2.csv", row.names = T)

# (from gff3 file)
# stra8.L: XBXL10_1g15025 XBmRNA27986
dmrt1L_counts_femalesonly[grep("XBmRNA27986", rownames(dmrt1L_counts_femalesonly)), ]
res[grep("XBmRNA27986", rownames(res)), ]
wtkofemsonly_dmrt1L_unfiltered[grep("XBmRNA27986", rownames(wtkofemsonly_dmrt1L_unfiltered)), ]

# How many sign DE genes are sex related?
significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L <- resOrdered[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(resOrdered)), ]
significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L_df <- as.data.frame(significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L)
significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L_df <- significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L_df[!(is.na(significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L_df$padj)), ]
significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L_df[(significantDE_that_are_sexrelated_wtkofemsonly_dmrt1L_df$padj < 0.05),]
#                                 baseMean log2FoldChange     lfcSE     stat       pvalue       padj
# gnl|XBXL10_1g34871|XBmRNA65898| 269.1589       1.191215 0.3503818 3.399763 0.0006744426 0.03248846
# this is androgen receptor; this suggests that dmrt1L may downregulate expression of AR in females 
# because it is more highly expressed in female ko compared to female wt...

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./wtkofemsonly_dmrt1L_rsquare.xls")




# dmrt1L males ko vs wt ----
dds <- DESeqDataSetFromMatrix(countData = round(dmrt1L_counts_malesonly),
                              colData = genotypez_malesonly,
                              design= ~ genotype)

# save the unfiltered log2FoldChange to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
wtkomalesonly_dmrt1L_unfiltered <- res;wtkomalesonly_dmrt1L_unfiltered
# Only sex related
sex_related_wtkomalesonly_dmrt1L <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_wtkomalesonly_dmrt1L, file="sex_related_wtkomalesonly_dmrt1L_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# get normalized counts as detailed here: https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
normalized_countz <- counts(dds, normalized=T)
# Get sex related countz
sex_related_wtkomalzonly_dmrt1L_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(normalized_countz)), ]
write.csv(sex_related_wtkomalzonly_dmrt1L_normalized_countz, file="Sex_related_wtkomalzonly_dmrt1L_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)

# Get germcell related countz
#germcell_related_wtkomalzonly_dmrt1L_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g35317\\||gnl\\|XBXL10_1g38173\\|", rownames(normalized_countz)), ]
#write.csv(germcell_related_wtkomalzonly_dmrt1L_normalized_countz, file="Germcell_related_wtkomalzonly_dmrt1L_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)



# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
# here we remove transcripts where the average count per sample is 1or more
#keep <- rowSums(counts(dds)) >=length(colnames(dds))
dds <- dds[keep,]
dim(dds)
# [1] 32042    11
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:37,];p
write.csv(p, file="wtkomalesonly_dmrt1L_DE_DeSeq2.csv", row.names = T)

# How many sign DE genes are sex related?
significantDE_that_are_sex_related_wtkomalesonly_dmrt1L <- resOrdered[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(resOrdered)), ]
significantDE_that_are_sex_related_wtkomalesonly_dmrt1L_df <- as.data.frame(significantDE_that_are_sex_related_wtkomalesonly_dmrt1L)
significantDE_that_are_sex_related_wtkomalesonly_dmrt1L_df <- significantDE_that_are_sex_related_wtkomalesonly_dmrt1L_df[!(is.na(significantDE_that_are_sex_related_wtkomalesonly_dmrt1L_df$padj)), ]
significantDE_that_are_sex_related_wtkomalesonly_dmrt1L_df[(significantDE_that_are_sex_related_wtkomalesonly_dmrt1L_df$padj < 0.05),]
# in males we do not detect any significantly differentially expressed sex related genes in ko compared to wt


# (from gff3 file)
# stra8.L: XBXL10_1g15025 XBmRNA27986
dmrt1L_counts_malesonly[grep("XBmRNA27986", rownames(dmrt1L_counts_malesonly)), ]
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA27986", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# sox9L XBXL10_1g39526 XBmRNA74501
dmrt1L_counts_malesonly[grep("XBmRNA74501", rownames(dmrt1L_counts_malesonly)), ]
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA74501", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# sox9S XBXL10_1g42722 XBmRNA80320
dmrt1L_counts_malesonly[grep("XBmRNA80320", rownames(dmrt1L_counts_malesonly)), ]
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA80320", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# vimentin.L (vim.L) XBmRNA50375 Thumfart.Mansuy.2023
dmrt1L_counts_malesonly[grep("XBmRNA50375", rownames(dmrt1L_counts_malesonly)), ]
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA50375", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# vimentin.S (vim.S) XBmRNA54220 Thumfart.Mansuy.2023
dmrt1L_counts_malesonly[grep("XBmRNA54220", rownames(dmrt1L_counts_malesonly)), ]
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA54220", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# gata1.L XBmRNA64366
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA64366", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# gata1.S XBmRNA72563
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA72563", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# follicle-stimulating hormone receptor (fshr.L) (Thumfart.Mansuy.2023) XBmRNA42008
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA42008", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# follicle-stimulating hormone receptor (fshr.S) (Thumfart.Mansuy.2023) XBmRNA46727
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA46727", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# GATA-binding protein 4 (Gata4.L)  (also Lydig and myoid cells but not germ cells; Thumfart.Mansuy.2023) XBmRNA45909
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA45909", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# GATA-binding protein 4 (Gata4.S)  (also Lydig and myoid cells but not germ cells; Thumfart.Mansuy.2023) XBmRNA49114
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA49114", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
# germ cells
#ddx25.L XBmRNA58557
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA58557", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
#dnd1.L XBmRNA22633  *** significantly differnetially expressed
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA22633", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
#nanos1.L XBmRNA56602
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA56602", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
#spire1.L XBmRNA52120
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA52120", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]
#spire1.S XBmRNA79956
wtkomalesonly_dmrt1L_unfiltered[grep("XBmRNA79956", rownames(wtkomalesonly_dmrt1L_unfiltered)), ]

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./wtkomalesonly_dmrt1L_rsquare.xls")


# MF dmrt1L ko_ko ----
# dmrt1L females ko vs wt ----
dds <- DESeqDataSetFromMatrix(countData = round(dmrt1L_counts_ko_MFonly),
                              colData = genotypez_ko_MFonly,
                              design= ~ sex)

# save the unfiltered log2FoldChange to a dataframe 
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
kokoMFonly_dmrt1L_unfiltered <- res;kokoMFonly_dmrt1L_unfiltered
# Only sex related
sex_related_kokoMFonly_dmrt1L <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_kokoMFonly_dmrt1L, file="Sex_related_kokoMFonly_dmrt1L_Kallisto_DeSeq2_unfiltered.csv", row.names = T)


# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
# here we remove transcripts where the average count per sample is 1or more
#keep <- rowSums(counts(dds)) >=length(colnames(dds))
dds <- dds[keep,]
dim(dds)
# [1] 32040    12
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:10153,];p
write.csv(p, file="kokoMFonly_dmrt1L_DE_DeSeq2.csv", row.names = T)

# (from gff3 file)
# stra8.L: XBXL10_1g15025 XBmRNA27986
kokoMFonly_dmrt1L_unfiltered[grep("XBmRNA27986", rownames(kokoMFonly_dmrt1L_unfiltered)), ]
res[grep("XBmRNA27986", rownames(res)), ]
kokoMFonly_dmrt1L_unfiltered[grep("XBmRNA27986", rownames(kokoMFonly_dmrt1L_unfiltered)), ]

# How many sign DE genes are sex related?
significantDE_that_are_sexrelated_kokoMFonly_dmrt1L <- resOrdered[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(resOrdered)), ]
significantDE_that_are_sexrelated_kokoMFonly_dmrt1L_df <- as.data.frame(significantDE_that_are_sexrelated_kokoMFonly_dmrt1L)
significantDE_that_are_sexrelated_kokoMFonly_dmrt1L_df <- significantDE_that_are_sexrelated_kokoMFonly_dmrt1L_df[!(is.na(significantDE_that_are_sexrelated_kokoMFonly_dmrt1L_df$padj)), ]
significantDE_that_are_sexrelated_kokoMFonly_dmrt1L_df[(significantDE_that_are_sexrelated_kokoMFonly_dmrt1L_df$padj < 0.05),]
#                                 baseMean log2FoldChange     lfcSE     stat       pvalue       padj
# gnl|XBXL10_1g7278|XBmRNA13700|  442.201      0.9200507 0.2465066 3.732358 0.0001896957 0.01828354
# this is rspo1.L; this suggests that absence of dmrt1L may lead to sex differences in rspo1.L expression
# because it is more highly expressed in female ko compared to female wt...

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./koko_MFonly_dmrt1L_rsquare.xls")





##################
###
###   DMRT1S
###
##################

# wtko dmrt1S ----
# we have two batches for dmrt1S that needs to be accounted for.
# load count data (from Kalisto)
# read the count data from Kallisto that was combined from each sample
# into a single file into a dataframe called "dds"
#dmrt1S_counts <- read.table("dmrt1S.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
dmrt1S_counts <- read.table("dmrt1S_NEWall.isoform.counts.matrix", header=T, row.names = 1)
dim(dmrt1S_counts)
# 44441    22
#coldata <- read.table("dmrt1S_sample_genotype.txt", header=T, row.names = 1)
coldata <- read.table("dmrt1S_sample_genotype_NEWall.txt", header=T, row.names = 1)
coldata$batch <- as.factor(coldata$batch)
coldata$genotype <- as.factor(coldata$genotype)
coldata$genotype <- relevel(coldata$genotype, ref="wt")
colnames(dmrt1S_counts) <- rownames(coldata)

# subset the males and females into separate dataframes
dmrt1S_counts_femalesonly <- dmrt1S_counts[,c(1,3,6:10,12:14,17,20,22)];colnames(dmrt1S_counts_femalesonly)
dmrt1S_counts_malesonly <- dmrt1S_counts[,c(2,4:5,11,15:16,18:19,21)];colnames(dmrt1S_counts_malesonly)
# subset the genotypes
genotypez_femsonly <- coldata[c(1,3,6:10,12:14,17,20,22),];genotypez_femsonly
genotypez_femsonly$genotype <- relevel(genotypez_femsonly$genotype, ref="wt")
genotypez_malesonly <- coldata[c(2,4:5,11,15:16,18:19,21),];genotypez_malesonly
genotypez_malesonly$genotype <- relevel(genotypez_malesonly$genotype, ref="wt")

# dmrt1S females ko vs wt ----
# batch effects dealt with as suggested here: https://rdrr.io/bioc/DESeq2/f/vignettes/DESeq2.Rmd
# and here https://support.bioconductor.org/p/121408/
dds <- DESeqDataSetFromMatrix(countData = round(dmrt1S_counts_femalesonly),
                              colData = genotypez_femsonly,
                              design= ~ batch + genotype)

# save the unfiltered log2FoldChange to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
wtkofemsonly_dmrt1S_unfiltered <- res;wtkofemsonly_dmrt1S_unfiltered
# Only sex related
sex_related_wtkofemsonly_dmrt1S <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_wtkofemsonly_dmrt1S, file="Sex_related_wtkofemsonly_dmrt1S_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# get normalized counts as detailed here: https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
normalized_countz <- counts(dds, normalized=T)
# Get sex related countz
sex_related_wtkofemsonly_dmrt1S_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(normalized_countz)), ]
write.csv(sex_related_wtkofemsonly_dmrt1S_normalized_countz, file="Sex_related_wtkofemsonly_dmrt1S_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)

# Get germcell related countz
#germcell_related_wtkofemsonly_dmrt1L_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g35317\\||gnl\\|XBXL10_1g38173\\|", rownames(normalized_countz)), ]
#write.csv(germcell_related_wtkofemsonly_dmrt1L_normalized_countz, file="Germcell_related_wtkofemsonly_dmrt1L_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)


# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
# here we remove transcripts where the average count per sample is 1or more
#keep <- rowSums(counts(dds)) >=length(colnames(dds))
dds <- dds[keep,]
dim(dds)
# [1] 32098    13
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:13,];p
write.csv(p, file="wtkofemsonly_dmrt1S_DE_DeSeq2_NEWall.csv", row.names = T)

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./wtkofemsonly_dmrt1S_rsquare_NEWall.xls")




# dmrt1S males ko vs wt ----
dds <- DESeqDataSetFromMatrix(countData = round(dmrt1S_counts_malesonly),
                              colData = genotypez_malesonly,
                              design= ~ batch + genotype)

# save the unfiltered log2FoldChange to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
wtkomalesonly_dmrt1S_unfiltered <- res;wtkomalesonly_dmrt1S_unfiltered
# Only sex related
sex_related_wtkomalesonly_dmrt1S <- res[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(res)), ]
write.csv(sex_related_wtkomalesonly_dmrt1S, file="sex_related_wtkomalesonly_dmrt1S_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# get normalized counts as detailed here: https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
normalized_countz <- counts(dds, normalized=T)
# Get sex related countz
sex_related_wtkomalzonly_dmrt1S_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(normalized_countz)), ]
write.csv(sex_related_wtkomalzonly_dmrt1S_normalized_countz, file="Sex_related_wtkomalzonly_dmrt1S_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)

# Get germcell related countz
#germcell_related_wtkomalzonly_dmrt1L_normalized_countz <- normalized_countz[grepl("gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g35317\\||gnl\\|XBXL10_1g38173\\|", rownames(normalized_countz)), ]
#write.csv(germcell_related_wtkomalzonly_dmrt1L_normalized_countz, file="Germcell_related_wtkomalzonly_dmrt1L_Kallisto_DeSeq2_unfiltered_normalized_countz.csv", row.names = T)


# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
# here we remove transcripts where the average count per sample is 1or more
#keep <- rowSums(counts(dds)) >=length(colnames(dds))
dds <- dds[keep,]
dim(dds)
# [1] 32262     9
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:1373,];p
write.csv(p, file="wtkomalesonly_dmrt1S_DE_DeSeq2_NEWall.csv", row.names = T)

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./wtkomalesonly_dmrt1S_rsquare.xls")


# permutations ----

# OK now do some permutations to assess significance of the difference between correlations
# between the log2FC of the ko:ko and each MF

# These permutations will randomly select 90 log2FoldChange from each MF and wtko
# 1000 times and calculate the correlation.  Then this will be compared to
# the observed

#MFcounts <- read.table("MF_NEW.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
MFcounts <- read.table("MF_NEW.isoform.counts.matrix", header=T, row.names = 1)


# get rownames of sexrelated transcripts
SL_rownames <- row.names(MFcounts[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(MFcounts)), ])

# the two functions below were developed by Ian Dworkin and colleagues
# https://github.com/DworkinLab/Trypoxylus_RNAseq/blob/master/analysis_scripts/re_analysis_scripts_jan_2018.Rmd

# This function calculates the Euclidean Distances (the L2 norm) 
# so we can use unit vectors in the estimation of the vector 
# correlation.
# These functions follow the logic of Kuruvilla et al 2002, 
# and were adapted from Pitchers et al 2013
PD <- function(x) { 
  sqrt(t(x)%*%x)}

#this function gives the vector correlation and angle, 
# and vector magnitude ratio, alpha, between two vectors
# alpha of 1 means that the length of vector 2 is the same 
# as the length as vector 1. 
# lower than one means that it is smaller, greater than 
# 1 means that vector 2 is larger
ang.vec.alph <- function(vec1, vec2) {
  vec1 <- vec1 - mean(vec1)
  vec2 <- vec2 - mean(vec2)
  vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
  vec.angle <- acos(vec.cor)*(180/pi)
  vec.alpha <- PD(vec1)/PD(vec2) 
  vec.ED <- PD(vec2-vec1) #Subtract vector one from vector two and then calculate PD for Euclidean Distance. 
  return(c(vector.cor=vec.cor, vec.angle=vec.angle, vec.alpha=vec.alpha, vector.ED=vec.ED))} 










# MF_ccdc vs dmrt1Lmales ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkomalesonly_dmrt1L_trim <- wtkomalesonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkomalesonly_dmrt1L_trim<- wtkomalesonly_dmrt1L_trim[-which(wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtkomalesonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkomalesonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L_trim[-which(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.2341081
print("pvalue: "); rank(correlations)[1001]/1001 # for males just use rank
# because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.8231768

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.07192807




# MF_dmrt1L vs dmrt1Lfems ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkofemsonly_dmrt1L_trim <- wtkofemsonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkofemsonly_dmrt1L_trim<- wtkofemsonly_dmrt1L_trim[-which(wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtkofemsonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkofemsonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L_trim[-which(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.1156141
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.978022

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.5074925





###############
##
##  BElow not used
##
###############


# OK now do some permutations to assess significance of the correlations
# between the log2FC of each wtko and each MF

# These permutations will randomly select 74 log2FoldChange from each MF and wtko
# 1000 times and calculate the correlation.  Then this will be compared to
# the observed

#MFcounts <- read.table("MF_NEW.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
MFcounts <- read.table("MF_NEW.isoform.counts.matrix", header=T, row.names = 1)


# get rownames of sexrelated transcripts
SL_rownames <- row.names(MFcounts[grepl("gnl\\|XBXL10_1g8729\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g34871\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g38948\\||gnl\\|XBXL10_1g42200\\||gnl\\|XBXL10_1g26898\\||gnl\\|XBXL10_1g28819\\||gnl\\|XBXL10_1g34124\\||gnl\\|XBXL10_1g38414\\||gnl\\|XBXL10_1g3180\\||gnl\\|XBXL10_1g5695\\||gnl\\|XBXL10_1g11466\\||gnl\\|XBXL10_1g458\\||gnl\\|XBXL10_1g3570\\||gnl\\|XBXL10_1g30978\\||gnl\\|XBXL10_1g12137\\||gnl\\|XBXL10_1g30063\\||gnl\\|XBXL10_1g42533\\|", rownames(MFcounts)), ])

# the two functions below were developed by Ian Dworkin and colleagues
# https://github.com/DworkinLab/Trypoxylus_RNAseq/blob/master/analysis_scripts/re_analysis_scripts_jan_2018.Rmd

# This function calculates the Euclidean Distances (the L2 norm) 
# so we can use unit vectors in the estimation of the vector 
# correlation.
# These functions follow the logic of Kuruvilla et al 2002, 
# and were adapted from Pitchers et al 2013
PD <- function(x) { 
  sqrt(t(x)%*%x)}

#this function gives the vector correlation and angle, 
# and vector magnitude ratio, alpha, between two vectors
# alpha of 1 means that the length of vector 2 is the same 
# as the length as vector 1. 
# lower than one means that it is smaller, greater than 
# 1 means that vector 2 is larger
ang.vec.alph <- function(vec1, vec2) {
  vec1 <- vec1 - mean(vec1)
  vec2 <- vec2 - mean(vec2)
  vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
  vec.angle <- acos(vec.cor)*(180/pi)
  vec.alpha <- PD(vec1)/PD(vec2) 
  vec.ED <- PD(vec2-vec1) #Subtract vector one from vector two and then calculate PD for Euclidean Distance. 
  return(c(vector.cor=vec.cor, vec.angle=vec.angle, vec.alpha=vec.alpha, vector.ED=vec.ED))} 




# MF_ccdc vs dmrt1Lfems ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkofemsonly_dmrt1L_trim <- wtkofemsonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkofemsonly_dmrt1L_trim<- wtkofemsonly_dmrt1L_trim[-which(wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtkofemsonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkofemsonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L_trim[-which(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.01327707
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.6353646

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.3506494



# MF_ccdc vs dmrt1Lmales ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkomalesonly_dmrt1L_trim <- wtkomalesonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkomalesonly_dmrt1L_trim<- wtkomalesonly_dmrt1L_trim[-which(wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtkomalesonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkomalesonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L_trim[-which(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.2341081
print("pvalue: "); rank(correlations)[1001]/1001 # for males just use rank
                                                # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.8311688

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.06393606



# MF_dmrt1L vs dmrt1Lfems ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkofemsonly_dmrt1L_trim <- wtkofemsonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkofemsonly_dmrt1L_trim<- wtkofemsonly_dmrt1L_trim[-which(wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtkofemsonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkofemsonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L_trim[-which(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.1156141
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.971029

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.5214785



# MF_dmrt1L vs dmrt1Lmales ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkomalesonly_dmrt1L_trim <- wtkomalesonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkomalesonly_dmrt1L_trim<- wtkomalesonly_dmrt1L_trim[-which(wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtkomalesonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkomalesonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L_trim[-which(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.367098
print("pvalue: "); rank(correlations)[1001]/1001 # for males just use rank
                                                  # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.3266733

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1778222



# MF_dmrt1S vs dmrt1Lfems ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkofemsonly_dmrt1L_trim <- wtkofemsonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkofemsonly_dmrt1L_trim<- wtkofemsonly_dmrt1L_trim[-which(wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'log2FoldChange'],
                         wtkofemsonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkofemsonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1S_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkofemsonly_dmrt1L_trim <- sex_related_wtkofemsonly_dmrt1L_trim[-which(sex_related_wtkofemsonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.005925032
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2907093

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkofemsonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.8131868



# MF_dmrt1S vs dmrt1Lmales ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkomalesonly_dmrt1L_trim <- wtkomalesonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkomalesonly_dmrt1L_trim<- wtkomalesonly_dmrt1L_trim[-which(wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'log2FoldChange'],
                         wtkomalesonly_dmrt1L_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkomalesonly_dmrt1L_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1S_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkomalesonly_dmrt1L_trim <- sex_related_wtkomalesonly_dmrt1L_trim[-which(sex_related_wtkomalesonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
#  -0.05217883
print("pvalue: "); rank(correlations)[1001]/1001  # for males just use rank
                                                  # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.4675325

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkomalesonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.5164835


# MF_ccdc vs dmrt1Sfems ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkofemsonly_dmrt1S_trim <- wtkofemsonly_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(wtkofemsonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkofemsonly_dmrt1S_trim<- wtkofemsonly_dmrt1S_trim[-which(wtkofemsonly_dmrt1S_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtkofemsonly_dmrt1S_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkofemsonly_dmrt1S_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkofemsonly_dmrt1S_trim <- sex_related_wtkofemsonly_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_wtkofemsonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkofemsonly_dmrt1S_trim <- sex_related_wtkofemsonly_dmrt1S_trim[-which(sex_related_wtkofemsonly_dmrt1S_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkofemsonly_dmrt1S_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.1086432
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2297702

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkofemsonly_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.02697303



# MF_ccdc vs dmrt1Smales ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkomalesonly_dmrt1S_trim <- wtkomalesonly_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(wtkomalesonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkomalesonly_dmrt1S_trim<- wtkomalesonly_dmrt1S_trim[-which(wtkomalesonly_dmrt1S_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtkomalesonly_dmrt1S_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkomalesonly_dmrt1S_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkomalesonly_dmrt1S_trim <- sex_related_wtkomalesonly_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_wtkomalesonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkomalesonly_dmrt1S_trim <- sex_related_wtkomalesonly_dmrt1S_trim[-which(sex_related_wtkomalesonly_dmrt1S_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkomalesonly_dmrt1S_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.1426504
print("pvalue: "); rank(correlations)[1001]/1001 # for males just use rank
                                                  # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.06693307

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkomalesonly_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.001998002




# MF_dmrt1L vs dmrt1Sfems ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkofemsonly_dmrt1S_trim <- wtkofemsonly_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(wtkofemsonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkofemsonly_dmrt1S_trim<- wtkofemsonly_dmrt1S_trim[-which(wtkofemsonly_dmrt1S_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtkofemsonly_dmrt1S_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkofemsonly_dmrt1S_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkofemsonly_dmrt1S_trim <- sex_related_wtkofemsonly_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_wtkofemsonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkofemsonly_dmrt1S_trim <- sex_related_wtkofemsonly_dmrt1S_trim[-which(sex_related_wtkofemsonly_dmrt1S_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkofemsonly_dmrt1S_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.1893069
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2187812

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkofemsonly_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1428571



# MF_dmrt1L vs dmrt1Smales ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkomalesonly_dmrt1S_trim <- wtkomalesonly_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(wtkomalesonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkomalesonly_dmrt1S_trim<- wtkomalesonly_dmrt1S_trim[-which(wtkomalesonly_dmrt1S_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtkomalesonly_dmrt1S_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkomalesonly_dmrt1S_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkomalesonly_dmrt1S_trim <- sex_related_wtkomalesonly_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_wtkomalesonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkomalesonly_dmrt1S_trim <- sex_related_wtkomalesonly_dmrt1S_trim[-which(sex_related_wtkomalesonly_dmrt1S_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkomalesonly_dmrt1S_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.1127484
print("pvalue: "); rank(correlations)[1001]/1001  # for males just use rank
                                                  # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.1118881

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkomalesonly_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.02997003



# MF_dmrt1S vs dmrt1Sfems ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkofemsonly_dmrt1S_trim <- wtkofemsonly_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(wtkofemsonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkofemsonly_dmrt1S_trim<- wtkofemsonly_dmrt1S_trim[-which(wtkofemsonly_dmrt1S_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'log2FoldChange'],
                         wtkofemsonly_dmrt1S_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkofemsonly_dmrt1S_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1S_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkofemsonly_dmrt1S_trim <- sex_related_wtkofemsonly_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_wtkofemsonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkofemsonly_dmrt1S_trim <- sex_related_wtkofemsonly_dmrt1S_trim[-which(sex_related_wtkofemsonly_dmrt1S_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkofemsonly_dmrt1S_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.6477949
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.07692308

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkofemsonly_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.4675325



# MF_dmrt1S vs dmrt1Smales ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtkomalesonly_dmrt1S_trim <- wtkomalesonly_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(wtkomalesonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtkomalesonly_dmrt1S_trim<- wtkomalesonly_dmrt1S_trim[-which(wtkomalesonly_dmrt1S_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'log2FoldChange'],
                         wtkomalesonly_dmrt1S_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtkomalesonly_dmrt1S_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1S_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtkomalesonly_dmrt1S_trim <- sex_related_wtkomalesonly_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_wtkomalesonly_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtkomalesonly_dmrt1S_trim <- sex_related_wtkomalesonly_dmrt1S_trim[-which(sex_related_wtkomalesonly_dmrt1S_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtkomalesonly_dmrt1S_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.5453617
print("pvalue: "); rank(correlations)[1001]/1001 # for males just use rank
                                                # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.1568432


# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtkomalesonly_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.3616384


# below not used ----
# MF_ccdc vs koko_MF (not used) ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 90, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from koko_MF
  kokoMFonly_dmrt1L_unfiltered_trim <- kokoMFonly_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(kokoMFonly_dmrt1L_unfiltered_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    kokoMFonly_dmrt1L_unfiltered_trim<- kokoMFonly_dmrt1L_unfiltered_trim[-which(kokoMFonly_dmrt1L_unfiltered_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         kokoMFonly_dmrt1L_unfiltered_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(kokoMFonly_dmrt1L_unfiltered_trim[,'log2FoldChange'], # ko:ko M:F first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:ko has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from koko
sex_related_kokoMFonly_dmrt1L_trim <- kokoMFonly_dmrt1L_unfiltered[SL_rownames,]
outliers <- boxplot(sex_related_kokoMFonly_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_kokoMFonly_dmrt1L_trim <- sex_related_kokoMFonly_dmrt1L_trim[-which(sex_related_kokoMFonly_dmrt1L_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_kokoMFonly_dmrt1L_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.3226452
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.08791209

# not significant

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_kokoMFonly_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]
# [1] 1.506053
print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.002997003

# Significant; because the ratio is >1, this means that the magnitude of the ko_ko_MF vector is larger than 
# the magnitude of the wt_wt_MF vector. This is the opposite of what I was expecting.
```
