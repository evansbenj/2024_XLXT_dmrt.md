# STAR, EdgeR
```
library(edgeR)
library(tximport)
library('edgeR')
library('rhdf5')
library('readxl')
library('ggplot2')
library(grid)
require('gridExtra')
library("org.Xl.eg.db")
library(PCAtools)
library("HTSFilter")
library(tidyverse)
library(purrr)
library(writexl)


setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_dmrt1_knockouts/2023_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_EdgeR_Analysis1")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_dmrt1_knockouts/2023_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_EdgeR_Analysis1"
list.files(dir)


f_files<- list.files(".", pattern = "_counts.out", full.names = T);f_files

read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names =T, comment = "#")
  cnt<- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)
  return(cnt)
}

raw_counts<- map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join) 

dim(raw_counts_df)
# [1] 44478    78
# View(raw_counts_df)

# get rid of any rows that have incomplete data
counts <- raw_counts_df[complete.cases(raw_counts_df), ]
dim(counts)
# [1]  44478    78 # 2023_STAR

colnames(counts) <- c("geneID","ccdc_11_wt_M","ccdc_12_wt_M","ccdc_13_wt_M",
                      "ccdc_14_ko_F","ccdc_2_wt_M","ccdc_25_wt_M","ccdc_3_wt_F",
                      "ccdc_30_ko_F","ccdc_32_ko_F","ccdc_34_wt_F","ccdc_35_ko_F",
                      "ccdc_36_ko_F","ccdc_42_ko_F","ccdc_9_wt_M","dmrt1L_11_ko_F",
                      "dmrt1L_13_wt_M","dmrt1L_15_wt_M","dmrt1L_17_ko_M",
                      "dmrt1L_19_ko_M","dmrt1L_24_ko_F","dmrt1L_25_ko_F",
                      "dmrt1L_26_ko_M","dmrt1L_27_ko_F","dmrt1L_30_wt_F",
                      "dmrt1L_34_wt_M","dmrt1L_35_ko_M","dmrt1L_41_ko_M",
                      "dmrt1L_43_wt_M","dmrt1L_50_wt_F","dmrt1L_55_wt_F",
                      "dmrt1L_59_ko_M","dmrt1L_6_ko_F","dmrt1L_7_ko_F",
                      "dmrt1L_8_wt_M","dmrt1S_1_wt_F_b2",
                      "dmrt1S_10_ko_F","dmrt1S_11_wt_F_b2","dmrt1S_13_ko_F_b2",
                      "dmrt1S_14_wt_M",
                      "dmrt1S_15_ko_F","dmrt1S_15_ko_F_b2","dmrt1S_16_ko_M_b2",
                      "dmrt1S_18_wt_M","dmrt1S_19_wt_M_b2",
                      "dmrt1S_20_wt_M",
                      "dmrt1S_23_ko_F","dmrt1S_24_wt_F","dmrt1S_25_ko_M_b2",
                      "dmrt1S_29_wt_F",
                      "dmrt1S_30_ko_F","dmrt1S_30_wt_M_b2",
                      "dmrt1S_4_wt_F","dmrt1S_5_ko_M","dmrt1S_6_wt_F_b2",
                      "dmrt1S_8_wt_M_b2","dmrt1S_9_ko_F_b2",
                      "dmw_12_wt_F","dmw_14_ko_F","dmw_15_wt_F","dmw_16_ko_F",
                      "dmw_17_wt_F","dmw_20_wt_F","dmw_22_wt_F","dmw_26_ko_F",
                      "dmw_28_ko_F","dmw_29_ko_F","dmw_35_ko_F","dmw_9_wt_F",
                      "scanw_14_wt_F","scanw_15_wt_F","scanw_19_ko_F",
                      "scanw_27_ko_F","scanw_28_ko_F","scanw_30_ko_F",
                      "scanw_31_ko_F","scanw_6_wt_F","scanw_9_wt_F")

# this is useful info about batch effects:
# https://support.bioconductor.org/p/96627/

# make design matrix
# Here we want to test for differential expression between KO and
# WT, while adjusting for differences between batches. In statistical
# terms, this "may be" an additive linear model with batch as the blocking factor:

# countz
countz <- counts[,-c(1)]
gene_names <- counts$geneID
row.names(countz) <- gene_names
# samples
samples <- read.table(file.path(dir, "samples_with_all_dmrt1S.txt"), header = T)
samples
# sexez
sexez <- factor(samples$sex)
sexez <- relevel(sexez, ref="F")
# batch
batchez <- factor(samples$batch) # for all genes except dmrt1S this refers to a single clutch 
                                # for dmrt1S it is also a single clutch but we can use the design
                                # to also control for lane effects because the samples for each
                                # sex and genotype are balanced on two lanes
                                # for dmw, dmrt1L, ccdc, the samples are not balanced on the lanes
                                # ie one treatment is missing from a lane, so we can't control
                                # for lane effects. for dmrt1L, 3 wt M were run on a different lane
                                # from the rest: dmrt1L_13_wt_M, dmrt1L_15_wt_M,dmrt1L_34_wt_M

# MF dmrt1S ----
# compare wt F to wt M ----
colnames(countz)
# save only wt F and wt M
new_counts <- as.data.frame(countz[,c(47,49,52,35,37,54,39,43,45,44,51,55) ])
row.names(new_counts) <- gene_names
# new_counts['XBXL10_1g4848',] # dmrt1.S
# new_counts['XBXL10_1g2070',] # dmrt1.L
new_samples <-as.data.frame(samples[c(47,49,52,35,37,54,39,43,45,44,51,55), ]);new_samples
new_sexez <- factor(samples$sex[c(47,49,52,35,37,54,39,43,45,44,51,55)])
new_sexez <- relevel(new_sexez, ref="F")
# batch
new_batchez <- factor(samples$batch[c(47,49,52,35,37,54,39,43,45,44,51,55)])
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 36004    12 # 2023 STAR dmrt1S

# d0$counts['XBXL10_1g4848',]
# d0$counts['XBXL10_1g2070',]
# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_batchez + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_dmrt1S_unfiltered <- exacttest$table;MF_dmrt1S_unfiltered
# Write sex_related to a file
sex_related_MF_dmrt1S <- data.frame(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
write.csv(sex_related_MF_dmrt1S, file="Sex_related_MF_dmrt1S_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MF_dmrt1S_counts <- new_counts[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(sex_related_MF_dmrt1S_counts, file="Sex_related_MF_dmrt1S_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# germcell-related
# Sox9.L: XBXL10_1g39526
# sox9.S: XBXL10_1g42722
# Cytokeratin.L: XBXL10_1g38948
# Cytokeratin.S: XBXL10_1g42200
# vimentin.L: XBXL10_1g26898
# vimentin.S: XBXL10_1g28819
# Gata1.L: XBXL10_1g34124
# Gata1.S: XBXL10_1g38414
# FSHr.L: XBXL10_1g22534
# FSHr.S: XBXL10_1g25047
# Gata4.L: XBXL10_1g24554
# Gata4.S: XBXL10_1g26280
# ddx4.L: XBXL10_1g3180
# ddx4.S: XBXL10_1g5695
# cyp17a1.L:XBXL10_1g30377
# foxl2.L: XBXL10_1g24241
# foxl2.S: XBXL10_1g26060
# vgll1.L: XBXL10_1g35317
# vgll1.S: XBXL10_1g38173
# ddx25.L: XBXL10_1g30978
# dnd1.L: XBXL10_1g12137
# nanos1.L: XBXL10_1g30063
# spire1.S: XBXL10_1g42533

# Write counts of germcell related to a file
germ_related_MF_dmrt1S_counts <- new_counts[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(germ_related_MF_dmrt1S_counts, file="Germ_related_MF_dmrt1S_STAR_edgeR_unfiltered.csv", row.names = T)


# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 28900    12 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",3),rep("F",3)),
#        col=c(rep("green",3),rep("blue",3)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_batchez + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") 
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down      77
# NotSig 28795
# Up        28
dmrt1S_DE <- as.data.frame(topTags(exacttest, n=105))
write.csv(dmrt1S_DE, file="MF_STAR_dmrt1S_DE_edgeR.csv", row.names = T)

# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
write_xlsx(rsquare, "./MF_rsquare_dmrt1S_STAR_edgeR.xls")

# negative logFC indicates higher expression in wt because wt is the
# denominator (reference)

#new_counts[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#exacttest$table[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#female related
# View(exacttest$table[c('XBXL10_1g34625','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g29076','XBXL10_1g30057','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g33473','XBXL10_1g35876','XBXL10_1g37293','XBXL10_1g5748','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8966'),])
# male related
# View(exacttest$table[c('XBXL10_1g2070','XBXL10_1g4848','XBXL10_1g11002','XBXL10_1g30252','XBXL10_1g32546','XBXL10_1g605','XBXL10_1g3639','XBXL10_1g37486'),])
# sox genes
# View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g35158','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g42158','XBXL10_1g39443','XBXL10_1g42662'),])
# testis differentiation
# View(exacttest$table[c('XBXL10_1g41173','XBXL10_1g43880','XBXL10_1g19698','XBXL10_1g22028','XBXL10_1g815','XBXL10_1g3800','XBXL10_1g8007','XBXL10_1g2154','XBXL10_1g4928','XBXL10_1g27310','XBXL10_1g29128','XBXL10_1g40425','XBXL10_1g43291','XBXL10_1g8118','XBXL10_1g10760','XBXL10_1g8117','XBXL10_1g10758','XBXL10_1g1634','XBXL10_1g4460'),])
# steroidogenic genes
# View(exacttest$table[c('XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g30377','XBXL10_1g13205','XBXL10_1g6054','XBXL10_1g29226','XBXL10_1g34871'),])
# Germ cells
View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])

# compare wt F to ko F dmrt1S ----
colnames(countz)
# save only wt F and wt M
new_counts <- as.data.frame(countz[,c(36,40,46,50,38,41,56,47,49,52,35,37,54) ])
row.names(new_counts) <- gene_names
# new_counts['XBXL10_1g4848',]
# new_counts['XBXL10_1g2070',]
new_samples <-as.data.frame(samples[c(36,40,46,50,38,41,56,47,49,52,35,37,54), ]);new_samples
new_genotypez <- factor(samples$genotype[c(36,40,46,50,38,41,56,47,49,52,35,37,54)])
new_genotypez <- relevel(new_genotypez, ref="wt")
# batch
new_batchez <- factor(samples$batch[c(36,40,46,50,38,41,56,47,49,52,35,37,54)])
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_genotypez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 36265    13 # 2023 STAR dmrt1S F to F

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_batchez + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
FF_dmrt1S_unfiltered <- exacttest$table;FF_dmrt1S_unfiltered
# Write sex_related to a file
sex_related_FF_dmrt1S <- data.frame(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
write.csv(sex_related_FF_dmrt1S, file="Sex_related_FF_dmrt1S_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_FF_dmrt1S_counts <- new_counts[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(sex_related_FF_dmrt1S_counts, file="Sex_related_FF_dmrt1S_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# d0$counts['XBXL10_1g29991',]
# d0$counts['XBXL10_1g2070',]
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29119    13 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
plotMDS(d0,labels=c(rep("ko",4),rep("ko2",3),rep("wt",3),rep("wt2",3)),
        col=c(rep("green",4),rep("lightgreen",3),rep("blue",3),rep("steelblue",3)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_batchez + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") 
summary(decideTests(object = exacttest, p.value = 0.1))
# ko-wt
# Down       9
# NotSig 29099
# Up        11
dmrt1S_FF_DE <- as.data.frame(topTags(exacttest, n=20))
write.csv(dmrt1S_FF_DE, file="FF_STAR_dmrt1S_DE_edgeR.csv", row.names = T)

# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
write_xlsx(rsquare, "./FF_rsquare_dmrt1S_STAR_edgeR.xls")

# negative logFC indicates higher expression in wt because wt is the
# denominator (reference)

#new_counts[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#exacttest$table[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#female related
# View(exacttest$table[c('XBXL10_1g34625','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g29076','XBXL10_1g30057','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g33473','XBXL10_1g35876','XBXL10_1g37293','XBXL10_1g5748','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8966'),])
# male related
# View(exacttest$table[c('XBXL10_1g2070','XBXL10_1g4848','XBXL10_1g11002','XBXL10_1g30252','XBXL10_1g32546','XBXL10_1g605','XBXL10_1g3639','XBXL10_1g37486'),])
# sox genes
# View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g35158','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g42158','XBXL10_1g39443','XBXL10_1g42662'),])
# testis differentiation
# View(exacttest$table[c('XBXL10_1g41173','XBXL10_1g43880','XBXL10_1g19698','XBXL10_1g22028','XBXL10_1g815','XBXL10_1g3800','XBXL10_1g8007','XBXL10_1g2154','XBXL10_1g4928','XBXL10_1g27310','XBXL10_1g29128','XBXL10_1g40425','XBXL10_1g43291','XBXL10_1g8118','XBXL10_1g10760','XBXL10_1g8117','XBXL10_1g10758','XBXL10_1g1634','XBXL10_1g4460'),])
# steroidogenic genes
# View(exacttest$table[c('XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g30377','XBXL10_1g13205','XBXL10_1g6054','XBXL10_1g29226','XBXL10_1g34871'),])
# check out differential expression of sex-related genes
View(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
# check out differential expression of germcell related genes
View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])



# compare wt M to ko M dmrt1S ----
colnames(countz)
# save only wt F and wt M
new_counts <- as.data.frame(countz[,c(53,42,48,39,43,45,44,51,55) ])
row.names(new_counts) <- gene_names
# new_counts['XBXL10_1g4848',]
# new_counts['XBXL10_1g2070',]
new_samples <-as.data.frame(samples[c(53,42,48,39,43,45,44,51,55), ]);new_samples
new_genotypez <- factor(samples$genotype[c(53,42,48,39,43,45,44,51,55)])
new_genotypez <- relevel(new_genotypez, ref="wt")
# batch
new_batchez <- factor(samples$batch[c(53,42,48,39,43,45,44,51,55)])
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_genotypez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 35668     9 # 2023 STAR dmrt1S F to F

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_batchez + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MM_dmrt1S_unfiltered <- exacttest$table;MM_dmrt1S_unfiltered
# Write sex_related to a file
sex_related_MM_dmrt1S <- data.frame(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
write.csv(sex_related_MM_dmrt1S, file="Sex_related_MM_dmrt1S_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MM_dmrt1S_counts <- new_counts[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(sex_related_MM_dmrt1S_counts, file="Sex_related_MM_dmrt1S_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29232     9 # 2023 STAR ccdc dmrt1L no dmrt1S
# check one
# d0$counts['XBXL10_1g9904',]
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
plotMDS(d0,labels=c(rep("ko",1),rep("ko2",2),rep("wt",3),rep("wt2",3)),
        col=c(rep("green",1),rep("lightgreen",2),rep("blue",3),rep("steelblue",3)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_batchez + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") 
summary(decideTests(object = exacttest, p.value = 0.1))
# ko-wt
# Down     654
# NotSig 27981
# Up       597
dmrt1S_MM_DE <- as.data.frame(topTags(exacttest, n=1251))
write.csv(dmrt1S_MM_DE, file="MM_STAR_dmrt1S_DE_edgeR.csv", row.names = T)

# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )


# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
write_xlsx(rsquare, "./MM_rsquare_dmrt1S_STAR_edgeR.xls")

# negative logFC indicates higher expression in wt because wt is the
# denominator (reference)

#new_counts[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#exacttest$table[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#female related
# View(exacttest$table[c('XBXL10_1g34625','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g29076','XBXL10_1g30057','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g33473','XBXL10_1g35876','XBXL10_1g37293','XBXL10_1g5748','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8966'),])
# male related
# View(exacttest$table[c('XBXL10_1g2070','XBXL10_1g4848','XBXL10_1g11002','XBXL10_1g30252','XBXL10_1g32546','XBXL10_1g605','XBXL10_1g3639','XBXL10_1g37486'),])
# sox genes
# View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g35158','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g42158','XBXL10_1g39443','XBXL10_1g42662'),])
# testis differentiation
# View(exacttest$table[c('XBXL10_1g41173','XBXL10_1g43880','XBXL10_1g19698','XBXL10_1g22028','XBXL10_1g815','XBXL10_1g3800','XBXL10_1g8007','XBXL10_1g2154','XBXL10_1g4928','XBXL10_1g27310','XBXL10_1g29128','XBXL10_1g40425','XBXL10_1g43291','XBXL10_1g8118','XBXL10_1g10760','XBXL10_1g8117','XBXL10_1g10758','XBXL10_1g1634','XBXL10_1g4460'),])
# steroidogenic genes
# View(exacttest$table[c('XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g30377','XBXL10_1g13205','XBXL10_1g6054','XBXL10_1g29226','XBXL10_1g34871'),])
# check out differential expression of sex-related genes
View(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
# check out differential expression of germcell related genes
View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
# check out differential expression of dmrt1.S
View(exacttest$table[c('XBXL10_1g4848'),])

# MF dmrt1L ----
# compare wt F to wt M ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(24,29,30,16,17,25,28,34) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[c(24,29,30,16,17,25,28,34), ]);new_samples
new_sexez <- factor(samples$sex[c(24,29,30,16,17,25,28,34)])
new_sexez <- relevel(new_sexez, ref="F")
# batch
new_batchez <- factor(samples$batch[c(24,29,30,16,17,25,28,34)])
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 35000     8 # 2023 STAR ccdc dmrt1L no dmrt1S

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_dmrt1L_unfiltered <- exacttest$table;MF_dmrt1L_unfiltered
# Write sex_related to a file
sex_related_MF_dmrt1L <- data.frame(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
write.csv(sex_related_MF_dmrt1L, file="Sex_related_MF_dmrt1L_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MF_dmrt1L_counts <- new_counts[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(sex_related_MF_dmrt1L_counts, file="Sex_related_MF_dmrt1L_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 28721     8 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",2),"F",rep("M",2),rep("F",2),"M"),
#        col=c(rep("blue",2),"red",rep("blue",2),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down      26
# NotSig 28689
# Up         6
topTags(exacttest, n=32)

dmrt1L_DE <- as.data.frame(topTags(exacttest, n=32))
write.csv(dmrt1L_DE, file="MF_STAR_dmrt1L_DE_edgeR.csv", row.names = T)


# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
library(writexl)
write_xlsx(rsquare, "./MF_rsquare_dmrt1L_STAR_edgeR.xls")

View(exacttest$table[c('XBXL10_1g8729','XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])

# compare wt F to ko F dmrt1L ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(24,29,30,15,20,21,23,32,33) ])
row.names(new_counts) <- gene_names
# new_counts['XBXL10_1g4848',]
new_samples <-as.data.frame(samples[c(24,29,30,15,20,21,23,32,33), ]);new_samples
new_genotypez <- factor(samples$genotype[c(24,29,30,15,20,21,23,32,33)])
new_genotypez <- relevel(new_genotypez, ref="wt")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_genotypez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 35455     9 # 2023 STAR ccdc dmrt1L no dmrt1S

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
FF_dmrt1L_unfiltered <- exacttest$table;FF_dmrt1L_unfiltered
# Write sex_related to a file
sex_related_FF_dmrt1L <- data.frame(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
write.csv(sex_related_FF_dmrt1L, file="Sex_related_FF_dmrt1L_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_FF_dmrt1L_counts <- new_counts[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(sex_related_FF_dmrt1L_counts, file="Sex_related_FF_dmrt1L_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 28633     9 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",2),"F",rep("M",2),rep("F",2),"M"),
#        col=c(rep("blue",2),"red",rep("blue",2),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# ko-wt
# Down     223
# NotSig 28296
# Up       114
topTags(exacttest, n=337)
dmrt1L_FF_DE <- as.data.frame(topTags(exacttest, n=337))
write.csv(dmrt1L_FF_DE, file="FF_STAR_dmrt1L_DE_edgeR.csv", row.names = T)

# Let's try to get the logFC for all the genes in MF dmrt1L 
# that are sigDE in dmrt1S
is.data.frame(exacttest$table)
# first read in the file from the dmrt1S
dmrt1_DE_FF_list <- read.csv(file.path(dir, "FF_STAR_dmrt1S_DE_edgeR.csv"), header = T)
dmrt1_DE_FF_list_names <- dmrt1_DE_FF_list$X
write.csv(exacttest$table[dmrt1_DE_FF_list_names,], file="FF_STAR_edgeR_dmrt1Lexpression_of_dmrt1L_DEs.csv", row.names = T)


# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
library(writexl)
write_xlsx(rsquare, "./FF_rsquare_dmrt1L_STAR_edgeR.xls")

# check out differential expression of sex-related genes
View(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
# check out differential expression of germcell related genes
View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])


# compare wt M to ko M dmrt1L ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(16,17,25,28,34,18,19,22,26,27,31) ])
row.names(new_counts) <- gene_names
# new_counts['XBXL10_1g4848',] # dmrt1.S
new_samples <-as.data.frame(samples[c(16,17,25,28,34,18,19,22,26,27,31), ]);new_samples
new_genotypez <- factor(samples$genotype[c(16,17,25,28,34,18,19,22,26,27,31)])
new_genotypez <- relevel(new_genotypez, ref="wt")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_genotypez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 35666    11 # 2023 STAR ccdc dmrt1L no dmrt1S

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MM_dmrt1L_unfiltered <- exacttest$table;MM_dmrt1L_unfiltered
# Write sex_related to a file
sex_related_MM_dmrt1L <- data.frame(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
write.csv(sex_related_MM_dmrt1L, file="Sex_related_MM_dmrt1L_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MM_dmrt1L_counts <- new_counts[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(sex_related_MM_dmrt1L_counts, file="Sex_related_MM_dmrt1L_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 28898    11 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",2),"F",rep("M",2),rep("F",2),"M"),
#        col=c(rep("blue",2),"red",rep("blue",2),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# ko-wt
# Down       2
# NotSig 28893
# Up         3
topTags(exacttest, n=5)
dmrt1L_MM_DE <- as.data.frame(topTags(exacttest, n=5))
write.csv(dmrt1L_MM_DE, file="MM_STAR_dmrt1L_DE_edgeR.csv", row.names = T)


# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
library(writexl)
write_xlsx(rsquare, "./MM_rsquare_dmrt1L_STAR_edgeR.xls")

# check out sex related
View(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])

# check out differential expression of germcell related genes
View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])


# MF ccdc ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(7,10,1,2,3,5,6,14) ])
row.names(new_counts) <- gene_names
# new_counts['XBXL10_1g4848',]

new_samples <-as.data.frame(samples[c(7,10,1,2,3,5,6,14), ]);new_samples
new_sexez <- factor(samples$sex[c(7,10,1,2,3,5,6,14)]);new_sexez
new_sexez <- relevel(new_sexez, ref="F")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 34467     8 # 2023 STAR ccdc only

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_ccdc_unfiltered <- exacttest$table;MF_ccdc_unfiltered
# Write sex_related to a file
sex_related_MF_ccdc <- data.frame(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])
write.csv(sex_related_MF_ccdc, file="Sex_related_MF_ccdc_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MF_ccdc_counts <- new_counts[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'), ]
write.csv(sex_related_MF_ccdc_counts, file="Sex_related_MF_dmrt1L_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29073     8 # 2023 STAR ccdc only
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",5),rep("F",2),"M"),
#        col=c(rep("blue",5),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down       2
# NotSig 29070
#  Up         1
topTags(exacttest, n=3)
# Coefficient:  new_sexezm 
# logFC   logCPM        F       PValue        FDR
# XBXL10_1g3050  -7.520457 -0.6587820 1.307588e-06 0.01585958
# XBXL10_1g43040 -8.439321 -1.6729722 1.477786e-06 0.01585958
# XBXL10_1g3473   9.507057  0.7988971 1.636526e-06 0.01585958
ccdc_DE <- as.data.frame(topTags(exacttest, n=3))
write.csv(ccdc_DE, file="MF_STAR_ccdc_DE_edgeR.csv", row.names = T)


# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
     }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
library(writexl)
write_xlsx(rsquare, "./MF_rsquare_ccdc_STAR_edgeR.xls")
## all comparisons are above 0.8


# check out sex related
View(exacttest$table[c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])

# check out differential expression of germcell related genes
View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g24554','XBXL10_1g26280','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g30377','XBXL10_1g24241','XBXL10_1g26060','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533'),])


# permutations ----

# OK now do some permutations to assess significance of the correlations
# between the log2FC of each wtko and each MF

# These permutations will randomly select 90 logFC from each MF and wtko
# 1000 times and calculate the correlation.  Then this will be compared to
# the observed

# get rownames of sexrelated transcripts
SL_rownames <- c('XBXL10_1g8729','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274','XBXL10_1g38948','XBXL10_1g42200','XBXL10_1g26898','XBXL10_1g28819','XBXL10_1g34124','XBXL10_1g38414','XBXL10_1g3180','XBXL10_1g5695','XBXL10_1g11466','XBXL10_1g458','XBXL10_1g3570','XBXL10_1g30978','XBXL10_1g12137','XBXL10_1g30063','XBXL10_1g42533')
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



# there are 12 permutations in total
# MF_ccdc vs dmrt1L_FF
# MF_ccdc vs dmrt1L_MM
# MF_dmrt1L vs dmrt1L_FF
# MF_dmrt1L vs dmrt1L_MM
# MF_dmrt1S vs dmrt1L_FF
# MF_dmrt1S vs dmrt1L_MM
# MF_ccdc vs dmrt1S_FF
# MF_ccdc vs dmrt1S_MM
# MF_dmrt1L vs dmrt1S_FF
# MF_dmrt1L vs dmrt1S_MM
# MF_dmrt1S vs dmrt1S_FF
# MF_dmrt1S vs dmrt1S_MM

# MF_ccdc vs dmrt1L_FF ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  FF_dmrt1L_trim <- FF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(FF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    FF_dmrt1L_trim<- FF_dmrt1L_trim[-which(FF_dmrt1L_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'logFC'],
                         FF_dmrt1L_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(FF_dmrt1L_trim[,'logFC'], # ko:wt first
             MF_ccdc_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_FF_dmrt1L_trim <- sex_related_FF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_FF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_FF_dmrt1L_trim <- sex_related_FF_dmrt1L_trim[-which(sex_related_FF_dmrt1L_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
                          sex_related_FF_dmrt1L_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# [1] 0.1043459
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.4505495

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_FF_dmrt1L_trim[SL_rownames,'logFC'],
           sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.3376623


# MF_ccdc vs dmrt1L_MM ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  MM_dmrt1L_trim <- MM_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MM_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MM_dmrt1L_trim<- MM_dmrt1L_trim[-which(MM_dmrt1L_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'logFC'],
                         MM_dmrt1L_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(MM_dmrt1L_trim[,'logFC'], # ko:wt first
             MF_ccdc_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_MM_dmrt1L_trim <- sex_related_MM_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MM_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MM_dmrt1L_trim <- sex_related_MM_dmrt1L_trim[-which(sex_related_MM_dmrt1L_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
                          sex_related_MM_dmrt1L_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# [1] 0.1687315
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.3036963

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_MM_dmrt1L_trim[SL_rownames,'logFC'],
           sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1418581



# MF_dmrt1L vs dmrt1LFF ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  FF_dmrt1L_unfiltered_trim <- FF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(FF_dmrt1L_unfiltered_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    FF_dmrt1L_unfiltered_trim<- FF_dmrt1L_unfiltered_trim[-which(FF_dmrt1L_unfiltered_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'logFC'],
                         FF_dmrt1L_unfiltered_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(FF_dmrt1L_unfiltered_trim[,'logFC'], # ko:wt first
             MF_dmrt1L_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_FF_dmrt1L_trim <- sex_related_FF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_FF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_FF_dmrt1L_trim <- sex_related_FF_dmrt1L_trim[-which(sex_related_FF_dmrt1L_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
                          sex_related_FF_dmrt1L_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001] 
# 0.0626931
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1]  0.984016

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_FF_dmrt1L_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.5464535



# MF_dmrt1L vs dmrt1LMM ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  MM_dmrt1L_unfiltered_trim <- MM_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MM_dmrt1L_unfiltered_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MM_dmrt1L_unfiltered_trim<- MM_dmrt1L_unfiltered_trim[-which(MM_dmrt1L_unfiltered_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'logFC'],
                         MM_dmrt1L_unfiltered_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(MM_dmrt1L_unfiltered_trim[,'logFC'], # ko:wt first
             MF_dmrt1L_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_MM_dmrt1L_trim <- sex_related_MM_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MM_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MM_dmrt1L_trim <- sex_related_MM_dmrt1L_trim[-which(sex_related_MM_dmrt1L_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
                          sex_related_MM_dmrt1L_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
#  -0.225467
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.3836164

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_MM_dmrt1L_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.06193806


# MF_dmrt1S vs dmrt1LFF ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  FF_dmrt1L_unfiltered_trim <- FF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(FF_dmrt1L_unfiltered_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    FF_dmrt1L_unfiltered_trim<- FF_dmrt1L_unfiltered_trim[-which(FF_dmrt1L_unfiltered_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'logFC'],
                         FF_dmrt1L_unfiltered_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(FF_dmrt1L_unfiltered_trim[,'logFC'], # ko:wt first
             MF_dmrt1S_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_FF_dmrt1L_trim <-  sex_related_FF_dmrt1L[SL_rownames,]
outliers <- boxplot( sex_related_FF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_FF_dmrt1L_trim <-  sex_related_FF_dmrt1L_trim[-which( sex_related_FF_dmrt1L_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
                          sex_related_FF_dmrt1L_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.01121297
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2607393

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_FF_dmrt1L_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1]  0.8791209


# MF_dmrt1S vs dmrt1LMM ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  MM_dmrt1S_unfiltered_trim <- MM_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MM_dmrt1S_unfiltered_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MM_dmrt1S_unfiltered_trim<- MM_dmrt1S_unfiltered_trim[-which(MM_dmrt1S_unfiltered_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'logFC'],
                         MM_dmrt1S_unfiltered_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(MM_dmrt1S_unfiltered_trim[,'logFC'], # ko:wt first
             MF_dmrt1S_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MM_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S_trim[-which(sex_related_MM_dmrt1S_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
                          sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001] 
#  -0.4297866
print("pvalue: "); rank(correlations)[1001]/1001 # for males just use rank
                                # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.5274725

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# 0.3816184



# MF_ccdc vs dmrt1SFF ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  FF_dmrt1S_trim <- FF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(FF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    FF_dmrt1S_trim<- FF_dmrt1S_trim[-which(FF_dmrt1S_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'logFC'],
                         FF_dmrt1S_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(FF_dmrt1S_trim[,'logFC'], # ko:wt first
             MF_ccdc_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_FF_dmrt1S_trim <- sex_related_FF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_FF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_FF_dmrt1S_trim <- sex_related_FF_dmrt1S_trim[-which(sex_related_FF_dmrt1S_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
                          sex_related_FF_dmrt1S_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.1363979
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# 0.7692308

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_FF_dmrt1S_trim[SL_rownames,'logFC'],
           sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# 0.1868132

# MF_ccdc vs dmrt1SMM ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  MM_dmrt1S_trim <- MM_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MM_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MM_dmrt1S_trim<- MM_dmrt1S_trim[-which(MM_dmrt1S_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'logFC'],
                         MM_dmrt1S_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(MM_dmrt1S_trim[,'logFC'], # ko:wt first
             MF_ccdc_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MM_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S_trim[-which(sex_related_MM_dmrt1S_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
                          sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.01176633
print("pvalue: "); rank(correlations)[1001]/1001  # for males just use rank
                        # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.2787213


# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'],
           sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.04195804





# MF_dmrt1L vs dmrt1S_FF ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  FF_dmrt1S_trim <- FF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(FF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    FF_dmrt1S_trim<- FF_dmrt1S_trim[-which(FF_dmrt1S_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'logFC'],
                         FF_dmrt1S_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(FF_dmrt1S_trim[,'logFC'], # ko:wt first
             MF_dmrt1L_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_FF_dmrt1S_trim <- sex_related_FF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_FF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_FF_dmrt1S_trim <- sex_related_FF_dmrt1S_trim[-which(sex_related_FF_dmrt1S_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
                          sex_related_FF_dmrt1S_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.114198
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.4495504

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_FF_dmrt1S_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1688312



# MF_dmrt1L vs dmrt1S_MM ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  MM_dmrt1S_trim <- MM_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MM_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MM_dmrt1S_trim<- MM_dmrt1S_trim[-which(MM_dmrt1S_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'logFC'],
                         MM_dmrt1S_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(MM_dmrt1S_trim[,'logFC'], # ko:wt first
             MF_dmrt1L_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MM_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S_trim[-which(sex_related_MM_dmrt1S_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
                          sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.068392
print("pvalue: "); rank(correlations)[1001]/1001 # for males just use rank
                          # because we expect a negative correlation
# [1] "pvalue: "
# [1] 0.2197802

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.2227772



# MF_dmrt1S vs dmrt1S_FF ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  FF_dmrt1S_trim <- FF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(FF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    FF_dmrt1S_trim<- FF_dmrt1S_trim[-which(FF_dmrt1S_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'logFC'],
                         FF_dmrt1S_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(FF_dmrt1S_trim[,'logFC'], # ko:wt first
             MF_dmrt1S_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_FF_dmrt1S_trim <- sex_related_FF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_FF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_FF_dmrt1S_trim <- sex_related_FF_dmrt1S_trim[-which(sex_related_FF_dmrt1S_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
                          sex_related_FF_dmrt1S_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# 0.465763
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.4225774

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_FF_dmrt1S_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.5094905



# MF_dmrt1S vs dmrt1S_MM ----
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 90, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  MM_dmrt1S_trim <- MM_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MM_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MM_dmrt1S_trim<- MM_dmrt1S_trim[-which(MM_dmrt1S_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'logFC'],
                         MM_dmrt1S_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(MM_dmrt1S_trim[,'logFC'], # ko:wt first
             MF_dmrt1S_trim[,'logFC'], # reference M:F second
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
outliers <- boxplot(sex_related_MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MM_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MM_dmrt1S_trim <- sex_related_MM_dmrt1S_trim[-which(sex_related_MM_dmrt1S_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
                          sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
correlations[1001]
# -0.4297866
print("pvalue: "); rank(correlations)[1001]/1001  # for males just use rank
                              # because we expect a negative correlation
# [1] "pvalue: "
# [1]  0.4885115


# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_MM_dmrt1S_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.3806194
```
