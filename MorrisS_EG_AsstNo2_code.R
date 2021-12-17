# ------------------------------------------------------
# EG Assignment 2 script
# 25 Oct 2021
# Sarah K. Morris
# ------------------------------------------------------

# only do the installs once: 

if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")
# 
BiocManager::install("DESeq2")
# 
BiocManager::install("vsn")
BiocManager::install("hexbin") # install before loading vsn library

## Import or install the libraries that we're likely to need

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
install.packages("ggpubr")
library(ggpubr)
install.packages("wesanderson")
library(wesanderson)
library(vsn)  

setwd("/Users/sarahmorris/Documents/PhD/fall_2021/ecological_genomics/transcriptomics_module/assignment/")

# import the counts matrix
# we copied the txt file from the server to our local machine
countsTable <- read.table("DE_counts_F3.txt", header=TRUE, row.names=1)
# counts generated from Salmon, so has decimals 
# bwa and others don't produce decimal counts...

head(countsTable)
dim(countsTable)
# [1] 24362    16
# how many reads mapped to a given sample 

countsTableRound <- round(countsTable) # b/c DESeq2 doesn't like decimals
head(countsTableRound)

# import sample description table
conds <- read.delim("RT_tonsa_F3_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

########################################################

# starting here on Monday, October 11, 2021

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound), 
        cex.names = 0.5, las = 3, ylim = c(0,20000000))
# las changes orientation to vertical
# cex.names changes font size
abline(h=mean(colSums(countsTableRound)), col = "blue", lwd = 2)
# lwd changes line width

# average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # 11930.81
median(rowSums(countsTableRound)) # 2226

apply(countsTableRound, 2, mean) 
# average number of reads across all the genes, for each of the samples
# 2 in the apply function does the action across columns

apply(countsTableRound, 1, mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound, 1, mean))
#histogram of average number of reads per gene
# long tail likely b/c of highly expressed genes or PCR amplification

hist(apply(countsTableRound, 1, mean), xlim=c(0,5000), breaks=1000, xlab = "average number of reads per gene", main = "Histogram of F1 reads")

# Create a DESeq object and define the experiment designed here with the tilda

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~ line + environment + line:environment)
# the DESeqDataSetFromMatrix command is just importing the data
# conds = conditions  
# line:environment --> interaction term
dim(dds) # check that something happened

# filter out genes with too few reads 
# keep reads that have an average > 10 reads per sample(gene?)
# 10 * 16 samples = 160

dds <- dds[rowSums(counts(dds)) > 160] # oberwriting the original object, can't go back...
dim(dds) # dimensions didn't change, so try a bigger threshold maybe?

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds) # main workhorse of DESeq
# 'estimating dispersions' = looking at variation in read size

# List the results you've generated
resultsNames(dds)
# [1] "Intercept"                  "line_combined_vs_ambient"   "environment_HH_vs_AA"      
# [4] "linecombined.environmentHH"

## above is line effect, environment effect, and the interaction effect

# let's look at these patterns with a PCA
# to visualize global gene expression patterns

vsd <- vst(dds, blind=F)

# make calculations to make the PCA plot
data <- plotPCA(vsd, intgroup=c("line","environment"), returnData=T)
percentVar <- round(100 * attr(data,"percentVar")) #variance 

ggplot(data, aes(PC1,PC2, color=environment, shape=line)) +
  geom_point(size = 4, alpha=0.85) +
  xlab(paste0("PC1: ", percentVar[1],"% variance"))+
  ylab(paste0("PC2: ", percentVar[2],"% variance"))+
  theme(plot.title = element_text(size = (18), face = "bold.italic"))+
  labs(title = "F3 Acartia tonsa")

# both line and environment are having a strong effect
# PCA2 16%  variance = 16% of the variance in the data is explained by PCA 2
# PCA1 will always explain the highest percentage of the variance, PCA 2 will explain the 2nd most percent of the variance in the dataset 

########################################################
# order and summarize results from specific contrasts

# DESeq doesn't handle interactions well...

resInteraction <- results(dds, alpha = 0.05) # the false positive rate we're willing to accept
resInteraction <- resInteraction[order(resInteraction$padj),] # adjusted p-values
head(resInteraction)
# baseMean = average gene expression
# log2FoldChange = to change the expression?
# most significant genes at the top

# log2 fold change (MLE): linecombined.environmentHH 
# Wald test p-value: linecombined.environmentHH 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN115950_c0_g1   2245.97        4.05237  0.358490   11.3040 1.25396e-29 2.99446e-25
# TRINITY_DN131561_c0_g1   3375.99        4.64570  0.439847   10.5621 4.46620e-26 5.33264e-22
# TRINITY_DN137662_c0_g1  16743.23        4.90200  0.474583   10.3291 5.20658e-25 4.14444e-21
# TRINITY_DN149842_c8_g4  25971.82        4.27274  0.420809   10.1536 3.19275e-24 1.90607e-20
# TRINITY_DN129565_c0_g3  24258.76        4.30553  0.426037   10.1060 5.19661e-24 2.48190e-20
# TRINITY_DN129401_c0_g5  11712.31        4.46355  0.446094   10.0059 1.43650e-23 5.71728e-20

summary(resInteraction)
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2839, 12%
# LFC < 0 (down)     : 1053, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 473, 1.9%
# (mean count < 18)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## 12% tested had a logfold increase in expression

########################################################

# starting here Wednesday, October 13

# DESeq uses a negative binomial GLM, which is used for modeling count variables, usually for over-dispersed count outcome variables like gene expression data. For example, considering one gene, one sample may have 10,000 read counts (i.e., reads that mapped to that gene from that sample) while another sample may have 10 read counts that map to that gene.
# aka 'overdispersion'

# our count data is the reads
# the model:
# counts/reads ~ line + envi + line:envi
# aka what effect does the line, envi a/o interaction have on the observed counts/gene expression?

# DESeq has two ways to test for significance, using the Wald test (standard) or using the Likelihood Ratio Test (LRT; useful for study designs where there may be an interaction between two factors, in this case line and environment). 

# Using the LRT (likelihood ratio test), however, we have to test for each effect separately

############### TEST FOR EFFECT OF ENVIRONMENT ########

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ line + environment)

dds <- DESeq(dds, test="LRT", reduced=~line)
# can run a full or reduced LRT
# "reduced=~line" --> b/c you take our the thing that you're testing for in the reduced
# then the LRT compares the results of with and w/o the thing your testing and shoots out 
# a p-value to show if it's significant

# List the results you've generated
resultsNames(dds)

# Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05) # the false-positive rate/significance cutoff we're comfortable with
resEnv <- resEnv[order(resEnv$padj),] # order the results by the p-value (get the most significant at the top)
head(resEnv) # we don't know the function (but Melissa has the annotation table, so we could look it up)

# log2 fold change (MLE): environment HH vs AA 
# LRT p-value: '~ line + environment' vs '~ line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN138549_c1_g2    582.410       -1.94341  0.173407  118.0825 1.66325e-27 3.96651e-23
# TRINITY_DN138549_c2_g12   773.349       -2.01757  0.203760   91.4210 1.16138e-21 1.38483e-17
# TRINITY_DN150696_c2_g3    297.068        1.31754  0.163636   63.1253 1.93963e-15 1.54188e-11
# TRINITY_DN123676_c0_g2    179.431       -2.51746  0.309813   59.1190 1.48423e-14 7.07917e-11
# TRINITY_DN131329_c1_g1    213.660       -1.23500  0.158361   59.4117 1.27903e-14 7.07917e-11
# TRINITY_DN105043_c0_g1    101.714       -3.94548  0.471847   57.1227 4.09446e-14 1.62741e-10

summary(resEnv)

# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 213, 0.87% --> log fold change >0 (upregulated - higher in HH vs. AA) - based on line 196 to know the order...
# 213 genes were more highly expressed in the HH
# LFC < 0 (down)     : 235, 0.96% --> downregulation - higher in AA vs HH
# 235 genes were more highly expressed in the AA
# outliers [1]       : 41, 0.17% --> not included
# low counts [2]     : 473, 1.9% --> < 18 (# this is beyond the low counts that we filtered out earlier)
# (mean count < 18) --> definition of low counts
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

dim(resEnv)
resEnv <- resEnv[!is.na(resEnv$padj),] # filters the data table to exclude everything that is N/A in padj column
dim(resEnv)
# [1] 23848     6 --> we lost a few (the N/As) - probably the low counts or outliers

degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) # save the genes with significance
# degs = differentially expressed genes from the environment
length(degsEnv)
# [1] 448
# should be the sum of LFC >0 & LFC <0
# aka 213 + 235, which it is!

############################# TEST FOR EFFECT OF LINE ###############

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment) # does including line better describe the data?
resultsNames(dds)

resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)

# log2 fold change (MLE): line combined vs ambient 
# LRT p-value: '~ environment + line' vs '~ environment' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN144155_c0_g8    55.1773       -1.73551  0.199458   75.7313 3.25019e-18 7.13710e-14
# TRINITY_DN116533_c0_g1   101.9489        3.49591  0.435406   53.9339 2.07352e-13 2.27662e-09
# TRINITY_DN142881_c2_g11 1414.7388       -1.43499  0.199146   49.9265 1.59616e-12 1.16834e-08
# TRINITY_DN140379_c0_g5    49.4278        1.69441  0.272190   37.9057 7.42487e-10 4.07607e-06
# TRINITY_DN140379_c0_g6   220.5736        1.86590  0.297107   37.1901 1.07155e-09 4.70602e-06
# TRINITY_DN138009_c0_g2    88.4193        2.12835  0.354530   33.7159 6.37772e-09 2.00069e-05

summary(resLine)

# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 72, 0.3%
# LFC < 0 (down)     : 154, 0.63%
# outliers [1]       : 41, 0.17%
# low counts [2]     : 2362, 9.7%
# (mean count < 25)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

resLine <- resLine[!is.na(resLine$padj),]

degsline <- row.names(resLine[resLine$padj < 0.05,])
# differentially expressed genes from the line
length(degsline)
############################## TEST FOR INTERACTION ################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)

resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)

# log2 fold change (MLE): environmentHH.linecombined 
# LRT p-value: '~ environment + line + environment:line' vs '~ environment + line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN115950_c0_g1   2245.97        4.05237  0.358490  118.7642 1.17951e-27 2.76100e-23
# TRINITY_DN131561_c0_g1   3375.99        4.64570  0.439847  101.4199 7.44127e-24 8.70926e-20
# TRINITY_DN137662_c0_g1  16743.23        4.90200  0.474583   95.9225 1.19473e-22 9.32205e-19
# TRINITY_DN149842_c8_g4  25971.82        4.27274  0.420809   94.8761 2.02680e-22 1.18608e-18
# TRINITY_DN129565_c0_g3  24258.76        4.30553  0.426037   93.8792 3.35378e-22 1.57010e-18
# TRINITY_DN129401_c0_g5  11712.31        4.46355  0.446094   91.5220 1.10356e-21 4.30535e-18

summary(resInt)

# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2802, 12%
# LFC < 0 (down)     : 1052, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 945, 3.9%
# (mean count < 20)

## we expected this based on PCA from earlier

resInt <- resInt[!is.na(resInt$padj),]

degsInt <- row.names(resInt[resInt$padj < 0.05,])

length(degsInt)
########################################################
################## DATA VISUALIZATION ##################
########################################################

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN115950_c0_g1", intgroup = (c("line","environment")), returnData=TRUE) # pulled from one of the results files above
# plotCounts pulls out gene from plot counts and gives us the data
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) # not a data point, a mean
p
# we see the crossing of the lines, as expected with the interaction (b/c we pulled this from the interaction model results, if it had been from envi or line model results, would've been the parralel lines)
# this gene might get downregulated in stressful conditions

################## PLOT OVERLAPPING DEGS IN VENN DIAGRAM ###########
# We’ll use the Eulerr package because we all know how nice it is to have the circle scaled

install.packages("eulerr")
library(eulerr)

# Total
length(degsEnv)  # 448
length(degsline)  # 226
length(degsInt)  # 3854

# Intersections
length(intersect(degsEnv,degsline))  # 37
length(intersect(degsEnv,degsInt))  # 44
length(intersect(degsInt,degsline))  # 34

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7

# Number unique
448-44-37-7 # 360
226-37-34-7 # 148
3854-44-34-7 # 3769


fit1 <- euler(c("Env" = 360, "Line" = 148, "Interaction" = 3769, "Env&Line" = 37, "Env&Interaction" = 44, "Line&Interaction" = 34, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE) 
# lty = 1:3 --> gives 3 different line weights

plot(fit1, quantities = TRUE,
     lty = 1:3,
     labels = list(font = 4),
     main = "F1 Acartia tonsa")

########################################################
# As with the PCA, we’ll use the vsd function again

# vst is a transformation implemented in DESeq2, which is “roughly similar to putting the data on the log2 scale, while also dealing with the sampling variability of low counts” (according the the package manual). 

# It uses the design formula to calculate the within-group variability (if blind=FALSE) or the across-all-samples variability (if blind=TRUE)

# Heatmap of top 20 differentially expressed genes sorted by pvalue

install.packages("pheatmap")
library(pheatmap)

# By interaction

topgenes <- head(rownames(resInt),20) 
mat <- assay(vsd)[topgenes,]
# run stabilization on topgenes
# transformation for visualization
mat <- mat - rowMeans(mat) # scale them for better visualization
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)
# dendrograms cluster by similiarity of genes
# heatmap colors indicate upregulation and downregulation

# By line (this one is the cleanest heatmap of the 3...)

topgenes <- head(rownames(resLine),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

# by envi
topgenes <- head(rownames(resEnv),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

########################################################
########################################################

########################################################
# F3
########################################################
setwd("/Users/sarahmorris/Documents/PhD/fall_2021/ecological_genomics/transcriptomics_module/assignment/")

# import the counts matrix
# we copied the txt file from the server to our local machine
countsTable <- read.table("DE_counts_F3.txt", header=TRUE, row.names=1)
# counts generated from Salmon, so has decimals 
# bwa and others don't produce decimal counts...

head(countsTable)
dim(countsTable)
# [1] 24362    16
# how many reads mapped to a given sample 

countsTableRound <- round(countsTable) # b/c DESeq2 doesn't like decimals
head(countsTableRound)

# import sample description table
conds <- read.delim("RT_tonsa_F3_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

########################################################

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound), 
        cex.names = 0.5, las = 3, ylim = c(0,20000000))
# las changes orientation to vertical
# cex.names changes font size
abline(h=mean(colSums(countsTableRound)), col = "blue", lwd = 2)
# lwd changes line width

# average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # 11930.81
median(rowSums(countsTableRound)) # 2226

apply(countsTableRound, 2, mean) 
# average number of reads across all the genes, for each of the samples
# 2 in the apply function does the action across columns

apply(countsTableRound, 1, mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound, 1, mean))
#histogram of average number of reads per gene
# long tail likely b/c of highly expressed genes or PCR amplification

hist(apply(countsTableRound, 1, mean), xlim=c(0,5000), breaks=1000, xlab = "average number of reads per gene", main = "Histogram of F3 reads")

# Create a DESeq object and define the experiment designed here with the tilda

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~ line + environment + line:environment)
# the DESeqDataSetFromMatrix command is just importing the data
# conds = conditions  
# line:environment --> interaction term
dim(dds) # check that something happened

# filter out genes with too few reads 
# keep reads that have an average > 10 reads per sample(gene?)
# 10 * 16 samples = 160

dds <- dds[rowSums(counts(dds)) > 160] # oberwriting the original object, can't go back...
dim(dds) # dimensions didn't change, so try a bigger threshold maybe?

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds) # main workhorse of DESeq
# 'estimating dispersions' = looking at variation in read size

# List the results you've generated
resultsNames(dds)
# [1] "Intercept"                  "line_combined_vs_ambient"   "environment_HH_vs_AA"      
# [4] "linecombined.environmentHH"

## above is line effect, environment effect, and the interaction effect

# let's look at these patterns with a PCA
# to visualize global gene expression patterns

vsd <- vst(dds, blind=F)

# make calculations to make the PCA plot
data <- plotPCA(vsd, intgroup=c("line","environment"), returnData=T)
percentVar <- round(100 * attr(data,"percentVar")) #variance 

ggplot(data, aes(PC1,PC2, color=environment, shape=line)) +
  geom_point(size = 4, alpha=0.85) +
  xlab(paste0("PC1: ", percentVar[1],"% variance"))+
  ylab(paste0("PC2: ", percentVar[2],"% variance"))+
  theme(plot.title = element_text(size = (18), face = "bold.italic"))+
  labs(title = "F3 Acartia tonsa")

# both line and environment are having a strong effect
# PCA2 16%  variance = 16% of the variance in the data is explained by PCA 2
# PCA1 will always explain the highest percentage of the variance, PCA 2 will explain the 2nd most percent of the variance in the dataset 

########################################################
# order and summarize results from specific contrasts

# DESeq doesn't handle interactions well...

resInteraction <- results(dds, alpha = 0.05) # the false positive rate we're willing to accept
resInteraction <- resInteraction[order(resInteraction$padj),] # adjusted p-values
head(resInteraction)
# baseMean = average gene expression
# log2FoldChange = to change the expression?
# most significant genes at the top

# log2 fold change (MLE): linecombined.environmentHH 
# Wald test p-value: linecombined.environmentHH 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN115950_c0_g1   2245.97        4.05237  0.358490   11.3040 1.25396e-29 2.99446e-25
# TRINITY_DN131561_c0_g1   3375.99        4.64570  0.439847   10.5621 4.46620e-26 5.33264e-22
# TRINITY_DN137662_c0_g1  16743.23        4.90200  0.474583   10.3291 5.20658e-25 4.14444e-21
# TRINITY_DN149842_c8_g4  25971.82        4.27274  0.420809   10.1536 3.19275e-24 1.90607e-20
# TRINITY_DN129565_c0_g3  24258.76        4.30553  0.426037   10.1060 5.19661e-24 2.48190e-20
# TRINITY_DN129401_c0_g5  11712.31        4.46355  0.446094   10.0059 1.43650e-23 5.71728e-20

summary(resInteraction)
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2839, 12%
# LFC < 0 (down)     : 1053, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 473, 1.9%
# (mean count < 18)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## 12% tested had a logfold increase in expression

########################################################

# starting here Wednesday, October 13

# DESeq uses a negative binomial GLM, which is used for modeling count variables, usually for over-dispersed count outcome variables like gene expression data. For example, considering one gene, one sample may have 10,000 read counts (i.e., reads that mapped to that gene from that sample) while another sample may have 10 read counts that map to that gene.
# aka 'overdispersion'

# our count data is the reads
# the model:
# counts/reads ~ line + envi + line:envi
# aka what effect does the line, envi a/o interaction have on the observed counts/gene expression?

# DESeq has two ways to test for significance, using the Wald test (standard) or using the Likelihood Ratio Test (LRT; useful for study designs where there may be an interaction between two factors, in this case line and environment). 

# Using the LRT (likelihood ratio test), however, we have to test for each effect separately

############### TEST FOR EFFECT OF ENVIRONMENT ########

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ line + environment)

dds <- DESeq(dds, test="LRT", reduced=~line)
# can run a full or reduced LRT
# "reduced=~line" --> b/c you take our the thing that you're testing for in the reduced
# then the LRT compares the results of with and w/o the thing your testing and shoots out 
# a p-value to show if it's significant

# List the results you've generated
resultsNames(dds)

# Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05) # the false-positive rate/significance cutoff we're comfortable with
resEnv <- resEnv[order(resEnv$padj),] # order the results by the p-value (get the most significant at the top)
head(resEnv) # we don't know the function (but Melissa has the annotation table, so we could look it up)

# log2 fold change (MLE): environment HH vs AA 
# LRT p-value: '~ line + environment' vs '~ line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN138549_c1_g2    582.410       -1.94341  0.173407  118.0825 1.66325e-27 3.96651e-23
# TRINITY_DN138549_c2_g12   773.349       -2.01757  0.203760   91.4210 1.16138e-21 1.38483e-17
# TRINITY_DN150696_c2_g3    297.068        1.31754  0.163636   63.1253 1.93963e-15 1.54188e-11
# TRINITY_DN123676_c0_g2    179.431       -2.51746  0.309813   59.1190 1.48423e-14 7.07917e-11
# TRINITY_DN131329_c1_g1    213.660       -1.23500  0.158361   59.4117 1.27903e-14 7.07917e-11
# TRINITY_DN105043_c0_g1    101.714       -3.94548  0.471847   57.1227 4.09446e-14 1.62741e-10

summary(resEnv)

# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 213, 0.87% --> log fold change >0 (upregulated - higher in HH vs. AA) - based on line 196 to know the order...
# 213 genes were more highly expressed in the HH
# LFC < 0 (down)     : 235, 0.96% --> downregulation - higher in AA vs HH
# 235 genes were more highly expressed in the AA
# outliers [1]       : 41, 0.17% --> not included
# low counts [2]     : 473, 1.9% --> < 18 (# this is beyond the low counts that we filtered out earlier)
# (mean count < 18) --> definition of low counts
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

dim(resEnv)
resEnv <- resEnv[!is.na(resEnv$padj),] # filters the data table to exclude everything that is N/A in padj column
dim(resEnv)

degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) # save the genes with significance
# degs = differentially expressed genes from the environment
length(degsEnv)

# should be the sum of LFC >0 & LFC <0


############################# TEST FOR EFFECT OF LINE ###############

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment) # does including line better describe the data?
resultsNames(dds)
resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)
summary(resLine)
resLine <- resLine[!is.na(resLine$padj),]
degsline <- row.names(resLine[resLine$padj < 0.05,])
length(degsline)
############################## TEST FOR INTERACTION ################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)

resultsNames(dds)
resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)
summary(resInt)
resInt <- resInt[!is.na(resInt$padj),]
degsInt <- row.names(resInt[resInt$padj < 0.05,])
length(degsInt)

########################################################
################## DATA VISUALIZATION ##################
########################################################

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN115950_c0_g1", intgroup = (c("line","environment")), returnData=TRUE) # pulled from one of the results files above
# plotCounts pulls out gene from plot counts and gives us the data
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) # not a data point, a mean
p
# we see the crossing of the lines, as expected with the interaction (b/c we pulled this from the interaction model results, if it had been from envi or line model results, would've been the parralel lines)
# this gene might get downregulated in stressful conditions


################## PLOT OVERLAPPING DEGS IN VENN DIAGRAM ###########
# We’ll use the Eulerr package because we all know how nice it is to have the circle scaled

#install.packages("eulerr")
library(eulerr)

# Total
length(degsEnv)  # 828
length(degsline)  # 1645
length(degsInt)  # 283

# Intersections
length(intersect(degsEnv,degsline))  # 141
length(intersect(degsEnv,degsInt))  # 14
length(intersect(degsInt,degsline))  # 32

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7

# Number unique
828-141-14-7 # 666
1645-141-32-7 # 1465
283-14-32-7 # 230


fit3 <- euler(c("Env" = 666, "Line" = 1465, "Interaction" = 230, "Env&Line" = 141, "Env&Interaction" = 14, "Line&Interaction" = 32, "Env&Line&Interaction" = 7))

plot(fit3,  lty = 1:3, quantities = TRUE) 
# lty = 1:3 --> gives 3 different line weights

plot(fit3, quantities = TRUE,
     lty = 1:3,
     labels = list(font = 4),
     main = "F3 Acartia tonsa")

########################################################
# As with the PCA, we’ll use the vsd function again

# vst is a transformation implemented in DESeq2, which is “roughly similar to putting the data on the log2 scale, while also dealing with the sampling variability of low counts” (according the the package manual). 

# It uses the design formula to calculate the within-group variability (if blind=FALSE) or the across-all-samples variability (if blind=TRUE)

# Heatmap of top 20 differentially expressed genes sorted by pvalue

# install.packages("pheatmap")
library(pheatmap)

# By interaction

topgenes <- head(rownames(resInt),20) 
mat <- assay(vsd)[topgenes,]
# run stabilization on topgenes
# transformation for visualization
mat <- mat - rowMeans(mat) # scale them for better visualization
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)
# dendrograms cluster by similiarity of genes
# heatmap colors indicate upregulation and downregulation

# By line (this one is the cleanest heatmap of the 3...)

topgenes <- head(rownames(resLine),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

# by envi
topgenes <- head(rownames(resEnv),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

########################################################
########################################################

