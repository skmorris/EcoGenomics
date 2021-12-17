# ------------------------------------------------------
# EG Assignment 3 analysis
# 22 Nov 2021
# Sarah K. Morris
# ------------------------------------------------------

# after command line analysis & downloading files to my local machine

setwd("/Users/sarahmorris/Documents/PhD/fall_2021/ecological_genomics/landscape_module/homework3")

library(ggplot2)
library(gridExtra)

# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Get the meta data:
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)

# merge them together:
meta_admx <- merge(meta, Admixed, by.x="ID", by.y="V1")
str(meta_admx)  

# Read in the Admixture coefficients for KBals that we made from the K=5 file:
KBals <- read.table("Admixed_KBals", sep="\t", header=F)
names(KBals) = c("ID","KBals")

# Second merge:
meta_admx_KBals <- merge(meta_admx,KBals,by="ID")


# Bring in climate data:
clim <- read.table("climDat.txt",sep="\t",header=T)

# Merge climate data with meta and KBals:
meta_admx_KBals_clim <- merge(meta_admx_KBals,clim,by="ID")
head(meta_admx_KBals_clim)


########################################################

#  mean final freeze
plotFF <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=mean_finalFreeze,
                                          color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("mean final freeze") 

plotFF

# average number of growing degree days
plotgDD <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=mean_cGDDfreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Avg. # of growing degree days") 

plotgDD

# mean number of chilling degree days across the year
plotDD0 <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=med_DD0, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Avg. chilling degree days ") 

plotDD0

grid.arrange(plotFF, plotgDD, plotDD0, nrow = 3)

# linear models testing climate variable ~ genome-wide admixture association
summary(lm(mean_finalFreeze~KBals + Transect.x, data=meta_admx_KBals_clim))

summary(lm(mean_cGDDfreeze~KBals + Transect.x, data=meta_admx_KBals_clim))

summary(lm(med_DD0~KBals + Transect.x, data=meta_admx_KBals_clim))

########################################################

######  Bring in Association results from Plink   ######

FF <- read.table("plink2.FF.glm.linear",skip=1,sep="\t",header=F)
names(FF) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
head(FF)
# each row is a SNP
# 2 tests: beta1, beta2

FF2 <- FF[which(FF$TEST=="ADD"),]

# Define association outliers as the upper 1% of p-values
# make snps
snps <- read.table("Chr01.kept.sites", sep="\t", header=T)

#########  average date of the last freezing event in spring (FF)  #########
FF2 <- cbind(snps, FF2[,-c(1:2)])
FF2$outlier = ifelse(FF2$P<quantile(FF2$P,0.01),2,1)

p1 <- ggplot(FF2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=FF2$outlier, color=FF2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Chromosome 1 avg. final freeze")

p1
# higher y-axis values are more significant
# the red zones are the first clue for where we should look for differential gene expression 
# which genes specifically 

#######  average number of growing degree days  #########
gDD <- read.table("plink2.GDD.glm.linear",skip=1,sep="\t",header=F)
names(gDD) = c("CHROM",  "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
gDD2 <- gDD[which(gDD$TEST=="ADD"),]
gDD2 <- cbind(snps, gDD2[,-c(1,2)])
gDD2$outlier = ifelse(gDD2$P<quantile(gDD2$P,0.01),2,1)

p2 <- ggplot(gDD2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=gDD2$outlier, color=gDD2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Chromosome 1 avg. # growing degree days")

p2

#########  mean number of chilling degree days across the year  #########
DD0 <- read.table("Plink2.DD0.glm.linear",skip=1,sep="\t",header=F)
names(DD0) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
DD0 <- DD0[which(DD0$TEST=="ADD"),]
DD02 <- cbind(snps, DD0[,-c(1,2)])
DD02$outlier = ifelse(DD02$P<quantile(DD02$P,0.01),2,1)

p3 <- ggplot(DD02,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=DD02$outlier, color=DD02$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Chromosome 1 avg. chilling degree days")

p3

grid.arrange(p1, p2, p3, nrow = 3)

########################################################

# Get outliers for a given trait association:

FF_outliers <- FF2[which(FF2$outlier==2),c(2,3,9)]
gDD_outliers <- gDD2[which(gDD2$outlier==2),c(2,3,9)]
DD0_outliers <- DD02[which(DD02$outlier==2),c(2,3,9)]

########################################################

# Read in list of positions
snps <- read.table("Chr01.kept.sites",sep="\t", header=T)

# Plot freq of LAI along chr
AF <- read.table("Chr01_LAI_freq.afreq", skip=1,sep="\t",header=F)
names(AF) = c("CHROM",  "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")
str(AF)

AF2 <- cbind(snps,AF)

windows <- seq(1,max(AF2$POS),5e4)
AF_windows <- numeric()

for(i in 1:length(windows)){
  tmp=AF2[which(AF2$POS>windows[i] & AF2$POS<windows[i+1]),"ALT_FREQS"]
  ancfreq=mean(tmp)
  AF_windows[i] = ancfreq
}

AF3 <- as.data.frame(cbind(windows,AF_windows))
names(AF3) = c("window","AvgAncFreq")

########################################################

upper = mean(AF3$AvgAncFreq,na.rm=T) + 2*sd(AF3$AvgAncFreq,na.rm=T)
lower = mean(AF3$AvgAncFreq,na.rm=T) - 2*sd(AF3$AvgAncFreq,na.rm=T)

outliers_upper = AF3[which(AF3$AvgAncFreq>upper),]
outliers_lower = AF3[which(AF3$AvgAncFreq<lower),]

# Print the outlier regions out
outliers_upper
outliers_lower

# make the 4-panel plot with the climate associations
p4 <- ggplot(AF3[,-3],aes(x=window,y=AvgAncFreq)) +
  geom_line(size=0.8, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("freq. trichocarpa ancestry") +
  geom_hline(yintercept=mean(AF2$ALT_FREQS), color = "red") + 
  geom_hline(yintercept=upper, linetype="dashed", color = "red") + 
  geom_hline(yintercept=lower, linetype="dashed", color = "red") +
  ggtitle("Chr01: Local ancestry")

p4


grid.arrange(p1, p2, p3, p4, nrow = 4)

########################################################

#########  Bud flush  #########
budflush <- read.table("Plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.01),2,1)

p6 <- ggplot(budflush2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=budflush2$outlier, color=budflush2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Chromosome 1 Bud flush")

p6

grid.arrange(p1, p2, p3, p6, p4, nrow = 5)

########################################################
######### QUESTION 2 ###################################
########################################################
# look for areas of chromosomal overlap b/n flush and climate
########################################################

library(GenomicRanges)
library(GenomicFeatures)

CHR="Chr01"

# define the window
budflushGR <- GRanges(CHR,IRanges(budflush2$POS-2.5e4,budflush2$POS+2.5e4),POS=budflush2$POS, P=budflush2$P, outlier=budflush2$outlier)

budflushGRout <- unlist(reduce(split(budflushGR, ~outlier)))
budflushGRout$outlier <- names(budflushGRout)
budflushGRCand <- subset(budflushGRout, outlier==2)

budflushGRCand # Print the candidate regions

########################################################
#FF
FF_GR <- GRanges(CHR,IRanges(FF2$POS-2.5e4,FF2$POS+2.5e4),POS=FF2$POS, P=FF2$P, outlier=FF2$outlier)

FF_GRout <- unlist(reduce(split(FF_GR, ~outlier)))
FF_GRout$outlier <- names(FF_GRout)
FF_GRCand <- subset(FF_GRout, outlier==2)

FF_GRCand # Print the candidate regions
########################################################
#gDD
gDD_GR <- GRanges(CHR,IRanges(gDD2$POS-2.5e4,gDD2$POS+2.5e4),POS=gDD2$POS, P=gDD2$P, outlier=gDD2$outlier)

gDD_GRout <- unlist(reduce(split(gDD_GR, ~outlier)))
gDD_GRout$outlier <- names(gDD_GRout)
gDD_GRCand <- subset(gDD_GRout, outlier==2)

gDD_GRCand # Print the candidate regions
########################################################
#DD0
DD0_GR <- GRanges(CHR,IRanges(DD02$POS-2.5e4,DD02$POS+2.5e4),POS=DD02$POS, P=DD02$P, outlier=DD02$outlier)

DD0_GRout <- unlist(reduce(split(DD0_GR, ~outlier)))
DD0_GRout$outlier <- names(DD0_GRout)
DD0_GRCand <- subset(DD0_GRout, outlier==2)

DD0_GRCand # Print the candidate regions
########################################################

# overlap b/n bud flush and final freeze
overlap_BF_FF <- subsetByOverlaps(budflushGRCand, FF_GRCand)
length(overlap_BF_FF) # 0
overlap_BF_FF # Print the overlapping regions

########################################################

# overlap b/n bud flush and gDD
overlap_BF_gDD <- subsetByOverlaps(budflushGRCand, gDD_GRCand)
length(overlap_BF_gDD) # 5
overlap_BF_gDD # Print the overlapping regions

########################################################

# overlap b/n bud flush and DD0
overlap_BF_DD0 <- subsetByOverlaps(budflushGRCand, DD0_GRCand)
length(overlap_BF_DD0) # 6
overlap_BF_DD0 # Print the overlapping regions

########################################################
########### QUESTION 3: functional enrichment ##########
########################################################

# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")
txdb

# How many chromosomes are present?
head(seqlevels(txdb))

seqlevels(txdb) <- CHR # subset for just chromosome 1

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) # collapse to make it one per gene
genes$geneID <- names(genes)

#FF
BFxFFcandGenes <- subsetByOverlaps(genes, overlap_BF_FF)
BFxFFcandGenes # 0 candidate genes

#gDD
BFxgDDcandGenes <- subsetByOverlaps(genes, overlap_BF_gDD)
BFxgDDcandGenes # 72 candidate genes

write.table(BFxgDDcandGenes$geneID, paste0("BFxgDDcandGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

#DD0
BFxDD0candGenes <- subsetByOverlaps(genes, overlap_BF_DD0)
BFxDD0candGenes # 78 candidate genes

write.table(BFxDD0candGenes$geneID, paste0("BFxDD0candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

########################################################

# open up the text file of candidate genes, copy the list, and then make a new “active” gene list on Popgenie
# then go to Analysis Tools > Enrichment

