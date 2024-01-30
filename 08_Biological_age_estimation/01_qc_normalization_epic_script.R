# QC and normalization script of EPIC array for the paper: "Cell-type-specific 
# ageing effects in the human OFC and implications for psychiatric disease" by 
# Fröhlich et al.
# ==============================================================================

# Author: Natan Yusupov, natan_yusupov@psych.mpg.de
# Supervision: Dr. Darina Czamara

# Raw data was preprocessed and normalized differently for each of the epigenetic
# clocks. Either using A. modified pipeline proposed by Maksimovic et al. 
# (Maksimovic J, Phipson B, Oshlack A. A cross-package Bioconductor workflow for 
# analysing methylation array data. F1000Res. 2016) or B. Bigmelon pipeline 
# by Gorrie-Stone et al. (Gorrie-Stone TJ, Smart MC, Saffari A, Malki K, Hannon E, 
# Burrage J, et al. Bigmelon: tools for analysing large DNA methylation datasets. 
# Bioinformatics 2019

################################################################################

# Load libraries
library(minfi)
library(minfiData)
library(minfiDataEPIC)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(sva)
library(RPMM)
library(ggplot2)
library(missMethyl)
library(matrixStats)
library(wateRmelon)
library(bigmelon)
library(reshape)
library(Hmisc)
library(tidyverse)
library(readxl)
library(magrittr)
library(tibble)

# Part A 
# Generate RGSet and phenotype data with iDAT files and sample sheet
targets <- read_xlsx("samplesheet_EPIC.xlsx")
RGSet	<- read.metharray.exp(targets = targets)
pd_orig <- pData(RGSet)

# Get minfi data types
Mset <- preprocessRaw(RGSet)
RatioSet <- ratioConvert(Mset, what = "both", keepCN = TRUE)
RawBetas <- getBeta(RatioSet)
gRatioSet <- mapToGenome(RatioSet,mergeManifest=TRUE) 

# QC
# Detection p values
detP <- detectionP(RGSet)

# Examine mean detection p-values across all samples 
pdf("detectionP.pdf", width = 8, height = 3)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Name)], las=2, cex.names=0.4,ylab="Mean detection p-values") 
abline(h=0.01,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal, bg="white") 
barplot(colMeans(detP), col=pal[factor(targets$Sample_Name)], las=2, cex.names=0.4, ylim = c(0,0.002), ylab="Mean detection p-values") 
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal, bg="white") 
dev.off()

# Remove poor quality samples (mean detection p-value>.05)
keep <- colMeans(detP) < 0.05
RGSet_qual <- RGSet[,keep]
RawBetas_qual <- RawBetas[,keep]
targets_qual <- targets[keep,] 
targets_qual[,c(1:5,8:17)]
detP_qual <- detP[,keep]

# Run minfi QC Report
qcReport(RGSet_qual, sampGroups=targets_qual$Slides, pdf="qcReport.pdf") 

# Remove distribution artefacts of individual densities (visual exploration)
pdf("beta_densities.pdf")
for (i in 1:ncol(RGSet_qual))
{titel<-paste(rownames(pData(RGSet_qual))[i])
densityPlot(as.matrix(RawBetas_qual[,i]),main=titel)
print(i)}
dev.off()

# Remove sex mismatches (predict sex with DNAm data and reported sex)
predictedSex<-getSex(gRatioSet,cutoff=-2)
sex<-cbind(sampleNames(RGSet),as.character(pd_orig$Sample_Name),
           predictedSex$predictedSex,predictedSex$yMed,predictedSex$xMed,
           predictedSex$yMed-predictedSex$xMed)
sex<-as.data.frame(sex)
gender <- targets %>% mutate(ArrayID = paste0(Slide, "_", Array)) %>%
                               dplyr::select(Sample_Name, ArrayID, gender)
names(sex)<-c("ArrayID", "predictedSex","yMed","xMed","yMed-xMed")
sex$yMed<-as.numeric(as.character(sex$yMed))
sex$xMed<-as.numeric(as.character(sex$xMed))
test<-merge(sex,gender,by.x="ArrayID",by.y="ArrayID")
NoIDMatchSexGender <- subset(sex$ArrayID, !(sex$ArrayID %in% gender$ArrayID))
table(test$Gender,test$predictedSex)
pd_qual<-pData(RGSet_qual)
annot <- getAnnotation(RGSet_qual)

# Stratified quantile normalization followed by BMIQ
source("addFiles/BMIQ_1.6_Teschendorff.R")
quantileN <- preprocessQuantile(RGSet_qual)
quantileNBetas = getBeta(quantileN)
quantileNMs <- getM(quantileN)
probeType <- as.data.frame(annot[rownames(quantileNBetas),c("Name","Type")])
probeType$probeType <- ifelse(probeType$Type %in% "I",1,2)
BMIQ.quantileN <- apply(quantileNBetas[,1:length(colnames(quantileNBetas))],2,function(a) BMIQ(a,probeType$probeType,plots=FALSE)$nbeta)
length(which(is.nan(BMIQ.quantileN)))

# Check distributions before and after normalization
pdf("BetaValue_Distributions_Norm_Quantile.pdf")
densityPlot(RawBetas_qual, sampGroups = pd_qual$Slide, legend=FALSE, main = "Raw Betas", xlab = "Beta", ylim=c(0,5))
densityPlot(quantileNBetas, sampGroups = pd_qual$Slide, legend=FALSE, main = "Quantile Adjusted Betas", xlab = "Beta", ylim=c(0,5))
densityPlot(BMIQ.quantileN, sampGroups = pd_qual$Slide, legend=FALSE, main = "Quantile-BMIQ Adjusted Betas", xlab = "Beta", ylim=c(0,5))
dev.off()

# Check for potential outliers in the distribution
names_RawBetas <- colnames(RawBetas_qual)
pdf("Beta_Densities_RawBetas.pdf")
for (i in 1:ncol(RawBetas_qual)) {
  i_mat <- as.matrix(RawBetas_qual[ ,i])
  densityPlot(i_mat, main=names_RawBetas[i])
  }
dev.off()

names_quantileBetas <- colnames(quantileNBetas)
pdf("Beta_Densities_QuantileBetas.pdf")
for (i in 1:ncol(quantileNBetas)) {
  i_mat <- as.matrix(quantileNBetas[ ,i])
  name <- colnames(quantileNBetas[,i])
  densityPlot(i_mat, main=names_quantileBetas[i])
}
dev.off()

# Filtering probes
# Remove probes that failed in one or more samples based on detection p < .01
detP_clean_f <- detP_qual[match(rownames(BMIQ.quantileN),rownames(detP_qual)),]
detP_clean_ff <- detP_qual[match(featureNames(quantileN),rownames(detP_qual)),]
keep <- rowSums(detP_clean_f < 0.01) == ncol(BMIQ.quantileN)
keep_ff <- rowSums(detP_clean_ff < 0.01) == ncol(quantileN)
BMIQ.quantileN_filtered <- BMIQ.quantileN[keep,]
quantileN_filtered <- quantileN[keep_ff,]

# Remove probes on sex chromosomes
keep <- !(rownames(BMIQ.quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,]
keep_ff <- !(featureNames(quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
quantileN_filtered <- quantileN_filtered[keep_ff,]

# Remove probes where common SNPs may affect the CpG, cross hybradizing and polymorph probes 
quantileN_filtered <- dropLociWithSnps(quantileN_filtered)
BMIQ.quantileN_filtered <-  BMIQ.quantileN_filtered[rownames(BMIQ.quantileN_filtered) %in% featureNames(quantileN_filtered),]
load("addFiles/ChenProbeIDs.rdata") # Data from Chen et al. 2013 PMID: 23314698 & McCartney et al. PMID: 27330998
annot2$SNPs <- annot2[,"EUR"]
index<-which(annot2$sex=="Exclude" | annot2$CRSH=="Exclude" | annot2$EUR=="Exclude")
exclude_Chen <-annot2[index,]
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_Chen$Name)
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,]
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_Chen$Name)
quantileN_filtered <- quantileN_filtered[keep_ff,]
exclude_crosshyb <- read.table("addFiles/CpGs_crosshybridizing_EPIC.txt",sep="",header=F)
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_crosshyb$V1)
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,]
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_crosshyb$V1)
quantileN_filtered <- quantileN_filtered[keep_ff,]
exclude_poly <- read.table("addFiles/CpGs_polymorphic_EPIC.txt",sep="",header=T)
index<-which(exclude_poly$EUR_AF>0.05)
exclude_polym <- exclude_poly[index,]
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_polym$IlmnID)
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,]
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_polym$IlmnID)
quantileN_filtered <- quantileN_filtered[keep_ff,]
Betas_quantileN_filtered <- getBeta(quantileN_filtered)
Ms_quantileN_filtered <- getM(quantileN_filtered)

# Check density plots after excluding the poorly-detected probes
pdf("BetaValue_Distributions_Norm_quantile_Filter.pdf")
densityPlot(BMIQ.quantileN_filtered,sampGroups = pd_qual$Slide, legend=FALSE, main = "Post Filtering - Normalized Beta", xlab = "Beta")
dev.off()

# Remove batch effects with Combat
# Check with PCA variation in array data associated with technical batches for 
# Slide, Array, Sample_Group, row, column, Brain_ph, Hemisphere, PMI, Freezer_storage
rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
} # function will calculate the variance of each row

mval <- apply(BMIQ.quantileN_filtered, 2, function(x) log2((x)/(1-x))) # calculate M values

vars <- as.matrix(rowVars(mval)) ## Calculate probe variance and remove if any of 0 prior to Combat
vars[vars == 0] <- NA
vars <- na.omit(vars)
intersect <- intersect(rownames(vars), rownames(mval))
BMIQ.quantileN_filtered_batch <- BMIQ.quantileN_filtered[intersect,]
mval <- mval[intersect,]
table(ifelse(rownames(pd_qual) == colnames(mval),"Match","Off"))

PCobj = prcomp(t(mval), retx = T, center = T, scale. = T)
boxplot(PCobj$x,col="grey",frame=F)
plot(PCobj,type="line",cex.lab=1.5, cex.main=1.5) 
PCs = PCobj$x
R = 4 # comment: chosen according to propvar, cummvar and screenplot of PCA
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R]
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R]

pdf("PCA_variance_explained.pdf")
par(mfrow=c(1,2))
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cumulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
dev.off()

PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_qual, by = "row.names",all=T)

pdf(file="PC_Variation_by_batch.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Group), xlab = "PC1", ylab = "PC2", main="SampleGroup")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

# Check for extreme outliers
o1 <- 3*sd(Prin.comp$PC1)
o2 <- 3*sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 & abs(Prin.comp$PC2) > o2) # comment: non

# Test for batch effects via ANOVA
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC1~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC2~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC3~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC4~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC1~Prin.comp$PMI))
anova(lm(Prin.comp$PC2~Prin.comp$PMI))
anova(lm(Prin.comp$PC3~Prin.comp$PMI))
anova(lm(Prin.comp$PC4~Prin.comp$PMI))
anova(lm(Prin.comp$PC1~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC2~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC3~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC4~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Row)))

# First correction for Array
# Comment: correcting first for Slide still left batch effects in array, but not the other way around, therefore
# we chose Array as first correction
mod <- model.matrix(~1, data=pd_qual)
M_combat_1array = ComBat(mval,batch = pd_qual$Array, mod = mod)

# Check for successful removal of batch effect
PCobj = prcomp(t(M_combat_1array), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs = PCs[,1:R]
Prin.comp <- merge(PCs,pd_qual, by = "row.names",all=T)
pdf(file="PC_Variation_by_batch_afterCombat1array.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Group), xlab = "PC1", ylab = "PC2", main="SampleGroup")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC1~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC2~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC3~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC4~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC1~Prin.comp$PMI))
anova(lm(Prin.comp$PC2~Prin.comp$PMI))
anova(lm(Prin.comp$PC3~Prin.comp$PMI))
anova(lm(Prin.comp$PC4~Prin.comp$PMI))
anova(lm(Prin.comp$PC1~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC2~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC3~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC4~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Row)))

# Second correction for Row
mod <- model.matrix(~1, data=pd_qual)
M_combat_2row = ComBat(M_combat_1array,batch = pd_qual$Row, mod = mod)

# Check for successful removal of batch effect
PCobj = prcomp(t(M_combat_2row), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs = PCs[,1:R]
Prin.comp <- merge(PCs,pd_qual, by = "row.names", all=T)

pdf(file="PC_Variation_by_batch_afterCombat2_row.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Group), xlab = "PC1", ylab = "PC2", main="Sample_Group")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC1~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC2~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC3~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC4~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC1~Prin.comp$PMI))
anova(lm(Prin.comp$PC2~Prin.comp$PMI))
anova(lm(Prin.comp$PC3~Prin.comp$PMI))
anova(lm(Prin.comp$PC4~Prin.comp$PMI))
anova(lm(Prin.comp$PC1~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC2~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC3~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC4~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Row)))

# comment: batch effect of brain_ph and freezer_storage left - to be included as covariates

# Convert batch-adjusted normalized and filtered M-values back into betas:
expit2 = function(x) 2^x/(1+2^x)
Betas_combated = expit2(M_combat_2row)

# Plot final densities
pdf("BetaValue_Distributions_afterNormCombat.pdf")
densityPlot(Betas_combated, sampGroups = pd_qual$Slide, legend=FALSE, main = "PostQC - Normalized and Batch Corrected Beta", xlab = "Beta")
dev.off()
all.equal(colnames(Betas_combated),rownames(pd_qual))
annotated_pd_qual <- new("AnnotatedDataFrame", data= as.data.frame(pd_qual)) # required for ExpressionSet
Betas_combated_ExprSet <- new("ExpressionSet", exprs= as.matrix(Betas_combated), phenoData=annotated_pd_qual)

# Remove batch effects of unfiltered data set with Combat (all steps repeated, see comments above)
mval2 <- apply(BMIQ.quantileN, 2, function(x) log2((x)/(1-x)))
vars2 <- as.matrix(rowVars(mval2))
vars2[vars2 == 0] <- NA
vars2 <- na.omit(vars2)
intersect2 <- intersect(rownames(vars2), rownames(mval2))
BMIQ.quantileN_batch <- BMIQ.quantileN[intersect2,]
mval2 <- mval2[intersect2,]

PCobj2 = prcomp(t(mval2), retx = T, center = T, scale. = T)
boxplot(PCobj2$x,col="grey",frame=F)
plot(PCobj2,type="line",cex.lab=1.5, cex.main=1.5)
PCs2 = PCobj2$x
R = 4
propvar = summary(PCobj2)$importance["Proportion of Variance", 1:R]
cummvar = summary(PCobj2)$importance["Cumulative Proportion", 1:R]

pdf("PCA_variance_explained_2.pdf")
par(mfrow=c(1,2))	
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cumulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
dev.off()

PCs2 = PCobj2$x
PCs2 = PCs2[,1:R]
Prin.comp2 <- merge(PCs2,pd_qual, by = "row.names",all=T) 

pdf(file="PC_Variation_by_batch_2.pdf")
par(mfrow=c(3,1))
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Sample_Group), xlab = "PC1", ylab = "PC2", main="Sample_Group")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

o1 <- 3*sd(Prin.comp2$PC1)
o2 <- 3*sd(Prin.comp2$PC2)
which(abs(Prin.comp2$PC1) > o1 & abs(Prin.comp2$PC2) > o2) # comment: non

anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC1~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC2~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC3~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC4~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC1~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC2~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC3~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC4~Prin.comp2$PMI))
anova(lm(Prin.comp$PC1~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC2~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC3~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC4~Prin.comp$Freezer_storage))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Array))) 
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Row)))

# First correction for Array
mod2 <- model.matrix(~1, data=pd_qual)
M_combat2_1array = ComBat(mval2,batch = pd_qual$Array, mod = mod2)
PCobj2 = prcomp(t(M_combat2_1array), retx = T, center = T, scale. = T)
PCs2 = PCobj2$x
PCs2 =PCs2[,1:R]
Prin.comp2 <-merge(PCs2,pd_qual, by = "row.names",all=T) 

pdf(file="PC_Variation_by_batch_2afterCombat1_array.pdf")
par(mfrow=c(3,1))
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Sample_Group), xlab = "PC1", ylab = "PC2", main="Sample_Group")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC1~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC2~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC3~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC4~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC1~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC2~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC3~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC4~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC1~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC2~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC3~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC4~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Array))) 
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Row)))

# Second correction for Row
mod2 <- model.matrix(~1, data=pd_qual)
M_combat2_2row = ComBat(M_combat2_1array,batch = pd_qual$Row, mod = mod2)
PCobj2 = prcomp(t(M_combat2_2row), retx = T, center = T, scale. = T)
PCs2 = PCobj2$x
PCs2 =PCs2[,1:R]
Prin.comp2 <- merge(PCs2,pd_qual, by = "row.names",all=T)

pdf(file="PC_Variation_by_batch2_afterCombat2_row.pdf")
par(mfrow=c(3,1))
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Sample_Group), xlab = "PC1", ylab = "PC2", main="Sample_Group")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Hemisphere)))
anova(lm(Prin.comp2$PC1~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC2~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC3~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC4~Prin.comp2$Brain_ph))
anova(lm(Prin.comp2$PC1~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC2~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC3~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC4~Prin.comp2$PMI))
anova(lm(Prin.comp2$PC1~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC2~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC3~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC4~Prin.comp2$Freezer_storage))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Slide)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Array)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Sample_Group)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Column)))
anova(lm(Prin.comp2$PC1~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC2~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC3~as.factor(Prin.comp2$Row)))
anova(lm(Prin.comp2$PC4~as.factor(Prin.comp2$Row)))

# comment: batch effect of brain_ph and freezer_storage left - to be included as covariates

# Convert batch-adjusted normalized M-values back into betas
expit2 <- function(x) 2^x/(1+2^x)
BMIQ.quantileN_combated <- expit2(M_combat2_2row)

# Plot final densities
pdf("BetaValue_Distributions_afterNormCombat2.pdf")
densityPlot(BMIQ.quantileN_combated, sampGroups = pd_qual$Slide, legend=FALSE, main = "PostQC - Normalized and Batch Corrected Beta", xlab = "Beta")
dev.off()
all.equal(colnames(BMIQ.quantileN_combated),rownames(pd_qual))
BMIQ.quantileN_combated_ExprSet = new("ExpressionSet", exprs= as.matrix(BMIQ.quantileN_combated), phenoData=annotated_pd_qual)

# Make RGSet with final samples
RGSet_cleaned_final <- RGSet_qual[ , (colnames(RGSet_qual) %in% sampleNames(Betas_combated_ExprSet))]

# Part B

# Steps performed as in Shireby et al. Brain 2020:
# (Shireby, G.L., et al. Recalibrating the epigenetic clock: implications for 
# assessing biological age in the human cortex. Brain 2020)
# 1. pfilter() function of wateRmelon package to exclude samples with >1% of probes 
# with detection P-value >0.05 and probes with >1% of samples with detection P-value >0.05
# 2. Exclude outliers with PCA
# 3. Removal of cross-hybridizing and SNP probes (Chen et al.,2013)
# 4. Normalization of the DNAm data with the dasen() function of bigmelon package

# Read in iDAT files to create gds.class file
gfile <- iadd2("/binder/common/methylation/raw_methylation/Projekt_M01085_brainbank_NSW/idat/", gds = "gfile.gds")

# Exclude samples with >1% of probes with detection P-value >0.05 & probes with 
# >1% of samples with detection P-value >0.05
pfilter(gfile)
# Removal report:
# 0 samples having1% of sites with a detection p-value greater than 0.05 were removed.
# 4514 sites were remove as beadcount <3 in 5% of samples.
# 9599 sites having 1% of samples with a detection p-value greater than 0.05 were removed.

# Perform normalization with dasen function
dasen(gfile, node = "betas_normalized") # store data in new node called "betas_normalized"

# Extract raw betas after first filtering
node <- betas(gfile) # first extract node
betas <- anode[ , ] # now extract all elements

# Check distributions before normalization
pdf("BetaValue_Distributions_raw.pdf")
densityPlot(betas, sampGroups = pd_qual$Slide, legend=FALSE, main = "Raw Betas", xlab = "Beta", ylim=c(0,5))
dev.off()

# Remove outlier samples in PCA
PCobj = prcomp(t(betas), retx = T, center = T, scale. = T)
PCs = PCobj$x
plot(PCobj,type="line",cex.lab=1.5, cex.main=1.5)
R = 5 # comment: chosen according to propvar, cummvar and scree plot of PCA
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R]
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R]
PCs = PCs[,1:R]

# Check for extreme outliers
o1 <- 3*sd(PCs[,1])
o2 <- 3*sd(PCs[,2])
which(abs(PCs[,1]) > o1 & abs(PCs[,2]) > o2) # comment: non

# Extract and save normalized betas as data.matrix
anode <- index.gdsn(gfile, 'betas_normalized') # first extract node
betas_normalized <- anode[ , ] # now extract all elements 
colnames(betas_normalized) <- str_sub(colnames(betas_normalized),14) # Correct colnames 
# by removing first part of colnames since doubled before merging with pData

# Check distributions after normalization
pdf("BetaValue_Distributions_norm_dasen.pdf")
densityPlot(betas_normalized, sampGroups = pd_qual$Slide, legend=FALSE, main = "bigmelon adjusted Betas", xlab = "Beta", ylim=c(0,5))
dev.off()

# Exclude probes according to Chen et al. 2013
keep <- !(rownames(betas_normalized) %in% exclude_Chen$Name)
betas_normalized_filtered <- betas_normalized[keep,]

# Check distributions after normalization and filtering
pdf("BetaValue_Distributions_norm_dasen_filtered.pdf")
densityPlot(betas_normalized_filtered, sampGroups = pd_qual$Slide, legend=FALSE, main = "Post Filtering - Normalized Beta", xlab = "Beta", ylim=c(0,5))
dev.off()

closefn.gds(gfile)

# Remove batch effects with Combat
# Check with PCA variation in array data associated with technical batches for 
# Slide, Array, Sample_Group, row, column, Brain_ph, Hemisphere, PMI, Freezer_storage
## function that will calculate the variance of each row
mval <- apply(betas_normalized_filtered, 2, function(x) log2((x)/(1-x))) # M values
vars = as.matrix(rowVars(mval))
# change order of rows for pd_qual according to colnames of mval
colnames <- colnames(mval)
epic_idx <- match(rownames(pd_qual), colnames(mval))
mval_ordered  <- mval[,epic_idx]
all(rownames(pd_qual) == colnames(mval_ordered))

PCobj = prcomp(t(mval_ordered), retx = T, center = T, scale. = T)
boxplot(PCobj$x,col="grey",frame=F)
plot(PCobj,type="line",cex.lab=1.5, cex.main=1.5) 
PCs = PCobj$x
R = 4 # comment: chosen according to propvar, cummvar and screenplot of PCA
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R]
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R]

pdf("PCA_variance_explained_normalized_filtered.pdf")
par(mfrow=c(1,2))
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cumulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
dev.off()

PCs = PCs[,1:R]
Prin.comp<-merge(PCs,pd_qual, by = "row.names",all=T)

pdf(file="PC_Variation_by_batch_normalized_filtered.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Group), xlab = "PC1", ylab = "PC2", main="SampleGroup")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

# Check for extreme outliers
o1 <- 3*sd(Prin.comp$PC1)
o2 <- 3*sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 & abs(Prin.comp$PC2) > o2) # comment: non

# Test for batch effects via ANOVA
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC1~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC2~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC3~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC4~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC1~Prin.comp$PMI))
anova(lm(Prin.comp$PC2~Prin.comp$PMI))
anova(lm(Prin.comp$PC3~Prin.comp$PMI))
anova(lm(Prin.comp$PC4~Prin.comp$PMI))
anova(lm(Prin.comp$PC1~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC2~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC3~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC4~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Row)))

# First correction for Array:
mod <- model.matrix(~1, data=pd_qual)
M_combat_1array = ComBat(mval_ordered,batch = pd_qual$Array, mod = mod)
save(M_combat_1array,file="M_combat_1array.Rdata")

## Check if batch effect was removed successfully
PCobj = prcomp(t(M_combat_1array), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs = PCs[,1:R]
Prin.comp <- merge(PCs,pd_qual, by = "row.names",all=T)

pdf(file="PC_Variation_by_batch_afterCombat1array.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Group), xlab = "PC1", ylab = "PC2", main="SampleGroup")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC1~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC2~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC3~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC4~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC1~Prin.comp$PMI))
anova(lm(Prin.comp$PC2~Prin.comp$PMI))
anova(lm(Prin.comp$PC3~Prin.comp$PMI))
anova(lm(Prin.comp$PC4~Prin.comp$PMI))
anova(lm(Prin.comp$PC1~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC2~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC3~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC4~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Row)))

# Second correction for Row
mod <- model.matrix(~1, data=pd_qual)
M_combat_2row = ComBat(M_combat_1array,batch = pd_qual$Row, mod = mod)
save(M_combat_2row,file="M_combat_2row.Rdata")

## Check if batch effect was removed successfully
PCobj = prcomp(t(M_combat_2row), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs = PCs[,1:R]
Prin.comp <- merge(PCs,pd_qual, by = "row.names", all=T)

pdf(file="PC_Variation_by_batch_afterCombat2_row.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Group), xlab = "PC1", ylab = "PC2", main="Sample_Group")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Column), xlab = "PC1", ylab = "PC2", main="Column")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Row), xlab = "PC1", ylab = "PC2", main="Row")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Hemisphere), xlab = "PC1", ylab = "PC2", main="Hemisphere")
dev.off()

anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Hemisphere)))
anova(lm(Prin.comp$PC1~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC2~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC3~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC4~Prin.comp$Brain_ph))
anova(lm(Prin.comp$PC1~Prin.comp$PMI))
anova(lm(Prin.comp$PC2~Prin.comp$PMI))
anova(lm(Prin.comp$PC3~Prin.comp$PMI))
anova(lm(Prin.comp$PC4~Prin.comp$PMI))
anova(lm(Prin.comp$PC1~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC2~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC3~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC4~Prin.comp$Freezer_storage))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Slide)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Array)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Sample_Group)))
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Column)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Column))) 
anova(lm(Prin.comp$PC1~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC2~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC3~as.factor(Prin.comp$Row)))
anova(lm(Prin.comp$PC4~as.factor(Prin.comp$Row)))  

# comment: batch effect of brain_ph and freezer_storage to be included as covariates

# Convert the batch-adjusted normalized and filtered M-values back into betas
expit2 = function(x) 2^x/(1+2^x)
betas_normalized_filtered_combated = expit2(M_combat_2row)

# Plot final densities
pdf("BetaValue_Distributions_afterNormalization_Filtering_Combat.pdf")
densityPlot(betas_normalized_filtered_combated, sampGroups = pd_qual$Slide, legend=FALSE, main = "PostQC - Normalized, Filtered and Batch Corrected Beta", xlab = "Beta")
dev.off()
all.equal(colnames(betas_normalized_filtered_combated),rownames(pd_qual))
betas_normalized_filtered_combated_ExprSet = new("ExpressionSet", exprs= as.matrix(betas_normalized_filtered_combated), phenoData=annotated_pd_qual)
save(betas_normalized_filtered_combated_ExprSet,file="betas_normalized_filtered_combated_ExprSet.Rdata")

# End of Script
