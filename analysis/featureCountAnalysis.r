##########################
### SETUP ################
##########################

setwd("/protected/projects/pulmseq/Genentech/coreData/")

#coreExprs <- read.table("coreExprs.txt")
#coreIpfVsNorm <- read.csv("coreIpfVsNorm.txt", header=T, sep="\t")
corePheno <- read.csv("corePheno.txt", header = T, sep="\t")

fileDir <- "/restricted/projectnb/pulmseq/fastq/K00205_30_H3LW3BBXX"
filenames <- data.frame("filenames"=as.character(list.files(fileDir)))
filenames$LIB <- unlist(lapply(strsplit(as.character(filenames$filenames), "_"), function(x){return(x[2])}))
filenames$sample <- unlist(lapply(strsplit(as.character(filenames$filenames), "_"), function(x){return(x[3])}))
filenames$Run <- unlist(lapply(strsplit(as.character(filenames$filenames), "_"), function(x){return(x[1])}))
filenames$FastQ_files <- file.path(fileDir, filenames$filenames)
filenames$sample_name <- paste(filenames$Run, filenames$sample, sep="_")
filenames <- filenames[, c("FastQ_files", "sample_name", "sample", "Run", "filenames", "LIB")]

# WHY ARE THERE DUPLICATES?!!!!
# The duplicates were re-run because the first run was not sequenced deeply enough
# Genentech has been combining the runs - secondary runs are 3310, 3336
# THERE ARE ONLY 90 SAMPLES!

# for now just remove the second samples - this is not a good final solution...
# appears that I've already dealt with the combining of samples into the VSD data
#filenames <- filenames[!duplicated(filenames$sample), ]
phenoData <- merge(filenames, corePheno, by.x = "sample", by.y="SAMPLE_ID")
phenoData <- phenoData[, c(2, 3, 1, 4:dim(phenoData)[2])]

setwd("/protected/projects/pulmseq/Genentech_IPF_Analysis/")

SAVED <- TRUE
filename <- "2016-11-30IPF_resequencedDataCombined.RData"

require(DESeq2)

if(!SAVED){
  
  featureCnts <- read.csv("results/hydra/deliverables/featureCount_raw_counts.txt", sep="\t")
  #sampleInfo <- read.csv("results/hydra/deliverables/sample_info.txt", sep="\t")
  sampleInfo <- phenoData
  counts2 <- as.matrix(featureCnts)
  row.names(counts2) <- counts2[, 1]
  counts2 <- counts2[, 2:dim(counts2)[2]]
  class(counts2) <- "numeric"
  
  sumCounts <- apply(counts2, 2, sum)
  plot(sumCounts, col=as.numeric(as.factor(sampleInfo$Run)), pch=as.numeric(as.factor(duplicated(sampleInfo$sample)))+14, cex=2)
  
  # remove the features not expressed in any of the population (at counts >= 1)
  # the function here returns TRUE/FALSE for features expressed in at least 1 sample
  keepWhich <- apply(counts2, 1, function(x){sum(x==0)<length(x)})
  counts3 <- counts2[keepWhich, ]
  
  sumCounts <- apply(counts3, 2, sum)
  plot(sumCounts, col=as.numeric(as.factor(sampleInfo$Run)), pch=as.numeric(as.factor(duplicated(sampleInfo$sample)))+14, cex=2)
  
  # Follow from DESeq2 Manual
  
  counts3_combined <- matrix(nrow = nrow(counts3), ncol=nrow(sampleInfo) - sum(duplicated(sampleInfo$sample)))
  colnames(counts3_combined) <- unique(sampleInfo$sample)
  rownames(counts3_combined) <- rownames(counts3)
  # for each duplicated sample, combine the counts from each sequencing run
  j <- 1
  for(i in unique(sampleInfo$sample)){
    #   print(grep(i, colnames(counts3)))
    inds <- grep(i, colnames(counts3))
    if(length(inds) > 1){
      counts3_combined[, j] <- apply(counts3[, inds], 1, sum)
    } else {
      counts3_combined[, j] <- counts3[, inds]
    }
    j <- j + 1
  }
  
  allData <- sampleInfo
  allData$duplicated <- allData$sample %in% allData$sample[duplicated(allData$sample)]
  allData <- allData[!duplicated(allData$sample), ]
  
  sumCounts <- apply(counts3_combined, 2, sum)
  plot(sumCounts, col=as.numeric(as.factor(allData$duplicated)), cex=2, pch=16)
  
  # to account for nested patients go to page 38/39
  allData$NESTED_PATIENT <- numeric(dim(allData)[1])
  for(i in levels(allData$TISSUE_DIAGNOSIS)){
    k <- 1
    for(j in unique(allData$PATIENT[allData$TISSUE_DIAGNOSIS==i])){
      allData$NESTED_PATIENT[allData$PATIENT==j] <- k
      k <- k + 1
    }
  }
  
  allData$NESTED_PATIENT <- as.factor(allData$NESTED_PATIENT)
  rownames(allData) <- allData$sample
  
  # make sure they're in the same order
  counts4 <- counts3_combined[, match(as.character(allData$sample), colnames(counts3_combined))]
  
  # for information about continuous variables go to page 50 (they can be included normally)
  
  dds <- DESeqDataSetFromMatrix(countData = counts4, colData = allData,
                                design = ~ 1)
  dds2 <- DESeq(dds) # this should normalize the data for use with pheatmap below...?
  
  # QUALITY CONTROL
  # DESeq 2.1.2 Extracting transformed values - i.e. for quality assurance/control and visualization
  #rld <- rlog(dds, blind=TRUE) # this takes a long time
  vsd <- varianceStabilizingTransformation(dds2, blind=TRUE)
  #vsd.fast <- vst(dds, blind=TRUE) # not available in version of DESeq I am using...? 
  
  # save data so far 
  save.image(file = paste(Sys.Date(), "IPF_resequencedDataCombined.RData", sep=""))
} else {
  load(filename)
}

# load the new data that was sent by Katrina M
newPheno <- read.csv("20161201_Hogg Lab MicroCT and Histology Data_combined.csv", header=T)
newPheno <- newPheno[1:90, ]
newPhenoOrig <- newPheno

newPheno <- read.csv("20170118_Hogg_Lab_MicroCT_and_Histology_Data_combined.csv", header=T)

remove0s <- function(str){
  temp <- str
  str_length <- nchar(temp)
  last_char <- substr(temp, str_length, str_length)
  if(last_char!="0" & last_char!=" ") {
    return(temp)
  } else{
    temp <- substr(temp, 1, str_length-1)
    return(remove0s(temp))
  }
}

# make an ID on which to merge the old and new data
colData(vsd)$core <- paste(gsub("PATIENT", "", colData(vsd)$PATIENT), apply(matrix(gsub("SAMPLE", "", colData(vsd)$SAMPLE)), 1, remove0s), sep=".")

newPheno$core2 <- as.character(newPheno$core)
newPheno$core2[grep(".", newPheno$core2, fixed = T, invert = T)] <- paste(newPheno$Individual[grep(".", newPheno$core2, fixed = T, invert = T)], newPheno$core2[grep(".", newPheno$core2, fixed = T, invert = T)], sep=".")
newPheno$core2 <- apply(matrix(newPheno$core2), 1, remove0s)
newPheno$core <- newPheno$core2
colData(vsd)$core <- apply(matrix(colData(vsd)$core), 1, remove0s)


# merge the new data with the old using core
tempPheno <- merge(colData(vsd), newPheno, by = "core")
tempPheno <- tempPheno[match(colData(vsd)$core, tempPheno$core), ]
rownames(tempPheno) <- rownames(colData(vsd))
colData(vsd) <- tempPheno
colData(dds2) <- colData(vsd)

pdf(paste(Sys.Date(), "PhenotypeFigures.pdf", sep="_"))
mar.default <- c(5.1, 4.1, 4.1, 2.1)
par(mar=mar.default + c(0, 2, 0, 0))
plot(vsd$AlveolarSurfaceDensityperMicro.y ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Alveolar Surface Density per Micro",
     main="Alveolar Surface Density per Micro by Diagnosis", )

plot(vsd$Percent_VolumeFraction_TotalCollagen_Histology ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="% Vol. Fraction Total Collagen",
     main="Vol. Fraction Total Collagen by Diagnosis", )

hist(vsd$Fibroblastic.Foci, main="Fibroblastic Foci in IPF", xlab="Fibroblastic Foci",
     cex.main=2, cex.lab=1.8, cex.axis=1.6)

hist(-log(vsd$Fibroblastic.Foci), main="Negative Log of Fibroblastic Foci in IPF", xlab="Fibroblastic Foci",
     cex.main=2, cex.lab=1.8, cex.axis=1.6)

hist(vsd$AshcroftFibrosisScore_Histology, main="Ashcroft Fibrosis Score in IPF", xlab="Ashcroft Fibrosis Score",
     cex.main=2, cex.lab=1.8, cex.axis=1.6)

plot(vsd$No.TB.permL_MicroCT ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Number of TB per mL",
     main="Number of TB per mL by Diagnosis")

plot(vsd$Tissue._MicroCT ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Tissue MicroCT",
     main="Tissue MicroCT by Diagnosis")

plot(vsd$Neutrophil ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Neutrophil",
     main="Neutrophil by Diagnosis")

plot(vsd$Macrophage.CD68. ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Macrophage CD68",
     main="Macrophage CD68 by Diagnosis")

plot(vsd$Bcell.CD79a. ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="B Cell CD79A",
     main="B Cell CD79A by Diagnosis")

plot(log(vsd$Bcell.CD79a.) ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Log of B Cell CD79A",
     main="Log of B Cell CD79A by Diagnosis")

plot(vsd$CD4 ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="CD4",
     main="CD4 by Diagnosis")

plot(log(vsd$CD4) ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Log of CD4",
     main="Log of CD4 by Diagnosis")

plot(vsd$CD8 ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="CD8",
     main="CD8 by Diagnosis")

plot(vsd$Eosinophils ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Eosinophils",
     main="Eosinophils by Diagnosis")

plot(vsd$Elastin ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Elastin",
     main="Elastin by Diagnosis")

plot(vsd$Col1 ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="Col1",
     main="Col1 by Diagnosis")

plot(vsd$VvCol3 ~ vsd$TISSUE_DIAGNOSIS, col=c("grey90", "grey20"),
     cex.main=2, cex.lab=1.8, cex.axis=1.6,
     xlab="Diagnosis", ylab="VvCol3",
     main="VvCol3 by Diagnosis")

dev.off()


### Generate a table 1
source("/protected/projects/pulmarray/Allegro/COPD_Cancer/functions/sourceDir.R")
sourceDir("/protected/projects/pulmarray/HoggCOPD/A1ATGrifols/R/tableONE/R/")
table1(colData(vsd), groupChar = "TISSUE_DIAGNOSIS", testGroupDiffs = F, saveResultsCSV = T, filename = "IPF_Table1.csv",
       characteristics = c("AlveolarSurfaceDensityperMicro.y",
                           "Percent_VolumeFraction_TotalCollagen_Histology",
                           "Fibroblastic.Foci",
                           "AshcroftFibrosisScore_Histology", "NESTED_PATIENT"))

table1(colData(vsd), groupChar = "TISSUE_DIAGNOSIS", testGroupDiffs = T, saveResultsCSV = T, filename = "IPF_Table1_pvals.csv",
       characteristics = c("AlveolarSurfaceDensityperMicro.y",
                           "Percent_VolumeFraction_TotalCollagen_Histology",
                           "No.TB.permL_MicroCT",
                           "Tissue._MicroCT",
                           "Neutrophil",
                           "Macrophage.CD68.",
                           "Bcell.CD79a.",
                           "CD4", "CD8", "Eosinophils",
                           "Elastin", "Col1", "VvCol3"))


####################################################################################
####################################################################################
###### Qualtity control 
####################################################################################
####################################################################################
saveImages <- T
if(saveImages)
{
  pdf(paste(Sys.Date(), "IPF_Analysis_QC_.pdf", sep=""))
}
# suggestions from Vinay
# visualize boxplot of all counts for each sample - color by condition or batch
boxplot(assay(vsd), col=vsd$TISSUE_DIAGNOSIS, main="Counts colored by Condition")

# median absolute deviation (MAD) across all genes
# take top 5000 or 10000, heatmap of those, look at clustering - are they batch or by condition
genes_mads <- apply(assay(vsd), 1, mad)
genes_mads <- genes_mads[order(genes_mads, decreasing = TRUE)]
top_genes_mads <- genes_mads[1:5000]


.libPaths(c("/unprotected/projects/cbmhive/R_packages/R-3.1.1", .libPaths()))
library(heatmap3)
bluered <- colorRampPalette(c("blue","white","red"))(256)

clabels <- cbind(as.factor(vsd$Run), as.factor(vsd$TISSUE_DIAGNOSIS))
heatmap3(assay(vsd)[match(names(top_genes_mads), rownames(vsd)), ],
         col=bluered,
         keep.dendro=TRUE,
         #col.clustering="unsupervised",
         hclustfun=function(d) hclust(d, method="ward.D2"),
         ColSideColors=clabels
)

# maybe set cook's cutoff to be FALSE - this will prevent DESEq dependent filtering
# independent filtering also maybe turn off

if(saveImages)
{
  dev.off()
}
####################################################################################
####################################################################################
# DESeq2 Quality Control section 2.2
# 2.2.1 Heatmap of the count matrix
#install.packages("pheatmap")

if(saveImages)
{
  pdf(paste(Sys.Date(), "IPF_Preliminary_Analysis_Results_QC_AlveolarDensity.pdf", sep=""))
}

df <- as.data.frame(colData(dds2)[, c("TISSUE_DIAGNOSIS", "Run", "AlveolarSurfaceDensityperMicro.y")])

library("pheatmap")
select <- order(rowMeans(counts(dds2, normalized=TRUE)), decreasing=TRUE)[1:20]
pheatmap(assay(vsd)[select, ], color = bluered, cluster_rows=T,
         show_rownames=T, cluster_cols=FALSE, annotation_col=df,
         main = "Top 20 Expressed Genes\nUnclustered", scale="row")

pheatmap(assay(vsd)[select, ], color = bluered, cluster_rows=TRUE,
         show_rownames=T, cluster_cols=TRUE, annotation_col=df,
         main = "Top 20 Expressed Genes\nClustered", scale="row")

# 2.2.2 Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$TISSUE_DIAGNOSIS, vsd$Run, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$TISSUE_DIAGNOSIS, vsd$Run, sep="-")

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, col=bluered)

# 2.2.3 Principal component plot of the samples
require(ggplot2)
p <- plotPCA(vsd, intgroup=c("TISSUE_DIAGNOSIS", "Run"))
p + labs(title="PCA of All Samples") + aes(shape=TISSUE_DIAGNOSIS)

####################################################################################
####################################################################################
###### Model Design Section
####################################################################################
####################################################################################
# First need to check on correlations between the variables of interest
varMat <- colData(vsd)[, c("AlveolarSurfaceDensityperMicro.y",
                           "Percent_VolumeFraction_TotalCollagen_Histology",
                           "No.TB.permL_MicroCT", 
                           "Tissue._MicroCT",
                           "Elastin",
                           "Col1",
                           "VvCol3",
                           "Macrophage.CD68.",
                           "Bcell.CD79a.",
                           "CD4", "CD8", "Neutrophil", "Eosinophils"
                           )]



colnames(varMat) <- c("AlveolarDensity", "PercentVvCollagen",
                      "NoTBPerMl", "MicroCT", "Elastin", "Col1", "VvCol3",
                      "MacrophageCD68",
                      "BCellCD79A",
                      "CD4", "CD8", "Neutrophil", "Eosinophils")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

require(corrplot)
varDat <- data.frame(varMat)
corMat <- cor(varDat)
corMat2 <- corr.test(varDat, adjust = "none")

pdf("VariableCorrelationsCorrPlot.pdf")
corrplot(corMat, method = "color", diag = T, p.mat = corMat2$p, insig="blank", order="hclust", hclust.method = "single")
dev.off()

pdf(paste(Sys.Date(), "VariableCorrelations.pdf", sep="_"))
pairs(varMat, upper.panel = panel.cor
dev.off()


# try incorporating the alveolar density
# probably don't want to use a diagnosis*patient interaction here as there actually isn't one
#dds_density <- dds2[, !is.na(dds2$AlveolarSurfaceDensityperMicro.y)]
m1 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT, colData(dds2))
m1 <- m1[ ,-which(apply(unname(m1), 2, sum)==0)]

m2 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + AlveolarSurfaceDensityperMicro.y, colData(dds2))
m2 <- m2[ ,-which(apply(unname(m2), 2, sum)==0)]

m3 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + AlveolarSurfaceDensityperMicro.y + TISSUE_DIAGNOSIS:AlveolarSurfaceDensityperMicro.y, colData(dds2))
m3 <- m3[ ,-which(apply(unname(m3), 2, sum)==0)]

m4 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Elastin, colData(dds2))
m4 <- m4[ ,-which(apply(unname(m4), 2, sum)==0)]

m5 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + VvCol3, colData(dds2))
m5 <- m5[ ,-which(apply(unname(m5), 2, sum)==0)]

# New models
m6 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Percent_VolumeFraction_TotalCollagen_Histology, colData(dds2))
m6 <- m6[ ,-which(apply(unname(m6), 2, sum)==0)]

m7 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + No.TB.permL_MicroCT, colData(dds2))
m7 <- m7[ ,-which(apply(unname(m7), 2, sum)==0)]

m8 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Tissue._MicroCT, colData(dds2))
m8 <- m8[ ,-which(apply(unname(m8), 2, sum)==0)]

m9 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Col1, colData(dds2))
m9 <- m9[ ,-which(apply(unname(m9), 2, sum)==0)]

m10 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Macrophage.CD68., colData(dds2))
m10 <- m10[ ,-which(apply(unname(m10), 2, sum)==0)]

m11 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Bcell.CD79a., colData(dds2))
m11 <- m11[ ,-which(apply(unname(m11), 2, sum)==0)]

m12 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Neutrophil, colData(dds2))
m12 <- m12[ ,-which(apply(unname(m11), 2, sum)==0)]

m13 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + CD4, colData(dds2))
m13 <- m13[ ,-which(apply(unname(m13), 2, sum)==0)]

m14 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + CD8, colData(dds2))
m14 <- m14[ ,-which(apply(unname(m14), 2, sum)==0)]

m15 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + Eosinophils, colData(dds2))
m15 <- m15[ ,-which(apply(unname(m15), 2, sum)==0)]


# for a likelihood ratio test can provide a reduced model to the 'reduced' argument
# likelihood ratio test for surface density
dds3 <- DESeq(dds2, full=m2, reduced=m1, test="LRT")
dds4 <- DESeq(dds2, full=m3, reduced=m2, test="LRT")
dds5 <- DESeq(dds2, full=m4, reduced=m1, test="LRT")
dds6 <- DESeq(dds2, full=m5, reduced=m1, test="LRT")
dds7 <- DESeq(dds2, full=m6, reduced=m1, test="LRT")
dds8 <- DESeq(dds2, full=m7, reduced=m1, test="LRT")
dds9 <- DESeq(dds2, full=m8, reduced=m1, test="LRT")
dds10 <- DESeq(dds2, full=m9, reduced=m1, test="LRT")
dds11 <- DESeq(dds2, full=m10, reduced=m1, test="LRT")
dds12 <- DESeq(dds2, full=m11, reduced=m1, test="LRT")
dds13 <- DESeq(dds2, full=m12, reduced=m1, test="LRT")
dds14 <- DESeq(dds2, full=m13, reduced=m1, test="LRT")
dds15 <- DESeq(dds2, full=m14, reduced=m1, test="LRT")
dds16 <- DESeq(dds2, full=m15, reduced=m1, test="LRT")

save.image(file = paste(Sys.Date(), "ddsResults.rdata"))

res1 <- results(dds3, name="AlveolarSurfaceDensityperMicro.y")
res2 <- results(dds4, name="TISSUE_DIAGNOSISIPF.AlveolarSurfaceDensityperMicro.y")
res3 <- results(dds5, name="Elastin")
res4 <- results(dds6, name="VvCol3")
res5 <- results(dds7, name="Percent_VolumeFraction_TotalCollagen_Histology")
res6 <- results(dds8, name="No.TB.permL_MicroCT")
res7 <- results(dds9, name="Tissue._MicroCT")
res8 <- results(dds10, name="Col1")
res9 <- results(dds11, name="Macrophage.CD68.")
res10 <- results(dds12, name="Bcell.CD79a.")
#res11 <- results(dds13, name="Neutrophil")
res12 <- results(dds14, name="CD4")
res13 <- results(dds15, name="CD8")
res14 <- results(dds16, name="Eosinophils")

save.image(file = paste(Sys.Date(), "ddsResults.rdata"))

deseqSummary <- function(result, p.val=0.05){
  #print(paste("Number of non converging:", summary(is.na(result$padj))[["TRUE"]]))
  res.ens <- rownames(result)[which(result$padj < p.val)]
  res.ens.up <- rownames(result)[which(result$padj < p.val & result$log2FoldChange > 0)]
  res.ens.dn <- rownames(result)[which(result$padj < p.val & result$log2FoldChange < 0)]
  print(paste("Number of significant genes:", length(res.ens)))
  #print(paste("Number at uncorrected p:", length(rownames(result)[which(result$pvalue < p.val)])))
  list("res"=res.ens, "res.ens.up"=res.ens.up, "res.ens.dn"=res.ens.dn)
}

fdr <- .1
res.ens1 <- deseqSummary(res1, 0.05)
deseqSummary(res2, fdr)
deseqSummary(res3, fdr)
deseqSummary(res4, fdr)
deseqSummary(res5, fdr)
deseqSummary(res6, fdr)
deseqSummary(res7, fdr)
deseqSummary(res8, fdr)
res.ens9 <- deseqSummary(res9, fdr)
deseqSummary(res10, fdr)
#deseqSummary(res11, fdr)
deseqSummary(res12, fdr)
deseqSummary(res13, fdr)
deseqSummary(res14, fdr)




fdrVal <- 0.05
res1.ens <- rownames(res1)[which(res1$padj < fdrVal)]
res1.ens.up <- rownames(res1)[which(res1$padj < fdrVal & res1$log2FoldChange > 0)]
res1.ens.dn <- rownames(res1)[which(res1$padj < fdrVal & res1$log2FoldChange < 0)]

if(FALSE){
  write.table(res1.ens, row.names=F, col.names=F, quote=F, file = paste(Sys.Date(), "IPF_AlveolarDensity_all_ensembl_IDS.txt"))
  write.table(res1.ens.dn, row.names=F, col.names=F, quote=F, file = paste(Sys.Date(), "IPF_AlveolarDensity_dn_ensembl_IDS.txt"))
  write.table(res1.ens.up, row.names=F, col.names=F, quote=F, file = paste(Sys.Date(), "IPF_AlveolarDensity_up_ensembl_IDS.txt"))
}

res2.ens <- rownames(res2)[which(res2$padj < .25)]

if(FALSE)
{
  pdf(paste(Sys.Date(), "IPF_Preliminary_Analysis_Results_Heatmaps_AlveolarDensity.pdf", sep=""))
}

# create a new colorbar that works better with pheatmap
sharpBlueRed <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(256)
# create a heatmap of the 'density' genes fdr < 0.05; clustered
pheatmap(assay(vsd)[match(res1.ens, rownames(vsd)), ], color=sharpBlueRed, scale="row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, main = "FDR Genes < 0.1 res1")

# create an ordered version of the data - ordered by density (as in Lm study)
oVsd <- vsd[, order(sampleInfo$AlveolarSurfaceDensityperMicro)]
oSampleInfo <- sampleInfo[order(sampleInfo$AlveolarSurfaceDensityperMicro), ]


# recreate the heatmap, supervised by density
pheatmap(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], color=sharpBlueRed, scale="row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, main = "FDR Genes < 0.1 main effect")
pheatmap(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], color=sharpBlueRed, scale="row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, main = "FDR Genes < 0.1 main effect", filename = paste(Sys.Date(), "_mainEffectfdr05.pdf", sep=""))

mainEffectMap <- heatmap3(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], col=sharpBlueRed, col.clustering = "supervised", keep.dendro = T)
mainClusters <- cutree(as.hclust(mainEffectMap$Rowv), 2)
mainEffectMap <- heatmap3(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], col=sharpBlueRed, col.clustering = "supervised", keep.dendro = T, RowSideColors = mainClusters)

# create a boxplot of the cluster means, colored by group, x-axis is density
# use gsva - option rnaseq=TRUE
#mainClusters

# create a heatmap of the genes for the interaction

pheatmap(assay(oVsd)[match(res2.ens, rownames(oVsd)), ], color=sharpBlueRed, scale="row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, main = "FDR Genes < 0.25 interaction")

results.overlapping <- intersect(res1.ens, res2.ens)

#save.image(paste(Sys.Date(), "_IPFResults.rds", sep=""))

ens2Symbol <- function(ensembl){
  require(org.Hs.eg.db)
  x <- org.Hs.egENSEMBL
  xx <- as.list(org.Hs.egENSEMBL2EG)
  
  # now go entrez to gene symbol
  gene_symbols_db <- org.Hs.egSYMBOL
  mapped_symbols <- mappedkeys(gene_symbols_db)
  mapped_symbols <- as.list(gene_symbols_db[mapped_symbols])
  
  inds <- match(ensembl, names(xx))
  inds <- inds[!is.na(inds)]
  
  entrez <- xx[inds] 
  entrez <- unlist(entrez)

  symbols <- mapped_symbols[match(entrez, names(mapped_symbols))]
  symbols <- unlist(symbols)
  return(as.character(symbols))
  
}


res9.symbols.up <- ens2Symbol(res.ens9$res.ens.up)
res9.symbols.dn <- ens2Symbol(res.ens9$res.ens.dn)

res1.symbols.up <- ens2Symbol(res.ens1$res.ens.up)
res1.symbols.dn <- ens2Symbol(res.ens1$res.ens.dn)
write.table(c("genesUpRes1", "genes increased with density", res1.symbols.up), quote=F, col.names=F, row.names=F, file = "genesUpDensity.gmx")
write.table(c("genesDnRes1", "genes decreased with density", res1.symbols.dn), quote=F, col.names=F, row.names=F, file = "genesDnDensity.gmx")

require(enrichR)
dbs <- c("Human_Gene_Atlas", "Biocarta_2016", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "GO_Biological_Process_2015", "GO_Cellular_Component_2015", "GO_Molecular_Function_2015")
res9.enrich <- enrichFullGeneList(up.genes = res9.symbols.up,
                                  dn.genes = res9.symbols.dn,
                                  databases = dbs)
res9.enrich <- res9.enrich[order(res9.enrich$qval), ]
res9.enrich <- res9.enrich[res9.enrich$qval < 0.25, ]

res1.enrich <- enrichFullGeneList(up.genes = res1.symbols.up,
                                  dn.genes = res1.symbols.dn,
                                  databases = dbs)
res1.enrich <- res1.enrich[order(res1.enrich$qval), ]
res1.enrich <- res1.enrich[res1.enrich$qval < 0.25, ]
View(res1.enrich)


write.table(ipfEnrichr, "ipfEnrichrfdr05_alveolardensity.txt", quote=F, row.names=F, col.names=T)
View(ipfEnrichr)
#require(biomaRt)
#ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

# Can paralellize if multiple cores available (would need to request an interactive session)
#library("BiocParallel")
#register(MulticoreParam(4))
#dds2 <- DESeq(dds, parallel=TRUE)
#res <- results(dds2, parallel=TRUE)

if(saveImages)
{
  dev.off()
}



####################################################################################
####################################################################################
###### Model Design Section
####################################################################################
####################################################################################
# try incorporating the alveolar density
# probably don't want to use a diagnosis*patient interaction here as there actually isn't one

sampleInfo_ash <- sampleInfo[!is.na(sampleInfo$AshcroftFibrosisScore_histology), ]

m1 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + AshcroftFibrosisScore_histology, sampleInfo_ash)
m1 <- m1[ ,-which(apply(unname(m1), 2, sum)==0)]

m3 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT + TISSUE_DIAGNOSIS:AshcroftFibrosisScore_histology, sampleInfo_ash)
m3 <- m3[ ,-which(apply(unname(m3), 2, sum)==0)]

m2 <- model.matrix(~TISSUE_DIAGNOSIS + TISSUE_DIAGNOSIS:NESTED_PATIENT, sampleInfo_ash)
m2 <- m2[ ,-which(apply(unname(m2), 2, sum)==0)]

# reassign the dummy "~ 1" design with one that allows 

# for a likelihood ratio test can provide a reduced model to the 'reduced' argument
# likelihood ratio test for surface density

dds2ash <- dds2[, !is.na(colData(dds2)$AshcroftFibrosisScore_histology)]
dds3 <- DESeq(dds2ash, full=m1, reduced=m2, test="LRT")
# likelihood ratio test for surface density:diagnosis interaction 
#dds4 <- DESeq(dds2ash, full=m3, reduced=m2, test="LRT")

res1 <- results(dds3, name="AshcroftFibrosisScore_histology")
#res2 <- results(dds4, name="TISSUE_DIAGNOSISIPF.AshcroftFibrosisScore_histology")

fdrVal <- 0.05
res1.ens <- rownames(res1)[which(res1$padj < fdrVal)]
res1.ens.up <- rownames(res1)[which(res1$padj < fdrVal & res1$log2FoldChange > 0)]
res1.ens.dn <- rownames(res1)[which(res1$padj < fdrVal & res1$log2FoldChange < 0)]

if(FALSE){
  write.table(res1.ens, row.names=F, col.names=F, quote=F, file = paste(Sys.Date(), "IPF_AshcroftFibrosisScore_histology_all_ensembl_IDS.txt"))
  write.table(res1.ens.dn, row.names=F, col.names=F, quote=F, file = paste(Sys.Date(), "IPF_AshcroftFibrosisScore_histology_dn_ensembl_IDS.txt"))
  write.table(res1.ens.up, row.names=F, col.names=F, quote=F, file = paste(Sys.Date(), "IPF_AshcroftFibrosisScore_histology_up_ensembl_IDS.txt"))
}

#res2.ens <- rownames(res2)[which(res2$padj < .01)]

if(TRUE)
{
  pdf(paste(Sys.Date(), "IPF_Preliminary_Analysis_Results_Heatmaps_Ashcroft.pdf", sep=""))
}

# create a new colorbar that works better with pheatmap
sharpBlueRed <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(256)
# create a heatmap of the 'density' genes fdr < 0.05; clustered
vsd_ash <- vsd[, !is.na(colData(vsd)$AshcroftFibrosisScore_histology)]
pheatmap(assay(vsd_ash)[match(res1.ens, rownames(vsd_ash)), ], color=sharpBlueRed, scale="row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, main = "FDR Genes < 0.05 IPF Ashcroft")

# create an ordered version of the data - ordered by density (as in Lm study)
oVsd <- vsd_ash[, order(sampleInfo_ash$AshcroftFibrosisScore_histology)]
oSampleInfo <- sampleInfo_ash[order(sampleInfo_ash$AshcroftFibrosisScore_histology), ]

df <- as.data.frame(colData(dds2)[, c("TISSUE_DIAGNOSIS", "Run", "AshcroftFibrosisScore_histology")])

# recreate the heatmap, supervised by density
pheatmap(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], color=sharpBlueRed, scale="row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, main = "FDR Genes < 0.05 IPF Aschroft")
pheatmap(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], color=sharpBlueRed, scale="row", cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, main = "FDR Genes < 0.05 IPF Ashcroft", filename = paste(Sys.Date(), "_IPFAshcroftmainEffectfdr05.pdf", sep=""))

mainEffectMap <- heatmap3(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], col=sharpBlueRed, col.clustering = "supervised", keep.dendro = T)
mainClusters <- cutree(as.hclust(mainEffectMap$Rowv), 2)
mainEffectMap <- heatmap3(assay(oVsd)[match(res1.ens, rownames(oVsd)), ], col=sharpBlueRed, col.clustering = "supervised", keep.dendro = T, RowSideColors = mainClusters)

# create a boxplot of the cluster means, colored by group, x-axis is density
# use gsva - option rnaseq=TRUE
#mainClusters


save.image(paste(Sys.Date(), "_IPFResults.rds", sep=""))

require(org.Hs.eg.db)
x <- org.Hs.egENSEMBL
xx <- as.list(org.Hs.egENSEMBL2EG)

#res.ensembl <- row.names(res)[which(res$padj < 0.05)] #782 genes
#tempinds <- match(results.overlapping, names(xx))
tempinds <- match(res1.ens, names(xx))
tempinds <- tempinds[!is.na(tempinds)]
res.entrez <- xx[tempinds] # lose 50 genes where there is no matching ensembl to entrez # 732 matches
res.entrez <- unlist(res.entrez) # 25 seem to have at least 2 multiple matches; # 757 IDs

# now go entrez to gene symbol
gene_symbols_db <- org.Hs.egSYMBOL
mapped_symbols <- mappedkeys(gene_symbols_db)
mapped_symbols <- as.list(gene_symbols_db[mapped_symbols])

res.symbol <- mapped_symbols[match(res.entrez, names(mapped_symbols))] # there are 757 gene symbols... bad news bears
res.symbol <- unlist(res.symbol)

if(FALSE){
  write.table(res.symbol, row.names=F, col.names=F, quote=F, file = paste(Sys.Date(), "IPF_AshcroftGenesfdr05_matched_geneSymbol.txt"))
  
}

require(enrichR)
dbs <- c("Biocarta_2016", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "GO_Biological_Process_2015", "GO_Cellular_Component_2015", "GO_Molecular_Function_2015")
ipfEnrichr <- enrichGeneList(res.symbol, databases = dbs)
ipfEnrichr <- ipfEnrichr[order(ipfEnrichr$qval), ]
ipfEnrichr <- ipfEnrichr[ipfEnrichr$qval < 0.25, ]

write.table(ipfEnrichr, "ipfEnrichrfdr25_AshcroftScore_modelfdr05.txt", quote=F, row.names=F, col.names=T)
View(ipfEnrichr)
#require(biomaRt)
#ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

# Can paralellize if multiple cores available (would need to request an interactive session)
#library("BiocParallel")
#register(MulticoreParam(4))
#dds2 <- DESeq(dds, parallel=TRUE)
#res <- results(dds2, parallel=TRUE)

if(saveImages)
{
  dev.off()
}