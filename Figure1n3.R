## clean working space
rm(list = ls())

## working dir
setwd("E:/microarray/onlinecode")

## load normalized gene expression profiles
asthma <- load(file="./ACT.rda")
dim(asthma.data)
anno_df <- asthma.data[,100]
asthma <- asthma.data[,-c(100)]

## phenotypes
pd <- read.csv("./act_phenotype.csv", sep=",")

## difeine group
groups <- as.factor(pd$act)
table(groups)



####################################
## PCA  analysis 
####################################
 exprdata <- asthma
## case/control order?
## PCA 
x <- t(exprdata)
pc <- prcomp(x)

## check Chip_ID 
pd$Chip_ID <- paste("X", pd$Chip_ID, sep="")
all(rownames(pd$Chip_ID) == colnames(exprdata))
 ###TRUE

## PCA calculate
pca_data=prcomp(t(exprdata))

## Do some calculation (% variance covered by component to show on X and Y axis for PC1 and PC2)
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)

## Create dataframe with PC1 and PC2 with some metadata:
df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], 
                       sample = colnames(exprdata), condition=groups)
library("ggplot2")
## Plot using ggplot2 (color by case/control):
ggplot(df_pca_data, aes(PC1,PC2, color = condition))+
  ggtitle("PCA")+
  geom_point(size=2.5, shape=1, stroke=2)+
  theme_classic()+
  theme(axis.text.x = element_text(size=18,colour="black"),
        axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title.x =element_text(size=20,colour="black"))+ 
  theme(axis.title.y =element_text(size=20,colour="black"))+ 
  labs(x=paste0("PC1 (",pca_data_perc[1],"%)"), 
       y=paste0("PC2 (",pca_data_perc[2],"%)"))+
  scale_colour_manual(name="condition",values=c("cyan3", "indianred1"))
#table(pd$act)



#######################################
# Differentially expressed gene analysis
#######################################
## design matrix
asthma.design <- model.matrix(~0+groups)
colnames(asthma.design) = c("Control", "Case")

## lmfit
library(limma)
library(lumi)
asthma.fit <- lmFit(asthma, asthma.design)

## contrast matrix
asthma.cont.matrix <- makeContrasts(Case-Control, levels=asthma.design)

## contrast fit
asthma.fit2 <- contrasts.fit(asthma.fit, asthma.cont.matrix)

## ebayes
asthma.ebfit <- eBayes(asthma.fit2)

## annotation
asthma.ebfit$genes <- anno_df

## result output
asthma.DF <- topTable(asthma.ebfit, number=15828, adjust.method="BH")
#write.csv(asthma.DF, file = "./results.csv")

## filter FDR<0.05
asthma_fdr <- subset(asthma.DF, asthma.DF$adj.P.Val<0.05)
dim(asthma_fdr)
 ###[1] 298   6  (298 DEGs)


####################################
## Heatmap
####################################
actfdr005 <- subset(asthma.DF, asthma.DF$adj.P.Val <0.05)

## find a specific gene list 
selGene <- actfdr005$probeList    
i = which(row.names(asthma) %in% selGene)
asthma_sel <- asthma[i,]

## find geneSymbol
i = which(anno_df$probeList %in% row.names(asthma_sel))
anno_df_sel <- anno_df[i,]
geneSymbol_sel <- anno_df_sel$geneSymbol

## chenck name of probeList and geneSymbol
all(rownames(asthma_sel) == rownames(anno_df_sel))
 ###TRUE
row.names(asthma_sel) <- geneSymbol_sel

act <- as.factor(pd$act)
## get heatmap plot  
library("gplots")
color.map <- function(act) { 
  if (act=="U") "indianred1" else "cyan3" }
patientcolors <- unlist(lapply(act, color.map))
tiff(file="./heatmap_gene.tiff", width = 10000, height = 6500, bg = "white",
     res=300, compression = 'lzw')
heatmap.2(as.matrix(asthma_sel), #matrix:row=gene_name; col=sample_id
          col=redgreen(75), 
          ColSideColors=patientcolors, 
          Rowv=TRUE, Colv=T, dendrogram="column", #c("both","row","column","none")
          key = TRUE, keysize = 0.5, 
          #denscol=tracecol, 
          density.info=c("density"),
          trace="none", scale="row",
          cexRow=0.9, cexCol = 0.9,
          main="heatmap")
dev.off()



##### NEU, EOS, LYM, PLT variables (heatmap)
### sample clustering by gene expression vs. subphenotypes(ex. neutrophil,lymphocyte)
var <- read.csv("./heatmap_var.csv", row.names = "id")
group <- var$ACT

### PLT_baseline
var <- data.frame(t(var))
xx  = c("PLT_baseline","Lym_baseline")
var2 <- var[xx, ]

## get heatmap plot  
library("gplots")
color.map <- function(group) { 
  if (group=="1") "indianred1" else "cyan3" }
patientcolors <- unlist(lapply(group, color.map))

tiff(file="./heatmap_PLT.tiff", width = 10000, height = 6500, bg = "white",
     res=300, compression = 'lzw')
library("RColorBrewer")
library("circlize")
cols <- colorRampPalette(c("white", "green4"))(n = 49) 
col_breaks = c(seq(151,300,length=25),  
               seq(301,495,length=25))
heatmap.2(as.matrix(var2), #matrix:row=gene_name; col=sample_id
          col=cols, 
          ColSideColors=patientcolors, 
          Rowv=F, Colv=F, #rerange rank by cluster¡Aneed same with dendrogram=
          dendrogram="none", #c("both","row","column","none")
          key = TRUE, keysize = 0.5, key.par = list(cex=1.5),
          #denscol=tracecol,
          breaks=col_breaks,
          density.info=c("none"),#c("histogram","density","none")
          trace="none", 
          scale="none", #normalization by:"both","row","column","none"
          cexRow=0.9, cexCol = 0.9,
          main="heatmap")
dev.off()


### LYMinc_baseline
xx  = c("LYMinc_baseline","Lym_baseline")
var2 <- var[xx, ]
## get heatmap plot  
library("gplots")
color.map <- function(group) { 
  if (group=="1") "indianred1" else "cyan3" }
patientcolors <- unlist(lapply(group, color.map))
tiff(file="./heatmap_LYMinc.tiff", width = 10000, height = 6500, bg = "white",
     res=300, compression = 'lzw')
library("RColorBrewer")
library("circlize")
cols <- colorRampPalette(c("white", "navy"))(n = 49) 
col_breaks = c(seq(930,2900,length=25),  
               seq(2901,6600,length=25))
heatmap.2(as.matrix(var2), #matrix:row=gene_name; col=sample_id
          col=cols, 
          ColSideColors=patientcolors, 
          Rowv=F, Colv=F, #rerange rank by cluster¡Aneed same with dendrogram=
          dendrogram="none", #c("both","row","column","none")
          key = TRUE, keysize = 0.5, key.par = list(cex=1.5),
          #denscol=tracecol,
          breaks=col_breaks,
          density.info=c("none"),#c("histogram","density","none")
          trace="none", 
          scale="none", #normalization by:"both","row","column","none"
          cexRow=0.9, cexCol = 0.9,
          main="heatmap")
dev.off()


### NEUinc_baseline
xx  = c("NEUinc_baseline","Lym_baseline")
var2 <- var[xx, ]
## get heatmap plot  
library("gplots")
color.map <- function(group) { 
  if (group=="1") "indianred1" else "cyan3" }
patientcolors <- unlist(lapply(group, color.map))
tiff(file="./heatmap_NEUinc.tiff", width = 10000, height = 6500, bg = "white",
     res=300, compression = 'lzw')
library("RColorBrewer")
library("circlize")
cols <- colorRampPalette(c("white", "red3"))(n = 49) 
col_breaks = c(seq(1300,9000,length=25),  
               seq(9001,19200,length=25))
heatmap.2(as.matrix(var2), #matrix:row=gene_name; col=sample_id
          col=cols, 
          ColSideColors=patientcolors, 
          Rowv=F, Colv=F, #rerange rank by cluster¡Aneed same with dendrogram=
          dendrogram="none", #c("both","row","column","none")
          key = TRUE, keysize = 0.5, key.par = list(cex=1.5),
          #denscol=tracecol,
          breaks=col_breaks,
          density.info=c("none"),#c("histogram","density","none")
          trace="none", 
          scale="none", #normalization by:"both","row","column","none"
          cexRow=0.9, cexCol = 0.9,
          main="heatmap")
dev.off()



### EOSinc_baseline
xx  = c("EOSinc_baseline","Lym_baseline")
var2 <- var[xx, ]
## get heatmap plot  
library("gplots")
color.map <- function(group) { 
  if (group=="1") "indianred1" else "cyan3" }
patientcolors <- unlist(lapply(group, color.map))
tiff(file="./heatmap_EOSinc.tiff", width = 10000, height = 6500, bg = "white",
     res=300, compression = 'lzw')
library("RColorBrewer")
library("circlize")
cols <- colorRampPalette(c("white", "orange"))(n = 49) 
col_breaks = c(seq(9,860,length=25),  
               seq(861,1729,length=25))
heatmap.2(as.matrix(var2), #matrix:row=gene_name; col=sample_id
          col=cols, 
          ColSideColors=patientcolors, 
          Rowv=F, Colv=F, #rerange rank by cluster¡Aneed same with dendrogram=
          dendrogram="none", #c("both","row","column","none")
          key = TRUE, keysize = 0.5, key.par = list(cex=1.5),
          #denscol=tracecol,
          breaks=col_breaks,
          density.info=c("none"),#c("histogram","density","none")
          trace="none", 
          scale="none", #normalization by:"both","row","column","none"
          cexRow=0.9, cexCol = 0.9,
          main="heatmap")
dev.off()



###############
# correlation -- NEU,EOS,LYM,PLT vs genes
###############
## select TH1, Th2, Th17 genes
a <- c("fSpP_tQvufNUyUoB3k","odF3XHR8CVl4SAUaUQ","fkSS9RKRCSeNY9eAac",
       "NSoMpEeiLQHuQqrHug","95LUl4qd9pA7JOpF0g","foL5LieM05V7HUl6uE",
       "Nic3g59cAJitXUBb7Q","07riWI_RT5ncl6kEEk","cSj1OhRehCQFQHWCH4",
       "EP5e.efnftcn6in1XQ","Hc1Lx5a.gDOXkN4lLk","lW4AaJ4EMp71TslV6U",
       "ot9VR9V1CFe70.rvXo","c1C92P0Tl9CHgRYcK8","Hfd3h3Yfu9VnXV4X10")
b <- c("IL1A","IL1B","IL4","IL5","IL6","IL10","IL12A","IL18","CCL2",
       "CCL5","IFNG","TNF","TGFB3","IL17A","TSLP")
genelist <- data.frame(a,b)
names(genelist) <- c("probeList","geneSymbol")

i = which(rownames(asthma) %in% genelist$probeList)
select <- asthma[i,] 
select2 <- merge(select, genelist, by.x="row.names", by.y="probeList")
row.names(select2) <- select2$geneSymbol 
select2 <- select2[,-c(1,101)]

t_select2 <- data.frame(t(select2))
t_select2$Chip_ID <- row.names(t_select2) 

## Neutrophils, Lymphocytes, Platelets 
var <- read.csv("./heatmap_var.csv", row.names = "id")
var$Chip_ID <- paste("X", row.names(var), sep="")

mydata <- merge(t_select2, var, by="Chip_ID")

mydata <- mydata[,c("IFNG","IL1A","IL1B","IL4","IL5","IL6","IL10",
                    "IL12A","IL17A","IL18","CCL2","CCL5","TGFB3","TNF","TSLP",
                    "NEUinc_baseline","LYMinc_baseline","PLT_baseline",
                    "EOSinc_baseline")]

# delete missing
mydata <- mydata[complete.cases(mydata), ]

# correlation
mydata.cor = cor(mydata, method = c("spearman"))

library("Hmisc")
mydata.rcorr = rcorr(as.matrix(mydata))
mydata.rcorr

mydata.coeff = mydata.rcorr$r
mydata.p = mydata.rcorr$P

library("corrplot")
corrplot(mydata.cor)

#palette = colorRampPalette(c("blue", "white", "red")) (20)
#heatmap(x = mydata.cor, col = palette, symm = TRUE)

##
part <- cor(mydata)[c("NEUinc_baseline","LYMinc_baseline","PLT_baseline","EOSinc_baseline"),
                    c("IFNG","IL1A","IL1B","IL4","IL5","IL6","IL10",
                      "IL12A","IL17A","IL18","CCL2","CCL5","TGFB3","TNF","TSLP")]
row.names(part) <- c("Neutrophils", "Lymphocytes", "Platelets","Eosinophils")
#corrplot(part, method ="circle", col = palette, tl.cex=1)

####
#install.packages("ggcorrplot")
library("ggcorrplot")
ggcorrplot(part, method="circle", 
           ggtheme = ggplot2::theme_classic, 
           tl.cex=12, tl.col="black",
           lab_col = "black", lab_size = 4)

dev.off()



##########
## WGCNA 
##########
row.names(asthma) <- anno_df$geneSymbol

## geneData: R=sample, C=gene; phenoData: R=sample, C=phenotype(trait)
asthma[1:4,1:4]
t(asthma)[1:4,1:4]  #need to transform for WGCNA 

# find a specific gene list  (!!! MUST check selgene number=2612!!!)
asthma_p005 <- subset(asthma.DF, asthma.DF$P.Value<0.05)
sigGene <- asthma_p005$geneSymbol    
i = which(row.names(asthma) %in% sigGene)
ACTdata_slegene <- asthma[i,]



## separate cases and controls
casedata <- ACTdata_slegene[, c("X200757700032_A",
                                "X3999927049_G",
                                "X3999538003_F",
                                "X3999890003_D",
                                "X3999538003_K",
                                "X9999593137_E",
                                "X3999890003_L",
                                "X200162410130_B",
                                "X3999927049_H",
                                "X200162410130_C",
                                "X3999545107_A",
                                "X200162410130_D",
                                "X200769980006_L",
                                "X200162410130_E",
                                "X3999927049_J",
                                "X200769890047_B",
                                "X3999545107_B",
                                "X3999927049_B",
                                "X3999545107_C",
                                "X9999593126_E",
                                "X200162410130_I",
                                "X9999593126_L",
                                "X200569440030_A",
                                "X200757700032_H",
                                "X200569440030_H",
                                "X200757700032_J",
                                "X200769890047_D",
                                "X200769890047_F",
                                "X200769890056_A",
                                "X200769890056_C",
                                "X200769890056_D",
                                "X200769890056_E",
                                "X200769890056_J",
                                "X3999927051_K",
                                "X3999927051_F",
                                "X3999927051_L",
                                "X3999927051_G",
                                "X200569440030_L",
                                "X3999927049_C",
                                "X9999593126_K",
                                "X9999593126_J",
                                "X200536850007_L",
                                "X200536850032_B")]

library(WGCNA)

### STEP 1. Choose a set of soft-thresholding powers
powers<-c(1:10,seq(12,20,2))
sft<-pickSoftThreshold(t(casedata),powerVector=powers, verbose=5)

# plot the results
par(mfrow=c(1,2)) #mfrow=c(nrows, ncols):combining two plots 
cex1 = 0.9
# scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n", main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
# this line corresponds to using a R^2 cutoff of h
abline(h=0.85,col="red")

# mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
     type="n", main=paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")

# choose beta value (power value)
best_beta=sft$powerEstimate
best_beta
#[1] 7

cor <- WGCNA::cor

### STEP 2:modules identification
nGenes=ncol(t(casedata))
power=best_beta
type = "signed" #unsigned, signed
corType = "pearson" 
maxPOutliers = ifelse(corType=="pearson",1,0.05) 
net = blockwiseModules(t(casedata), power = power,
                       maxBlockSize = nGenes, TOMType = type, 
                       minModuleSize = 100, reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,corType = corType,
                       maxPOutliers=maxPOutliers,
                       saveTOMFileBase = "AS-yellow-FPKM-TOM",
                       verbose = 3)
table(net$colors)

cor<-stats::cor
#STEP 2:modules visualization
# Convert labels to colors for plotting
moduleColors = labels2colors(net$colors)
table(moduleColors)
moduleLabelsAutomatic <- net$colors

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",main = "Gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# A data frame with module eigengenes can be obtained as follows
MEsAutomatic <- net$MEs
#write.csv(MEsAutomatic, file = "./MEsAutomatic.csv")

## extrat genes by color
b<- as.data.frame(moduleLabelsAutomatic)
colnames(b)<- c("colornum")
b$gene <- row.names(b)
browngene <-subset(b, colornum==3)
#write.csv(browngene, "./case_browngene.csv")


######
## drug categories 
## import phenotype data
pd2 <- read.csv("./druguse.csv", sep=",")
pd2$Chip_ID <- paste("X", pd2$Chip_ID, sep="")
row.names(pd2) <- pd2$Chip_ID
pd2 <- pd2[,-c(1)]

all(colnames(casedata) %in% row.names(pd2))
 # TRUE

## check No. of genes and samples
nGenes <- ncol(t(casedata))
nSamples <- nrow(t(casedata))


### STEP5: associations between modules and traits
dataExpr=t(casedata)
traitData= pd2[,c(1:4)] 
######################### 
sampleName = rownames(dataExpr)
traitData = traitData[match(sampleName, rownames(traitData)), ]

## Sample dendrogram and trait heatmap
datExpr_tree<-hclust(dist(t(casedata)), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData, signed = FALSE);

## Eigengene adjacency heatmap
MEs_col= moduleEigengenes(t(casedata), moduleColors)$eigengenes
ME_pheno <- cbind(MEs_col, traitData)
MEs_colpheno = orderMEs(ME_pheno)


###  Module-trait relationships
modTraitCor = cor(MEs_col, traitData, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

newmatrix=-log10(modTraitP)*modTraitCor

dev.off()
#sizeGrWindow(12,6)
#par(mar = c(3, 6, 1, 1))
labeledHeatmap(Matrix = newmatrix, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               #textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 1, zlim = c(-0.4,0.4),
               horizontalSeparator.lty=1,
               main = paste("Module-trait relationships"))




#############
## Controls
#############
controldata <-ACTdata_slegene[, c("X3999927049_E",
                                  "X200536850027_F",
                                  "X3999890003_G",
                                  "X200757700032_B",
                                  "X3999927049_L",
                                  "X200536850027_J",
                                  "X200536850027_K",
                                  "X3999538003_D",
                                  "X3999890003_E",
                                  "X3999538003_E",
                                  "X9999593137_J",
                                  "X200544710005_B",
                                  "X200757700032_C",
                                  "X200544710005_C",
                                  "X3999538003_G",
                                  "X9999593126_G",
                                  "X200757700032_D",
                                  "X3999927049_A",
                                  "X9999593137_G",
                                  "X200544710005_D",
                                  "X200544710005_E",
                                  "X200544710005_J",
                                  "X3999538003_I",
                                  "X200544710005_K",
                                  "X200544710005_L",
                                  "X9999593126_H",
                                  "X200162410130_A",
                                  "X200757700032_F",
                                  "X3999890003_F",
                                  "X3999927049_D",
                                  "X3999890003_H",
                                  "X9999593126_A",
                                  "X200162410130_H",
                                  "X3999545107_D",
                                  "X9999593126_D",
                                  "X3999545107_E",
                                  "X9999593126_B",
                                  "X200162410130_J",
                                  "X200162410130_K", 
                                  "X200569440030_B",
                                  "X200569440030_C",
                                  "X200569440030_D",
                                  "X200757700032_I",
                                  "X200569440030_E",
                                  "X200569440030_F",
                                  "X200569440030_G",
                                  "X200569440030_I",
                                  "X3999927051_H",
                                  "X3999545107_G",
                                  "X200569440030_K",
                                  "X3999545107_H",
                                  "X3999545107_I",
                                  "X200536850007_A",
                                  "X3999545107_L",
                                  "X200536850007_I",
                                  "X200536850032_C")]


library(WGCNA)
### STEP 1. Choose a set of soft-thresholding powers
powers<-c(1:10,seq(12,20,2))
sft<-pickSoftThreshold(t(controldata),powerVector=powers, verbose=5)

# plot the results
par(mfrow=c(1,2)) #mfrow=c(nrows, ncols):combining two plots 
cex1 = 0.9
# scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n", main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
# this line corresponds to using a R^2 cutoff of h
abline(h=0.85,col="red")

# mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
     type="n", main=paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")

# choose beta value (power value)
best_beta=sft$powerEstimate
best_beta
#[1] 7

cor <- WGCNA::cor

### STEP 2:modules identification
nGenes=ncol(t(controldata))
power=best_beta
type = "signed" #unsigned, signed
corType = "pearson" 
maxPOutliers = ifelse(corType=="pearson",1,0.05) 
net = blockwiseModules(t(controldata), power = power,
                       maxBlockSize = nGenes, TOMType = type, 
                       minModuleSize = 100, reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,corType = corType,
                       maxPOutliers=maxPOutliers,
                       saveTOMFileBase = "AS-yellow-FPKM-TOM",
                       verbose = 3)
table(net$colors)

cor<-stats::cor
#STEP 2:modules visualization
# Convert labels to colors for plotting
moduleColors = labels2colors(net$colors)
table(moduleColors)
moduleLabelsAutomatic <- net$colors

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",main = "Gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)




