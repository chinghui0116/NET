## clean working space
rm(list = ls())
## working dir
setwd("E:/microarray/onlinecode")

library("DESeq2")
library("ggplot2")

#########
# NEUseq
#########
## load gene count file
genecount <- read.csv("gene_count_matrix.csv", row.names="gene_id")

## load gene type file
type <- read.table("./mart_export.txt",header=T,sep="\t")

## separate "A|B" to "A" and "B"
genecount$geneid <- row.names(genecount)
as.character(genecount$geneid)
a <-strsplit(genecount$geneid, "|", fixed = TRUE)
#head(a)
#a[[1]][2]
a2 <- lapply(a, "[", 2)
genecount$id <- lapply(a, "[", 2)
genecount <- genecount[order(as.character(genecount$id)),]
genecount = genecount[!is.na(genecount[,45]),]  #delete: id=NA  
genecount$id <- unlist(genecount$id, use.names = FALSE)

id <- data.frame(genecount[,c("id")])
colnames(id)<- c("id")

## calculate average
all <- genecount[,1:43]
all$mean <- apply(all, 1, mean) #1:by row 

## combine gene_symbol
all2 <- cbind(all,id)

## delete low count 
all3 <- subset(all2, abs(mean)>0)## delete mean count 

## delete duplication, keep higher mean gene
all4 <- all3[order(all3$id, all3$mean,decreasing=TRUE),]
all5 <- all4[!duplicated(all4$id), ]
row.names(all5)<- all5$id
#all5 <- all5[,-c(49,50)] #mean, id

## merge gene type file
all6 <- merge(all5,type,by.x="id",by.y="Gene.name",all.x=T)
all7 <- subset(all6, Gene.type=="protein_coding")
protein_anno<- all7[,c("Gene.stable.ID", "Gene.type", "id")]
row.names(all7)<- all7$Gene.stable.ID

all8 <- all7[order(all7$id, all7$mean,decreasing=TRUE),]
all9 <- all8[!duplicated(all8$id), ]
row.names(all9)<- all9$id
all9 <- all9[,-c(1,45,46,47)]  #delete

############
## import NE, NET ratio file
#####
library(xlsx)
nenet <- read.xlsx("./NEUseq/phenotype.xlsx", sheetName = "import")

###### 
#     DEG analysis  - b_NETratio 
#####
## load phenotype data (delete )
delete <- c("UM05","UM07","UM26","UM33","UM36","UM42","UM48","UM57") 
i = -which(nenet$id %in% delete)
del_pheno <- nenet[i,]
condition <- (del_pheno$b_NETratio) 

## load gene file  
all10 <- all9[,c("UM01","UM06","UM09","UM10","UM13","UM15",
                 "UM16","UM17","UM18","UM19","UM20","UM21","UM22","UM23",
                 "UM24","UM25","UM27","UM28","UM29","UM30","UM31",
                 "UM32","UM35","UM37","UM39","UM40","UM41",
                 "UM43","UM44","UM45","UM46","UM47","UM49","UM50",
                 "UM59")]  ## delete missing
genecount <- as.matrix(all10)

colData <- data.frame(row.names = colnames(genecount), condition)
countData <- genecount

# Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))
#  TRUE
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
# > TRUE

# Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
#head(dds)

# filter low: at least 21 samples more than 5  counts
idx <- rowSums(counts(dds, normalized=F) >= 1 ) >= 18
dds <- dds[idx,]

# Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds) #normalize
res <- results(dds)
# Sort by adjusted p-value and display
resordered <- res[order(res$padj),]
summary(res)
b_NETratio <- as.data.frame(resordered)

# filter
b_NETratio$geneSymbol <- row.names(b_NETratio)
net005_fc585 <- subset(b_NETratio, b_NETratio$pvalue<0.05 & abs(b_NETratio$log2FoldChange)>=0.585)
#write.csv(net_005_fc585, file="./net_nolowqc005_fc585.csv",row.names=T)


#####
## NEUseq - ICS response
#####
#########
delete <- c("UM42") 

i = -which(nenet$id %in% delete)
del_pheno <- nenet[i,]

all10 <- all9[,c("UM01","UM05","UM06","UM07","UM09","UM10","UM13","UM15",
                 "UM16","UM17","UM18","UM19","UM20","UM21","UM22","UM23",
                 "UM24","UM25","UM26","UM27","UM28","UM29","UM30","UM31",
                 "UM32","UM33","UM35","UM36","UM37","UM39","UM40","UM41",
                 "UM43","UM44","UM45","UM46","UM47","UM48","UM49","UM50",
                 "UM57","UM59")]  ## delete missing

condition <- as.factor(del_pheno$F1P)
#table(condition)

## load gene file
genecount <- as.matrix(all10)

colData <- data.frame(row.names = colnames(genecount), condition)
countData <- genecount

# Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))
#  TRUE
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
# > TRUE

## chose reference group  (x: categories)
colData$condition <- relevel(colData$condition, ref = "R")

# Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)


# Run the default analysis for DESeq2 and generate results table
#DESeq(): Differential expression analysis based on the Negative Binomial (a.k.a.Gamma-Poisson) distribution
dds <- DESeq(dds) #normalize
res <- results(dds)
# Sort by adjusted p-value and display
resordered <- res[order(res$padj),]
summary(res)


#filter 94 NET related gene
out <- as.data.frame(resordered)
selgene <- net005_fc585$geneSymbol  
i = which(rownames(out) %in% selgene)
out94 <- out[i,] 
#write.csv(out94, file="./NET94_NEUseqR.csv")


library("vsn")
#vsd <- vst(dds, blind = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

### extract normalized values
sampleDists <- assay(vsd)  #R:genename, C:sample, value: normalized gene expression
###
selgene <- net005_fc585$geneSymbol  
i = which(rownames(sampleDists) %in% selgene)
select <- sampleDists[i,] 
#write.csv(select, file="./NET94gene_normalvalue.csv") 


#####
###heatmap 
#####
## delete missing & order by NR/R
select <- select[,c("UM05","UM07","UM10","UM13","UM17","UM18","UM24",
                    "UM27","UM29","UM31","UM32","UM33","UM35","UM37",
                    "UM39","UM41","UM43","UM44","UM48","UM49","UM57","UM59",
                    "UM01","UM06","UM09","UM15","UM16","UM19","UM20",
                    "UM21","UM22","UM23","UM25","UM28","UM30","UM36",
                    "UM40","UM45","UM46","UM47")]
select <- select[,c("UM43","UM49","UM59","UM18","UM31","UM44","UM48",
                    "UM33","UM35","UM17","UM07","UM29","UM27","UM39",
                    "UM13","UM10","UM24","UM32","UM05","UM57","UM41","UM37",
                    "UM30","UM19","UM22","UM06","UM21","UM23","UM09",
                    "UM40","UM16","UM36","UM15","UM46","UM01","UM47",
                    "UM25","UM20","UM45","UM28")]
df <- as.data.frame(colData(dds)[,c("condition")]) #"sizeFactor"
rownames(df) <- del_pheno$id
colnames(df) <- "condition"
library("pheatmap")
#tiff(file="./pheatmap_b_NENETratio_d_F1_good.tiff", width = 10000, height = 6500, bg = "white",
#     res=300, compression = 'lzw')
pheatmap(select, cluster_rows=T, show_rownames=T,
         cluster_cols=F, annotation_col=df)
#dev.off()


### method2 : heatmap
library(tidyverse)
select <- select %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
norm_OEsig <- select[,] %>% 
  filter(select$gene %in% net005_fc585$geneSymbol) %>%  #MUST change filename
  data.frame() %>%
  column_to_rownames(var = "gene") 

###plot
mov10_meta <- data.frame(t(norm_OEsig)) %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()
mov10_meta$sampletype <- as.factor(c(rep("NR",22),rep("R",18)))
annotation <- mov10_meta %>% 
  dplyr::select(samplename, sampletype) %>% 
  data.frame(row.names = "samplename")

#drows = dist(norm_OEsig, method = "minkowski")
# Set a color palette
heat_colors <- colorRampPalette(c("green", "black", "red"))(n = 1000) 
ann_colors = list(
  sampletype = c(NR = "indianred1", R = "cyan3"))
# Run pheatmap
#tiff(file="./pheatmap_F1P.tiff", width = 9000, height = 6500, bg = "white",
#     res=300, compression = 'lzw')
pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_cols = F, 
         cluster_rows = T, 
         clustering_distance_rows = "euclidean", #correlation, euclidean, maximum, manhattan"
         show_rownames = T,
         annotation = annotation, 
         annotation_colors = ann_colors,
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)
#dev.off()


######
##  boxplot - vst - normalized values
#####
### extract normalized values
sampleDists <- assay(vsd)  
###
normalvalue <- net005_fc585$geneSymbol  
i = which(rownames(sampleDists) %in% normalvalue)
normalvalue <- sampleDists[i,] 

gene <- data.frame(t(normalvalue))
gene$id <- row.names(gene)
pheno_94gene <- merge(del_pheno, gene, by="id")

a <-glm(b_NETratio ~ CCL4L2, data=pheno_94gene, 
        family=gaussian(link="identity"))
summary(a)

library(ggplot2)
p<- ggplot(pheno_94gene, aes(x=CCL4L2, y=b_NETratio)) +
  geom_point(size=4)+
  theme_classic()+
  stat_smooth(method='glm',        #lm, glm, gam, loess, rlm
              se=FALSE, alpha=0.3,  # alpha: scale of se if have
              colour="black")+    
  scale_x_continuous("CCL4L2 normalized expression") +
  scale_y_continuous(name="NET ratio", limits = c(0, 7))+
  theme(axis.text.x = element_text(face="bold",size=28,colour="black"),
        axis.text.y = element_text(face="bold",size=28,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ 
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))
p
#ggsave("./CCL4L2_b_NETratio.tiff", width = 9, height = 9.5) #save plot
#dev.off()



#####
# boxplot  - NEUseq
#####
table(pheno_94gene$F1P)

#####
# CCL4L2
#####
f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
p<-ggplot(pheno_94gene, aes(x=F1P, y=CCL4L2)) + 
  geom_boxplot(width=0.1, outlier.shape=NA, lwd=1) +
  stat_summary(fun.data=f, geom="boxplot",size=1.75) + #,width=0.5
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(limits=c("NR","R"), breaks=c("NR", "R"),labels=c("NR(N=22)", "R(N=20)"))+
  theme(axis.text.x = element_text(face="bold",size=14,colour="black"),
        axis.text.y = element_text(face="bold",size=14,colour="black"))
p
#ggsave("./NEUseq_CCL4L2d.tiff", width = 2.5, height = 3.5) #save plot , small plot










############
## Bloodseq
###########
## clean working space
rm(list = ls())

## working dir
setwd("E:/microarray/onlinecode")

library("DESeq2")
library("ggplot2")

## load gene count file
genecount <- read.csv("gene_count_matrix.csv", row.names="gene_id")

## load gene type file
type <- read.table("./mart_export.txt",header=T,sep="\t")

## separate "A|B" to "A" and "B"
genecount$geneid <- row.names(genecount)
as.character(genecount$geneid)
a <-strsplit(genecount$geneid, "|", fixed = TRUE)
#head(a)
#a[[1]][2]
a2 <- lapply(a, "[", 2)
genecount$id <- lapply(a, "[", 2)
genecount <- genecount[order(as.character(genecount$id)),]
genecount = genecount[!is.na(genecount[,48]),]  #delete: id=NA  
genecount$id <- unlist(genecount$id, use.names = FALSE)

id <- data.frame(genecount[,c("id")])
colnames(id)<- c("id")
case <- genecount[,c("QU01","QU08","QU14","QU17","QU21","QU22","QU29",
                     "QU33","QU35","QU37","QU40","QU41","QU48","QU49",
                     "QU61","QU62","QU66")]
control <- genecount[,c("QU03","QU04","QU05","QU10","QU13","QU15",
                        "QU16","QU19","QU20","QU23","QU24","QU25",
                        "QU26","QU30","QU31","QU32","QU34","QU38",
                        "QU39","QU43","QU44","QU45","QU46","QU47",
                        "QU50","QU63","QU64","QU65","QU67")]
## calculate average by case and control
case$mean_ca <- apply(case, 1, mean) 
control$mean_co <- apply(control, 1, mean) 
## combine id, case, control
all <- cbind(id, case, control)
all$meandiff <- (all$mean_ca + all$mean_co)/2
all2=all
all3 <- all2[,c("id","meandiff","QU01","QU03","QU04","QU05","QU08",
                "QU10","QU13","QU14","QU15","QU16","QU17","QU19",
                "QU20","QU21","QU22","QU23","QU24","QU25","QU26","QU29",
                "QU30","QU31","QU32","QU33","QU34","QU35","QU37","QU38","QU39",
                "QU40","QU41","QU43","QU44","QU45","QU46","QU47","QU48","QU49",
                "QU50","QU61","QU62","QU63","QU64","QU65","QU66","QU67")]

## delete duplication, keep higher meandiff gene
all4 <- all3[order(all3$id, all3$meandiff,decreasing=TRUE),]
all5 <- all4[!duplicated(all4$id), ]
row.names(all5)<- all5$id
#all5 <- all5[,-c(2)]

## merge gene type file
all6 <- merge(all5,type,by.x="id",by.y="Gene.name",all.x=T)
all7 <- subset(all6, Gene.type=="protein_coding")
protein_anno<- all7[,c("Gene.stable.ID", "Gene.type", "id")]
row.names(all7)<- all7$Gene.stable.ID

all8 <- all7[order(all7$id, all7$meandiff,decreasing=TRUE),]
all9 <- all8[!duplicated(all8$id), ]
row.names(all9)<- all9$id
all9 <- all9[,-c(1,2,49,50)]

## load gene file
genecount <- as.matrix(all9)

### load phenotype file 
pheno_data = read.csv("./phenodata.csv",stringsAsFactors=FALSE)
condition <- factor(pheno_data$F1P) 
#table(pheno_data$F1P)

colData <- data.frame(row.names = colnames(genecount), condition)
countData <- genecount

# Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))
#  TRUE
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
# > TRUE

## chose reference group  (x: categories)
colData$condition <- relevel(colData$condition, ref = "R")

# Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
#head(dds)

# Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds) #normalize
res <- results(dds, contrast = c("condition","NR","R"))
# Sort by adjusted p-value and display
resordered <- res[order(res$padj),]
summary(res)
#write.csv(as.data.frame(resordered), file="./ICS2_F1P.csv")


##### extract genename, p-values 
ICS_F1P <- as.data.frame(resordered)
ICS_F1P$id <- row.names(ICS_F1P)

#filter 94 NET related gene
net005_fc585<- read.csv("./NET94_NEUseqR.csv", row.names="geneSymbol")

out <- as.data.frame(resordered)
selgene <- row.names(net005_fc585) 
i = which(rownames(out) %in% selgene)
out94 <- out[i,] 
#write.csv(as.data.frame(out94), file="./NET94_BloodseqR.csv")

library("vsn")
#vsd <- vst(dds, blind = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

### extract normalized values
sampleDists <- assay(vsd)  #R:genename, C:sample, value: normalized gene expression
###
i = which(rownames(sampleDists) %in% selgene)
ICS2_94g <- sampleDists[i,] 
#write.csv(ICS2_94g, file="./ICS2_NET94gene_normalvalue.csv") 


#####
## boxplot -CCL4L2 (Bloodseq)
#####
## import normalized gene expression values
ICS2gene <- data.frame(t(ICS2_94g))
ICS2gene$id <- row.names(ICS2gene)

## combine gene and phenotype
ICS2 <- merge(pheno_data, ICS2gene, by.x="ids", by.y="id")

table(ICS2$F1P)

f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
ggplot(ICS2, aes(x=F1P, y= CCL4L2)) +   
  geom_boxplot(width=0.1, outlier.shape=NA, lwd=1) +  
  stat_summary(fun.data=f, geom="boxplot",size=1.75) + 
  theme_classic()+
  labs(title="ICS2",
       x ="ICS response", y = "CCL4L2")+
  scale_x_discrete(limits=c("NR","R"),breaks=c("NR","R"),labels=c("NR(N=22)", "R(N=24)"))+  
  theme(axis.text.x = element_text(face="bold",size=14,colour="black"),
        axis.text.y = element_text(face="bold",size=14,colour="black"))
#ggsave("./Bloodseq_CCL4L2d.tiff", width = 2.5, height = 3.5) #save plot , small plot




