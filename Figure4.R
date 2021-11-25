## clean working space
rm(list = ls())

## working dir
setwd ("E:/")

## import data
plasma <- read.csv("./plasma.csv", sep=",")

####  NETratio
library("ggsignif")
library("ggplot2")
mean_se(plasma$b_NETratio)
## baseline NET  
ggplot(plasma, aes(x=response, y=b_NETratio, shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.1)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "NET ratio")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./NETratio.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()



####  SIRPA
library("ggsignif")
library("ggplot2")
mean_se(plasma$SIRPA)
ggplot(plasma, aes(x=response, y=log10(SIRPA), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.015)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(SIRPA, ng/mL)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./SIRPA.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()



####  S100A12
library("ggsignif")
library("ggplot2")
mean_se(plasma$S100A12)
ggplot(plasma, aes(x=response, y=log10(S100A12), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.02)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(S100A12, ng/mL)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./S100A12.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()


####  S100A8
library("ggsignif")
library("ggplot2")
mean_se(plasma$S100A8)
ggplot(plasma, aes(x=response, y=log10(S100A8), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.015)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(S100A8, ng/mL)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./S100A8.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()




#####  MIF
library("ggsignif")
library("ggplot2")
mean_se(plasma$MIF)
ggplot(plasma, aes(x=response, y=log10(MIF), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.05)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(MIF, ng/mL)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./MIF.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()



#####  S100A9
library("ggsignif")
library("ggplot2")
mean_se(plasma$S100A9)
ggplot(plasma, aes(x=response, y=log10(S100A9), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.025)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(S100A9, ng/mL)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./S100A9.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()


##### TIMP2
library("ggsignif")
library("ggplot2")
mean_se(plasma$TIMP2)
ggplot(plasma, aes(x=response, y=TIMP2, shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.013)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "TIMP2, ng/mL")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./TIMP2.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()


#####  PECAM1
library("ggsignif")
library("ggplot2")
mean_se(plasma$PECAM1)
ggplot(plasma, aes(x=response, y=log10(PECAM1), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.025)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(PECAM1, ng/mL)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./PECAM1.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()


##### NEUinc_baseline
library("ggsignif")
library("ggplot2")
mean_se(plasma$NEUinc_baseline)
ggplot(plasma, aes(x=response, y=log10(NEUinc_baseline), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.021)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR", "R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(Neutrophil count, #/uL)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./NEUinc.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()


##### neu_baseline
library("ggsignif")
library("ggplot2")
mean_se(plasma$neu_baseline)
ggplot(plasma, aes(x=response, y=log10(neu_baseline), shape=response)) + 
  #Creating your own function within geom_errorbar()
  geom_errorbar(stat="summary",  width=0.3,
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))})+
  geom_dotplot(binaxis='y', stackdir='center', binwidth=0.013)+
  stat_summary(fun=mean, fun.min=mean, fun.max=mean, 
               geom="crossbar", width=0.3)+
  theme_classic()+  #theme_classic : theme with axis lines and no grid lines
  scale_x_discrete(breaks=c("NR","R"),labels=c("NR","R"))+
  theme(axis.text.x = element_text(face="bold",size=32,colour="black"),
        axis.text.y = element_text(face="bold",size=24,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=32,colour="black"))+ #X
  theme(axis.title.y =element_text(face="bold",size=32,colour="black"))+ #Y
  theme(legend.position='none')+
  labs(title="", x ="", y = "log(Neutrophil %)")+
  #wilcox.test (p=), t.test (p=) 
  geom_signif(comparisons = list(c("NR", "R")), test= t.test, 
              textsize = 8, tip_length=0,
              map_signif_level=function(p)sprintf("p = %.2g", p)) 
ggsave("./neu.tiff", 
       width = 5, height = 7.5) #save plot
dev.off()




#############################
##   ROC
############################
#####
nenetROC=plasma
nenetROC$response2[which(nenetROC$response =="R")] <- 0
nenetROC$response2[which(nenetROC$response =="NR")] <- 1


#####
# NET unadjusted and adjusted ROC
##### 
##unadj
unadjglm <- glm(formula = response2 ~ b_NETratio , family = "binomial", data = nenetROC)
summary(unadjglm)
resultunadj <- predict(unadjglm, type = "response")

##adj
adj_glm <- glm(formula = response2 ~ age+gender+b_NETratio , family = "binomial", data = nenetROC)
summary(adj_glm)
resultadj <- predict(adj_glm, type = "response")

##combine unadj & adj
library(dplyr)
nenetROC_unadj <- na.omit(select(nenetROC,response2,b_NETratio))
roc_unadj <- data.frame(cbind(resultunadj, nenetROC_unadj$response2))
roc_unadj$group <- "unadj"
colnames(roc_unadj) <- c("result","response","group")

nenetROC_adj <- na.omit(select(nenetROC,response2,b_NETratio))
roc_adj <- data.frame(cbind(resultadj, nenetROC_adj$response2))
roc_adj$group <- "adj"
colnames(roc_adj) <- c("result","response","group")

roc <- rbind(roc_unadj,roc_adj)

### plot ROC curve
library(ggplot2)
library(plotROC)

roc$response[which(roc$response =="2")] <- 0

basicplot <- ggplot(roc, aes(d = response, m = result, color = group)) + 
  geom_roc(n.cuts = 0)+
  style_roc(xlab="1 - Specificity", ylab="Sensitivity",
            theme=theme_classic )+
  theme(axis.text.x = element_text(size=18,colour="black"),
        axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title.x =element_text(face="bold",size=24,colour="black"))+ #X sixe
  theme(axis.title.y =element_text(face="bold",size=24,colour="black"))+ #Y size
  labs(title="NET, unadj and adj.age, gender")
basicplot
#ggsave("./ROC_NET_unadj_adj.tiff", width=10, height=6.98) #save plot



