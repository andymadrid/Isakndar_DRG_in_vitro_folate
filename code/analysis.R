### In Vitro Folate (Iskandar)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(dplyr)
library(magrittr)
library(lmerTest)
library(wesanderson)

######### New method
setwd("~/Desktop/Postdoc Research/in vitro folate/data2/")

# load in all data
sci.ddi <- read.table("DDI_SCI.txt", header = T)
sni.ddi <- read.table("DDI_SNI.txt", header = T)
ui.ddi <- read.table("DDI_UI.txt", header = T)
sci.fa <- read.table("FA_SCI.txt", header = T)
sni.fa <- read.table("FA_SNI.txt", header = T)
ui.fa <- read.table("FA_UI.txt", header = T)
contra <- read.table("DDI_SNI_Contra.txt", header = T)
myelin.fa <- read.table("Myelin_FA.txt", header = T)
myelin.ddi <- read.table("Myelin_DDI.txt", header = T)
myelin.embryo <- read.table("Myelin_Embryo.txt", header = T)
myelin.embryo$Hour <- gsub(12, 13, myelin.embryo$Hour)

# DDI SCI vs DDI UI
data <- rbind(sci.ddi, ui.ddi)
data$Group <- gsub("DDI_SCI", "SCI", data$Group)
data$Group <- gsub("DDI_UI", "UI", data$Group)
data$Group <- factor(data$Group, levels = c("SCI", "UI"))
new <- c()
for (i in 1:nrow(data)) {
if (data[i,"AxonLength"] > 0) {
new <- rbind(new,data[i,])}}
data <- new
data$logged <- log(data$AxonLength)
#res.mme <- lmer(log(AxonLength) ~ Group + Hour + (1|replicate), data = data, REML = F)
res.mme <- lmer(logged ~ Group + Hour + (1|replicate), data = data, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.004626254
pdf("figures/boxPlot.ddi.sci.ddi.ui.pdf")
ggboxplot(data,x="Hour",y="logged",fill="Group",palette=c("firebrick3","#00AFBB"),ylab="log(Axon Length (um))") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()
data %>% count(Hour, replicate, sort = T)
as.data.frame(data %>% group_by(Hour, Group) %>% summarise(across(.cols = AxonLength,list(mean = mean, sd = sd))))

# DDI SNI vs DDI UI
data <- rbind(sni.ddi, ui.ddi)
data$Group <- gsub("DDI_SNI", "SNI", data$Group)
data$Group <- gsub("DDI_UI", "UI", data$Group)
data$Group <- factor(data$Group, levels = c("SNI", "UI"))
new <- c()
for (i in 1:nrow(data)) {
if (data[i,"AxonLength"] > 0) {
new <- rbind(new,data[i,])}}
data <- new
data$logged <- log(data$AxonLength)
res.mme <- lmer(logged ~ Group + Hour + (1|replicate), data = data, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.004626254
pdf("figures/boxPlot.ddi.sni.ddi.ui.pdf")
ggboxplot(data,x="Hour",y="logged",fill="Group",palette=c("#E7B800","#00AFBB"),ylab="log(Axon Length (um))") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()
data %>% count(Hour, replicate, sort = T)
as.data.frame(data %>% group_by(Hour, Group) %>% summarise(across(.cols = AxonLength,list(mean = mean, sd = sd))))

# FA SCI vs DDI UI
data <- rbind(sci.fa, ui.ddi)
data$Group <- gsub("FA_SCI", "FA SCI", data$Group)
data$Group <- gsub("DDI_UI", "UI", data$Group)
data$Group <- factor(data$Group, levels = c("FA SCI", "UI"))
new <- c()
for (i in 1:nrow(data)) {
if (data[i,"AxonLength"] > 0) {
new <- rbind(new,data[i,])}}
data <- new
data$logged <- log(data$AxonLength)
res.mme <- lmer(logged ~ Group + Hour + (1|replicate), data = data, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.004626254
pdf("figures/boxPlot.fa.sci.ddi.ui.pdf")
ggboxplot(data,x="Hour",y="logged",fill="Group",palette=c("chocolate2","#00AFBB"),ylab="log(Axon Length (um))") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()
data %>% count(Hour, replicate, sort = T)
as.data.frame(data %>% group_by(Hour, Group) %>% summarise(across(.cols = AxonLength,list(mean = mean, sd = sd))))

# Contra vs Ipsi
data <- rbind(sni.ddi, contra)
data$Group <- gsub("DDI_SNI", "Ipsi", data$Group)
data$Group <- factor(data$Group, levels = c("Contra", "Ipsi"))
new <- c()
for (i in 1:nrow(data)) {
if (data[i,"AxonLength"] > 0) {
new <- rbind(new,data[i,])}}
data <- new
data$logged <- log(data$AxonLength)
res.mme <- lmer(logged ~ Group + Hour + (1|replicate), data = data, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.004626254
pdf("figures/boxPlot.contra.ipsi.sni.pdf")
ggboxplot(data,x="Hour",y="logged",fill="Group",palette=c("forestgreen", "#E7B800"),ylab="log(Axon Length (um))") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ", signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()
data %>% count(Hour, replicate, sort = T)
as.data.frame(data %>% group_by(Hour, Group) %>% summarise(across(.cols = AxonLength,list(mean = mean, sd = sd))))

# Myelin data
data <- rbind(myelin.embryo, myelin.fa, myelin.ddi)
data$Group <- gsub("FA", "FA_Adult", data$Group)
data$Group <- gsub("Embryo", "Embryonic", data$Group)
data$Group <- gsub("DDI", "Untrt_Adult", data$Group)
data$Group <- factor(data$Group, levels = c("Embryonic", "FA_Adult", "Untrt_Adult"))
new <- c()
for (i in 1:nrow(data)) {
if (data[i,"AxonLength"] > 0) {
new <- rbind(new,data[i,])}}
data <- new
s <- which(data$Hour == 72)
data <- data[ !(rownames(data) %in% s),]
data$Group <- factor(data$Group)
data$Hour <- factor(data$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
data <- data[!is.na(data$Hour),]
data$logged <- log(data$AxonLength)
res.aov <- aov(logged ~ Group * Hour, data = data)
x <- TukeyHSD(res.aov)
x <- x$'Group:Hour'
x <- as.data.frame(x)
x[which(x$'p adj'<0.05),]
pal <- as.vector(wes_palette("Darjeeling1",3,type="continuous"))
pdf("figures/boxPlot.myelin.emb.fa.ddi.pdf")
ggboxplot(data,x="Hour",y="logged",fill="Group",palette=pal,ylab="Axon Length (um)") + theme(text=element_text(size=16,color="black")) + theme(plot.title = element_text(hjust=0.5))
dev.off()
data %>% count(Hour, replicate, sort = T)
as.data.frame(data %>% group_by(Hour, Group) %>% summarise(across(.cols = AxonLength,list(mean = mean, sd = sd))))

# DNMT3B analysis
control_noFA <- read.table("dnmt3b_48h_controlSI_noFA.txt", header = F)
control_FA <- read.table("dnmt3b_48h_controlSI_FA.txt", header = F)
dnmt3b_noFA <- read.table("dnmt3b_48h_dnmt3bSI_noFA.txt", header = F)
dnmt3b_FA <- read.table("dnmt3b_48h_dnmt3bSI_FA.txt", header = F)
control_noFA$Group <- "Control siRNA" 
control_FA$Group <- "Control siRNA + FA"
dnmt3b_noFA$Group <- "Dnm3b siRNA"
dnmt3b_FA$Group <- "Dnmt3b siRNA + FA"
data <- rbind(control_noFA, control_FA, dnmt3b_noFA, dnmt3b_FA)
data$Group <- factor(data$Group, levels = c("Control siRNA", "Dnm3b siRNA", "Control siRNA + FA", "Dnmt3b siRNA + FA"))
new <- c()
for (i in 1:nrow(data)) {
if (data[i,"V1"] > 0) {
new <- rbind(new,data[i,])}}
data <- new
data$logged <- log(data$V1)
as.data.frame(data %>% group_by(Group) %>% summarise(across(.cols = logged,list(mean = mean, sd = sd))))
data$Group <- factor(data$Group, levels = c("Control siRNA" , "Dnm3b siRNA", "Control siRNA + FA", "Dnmt3b siRNA + FA"))
pdf("figures/boxPlot.dnmt3b.pdf")
ggboxplot(data,x="Group",y="logged",fill="Group",palette=c("firebrick3", "darkgoldenrod3", "forestgreen", "mediumblue"),ylab="log(Axon Length (um))") + theme(text=element_text(size=16,color="black")) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") + xlab("")
dev.off()

