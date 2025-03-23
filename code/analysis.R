### In Vitro Folate (Iskandar)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(dplyr)
library(magrittr)
library(lmerTest)
library(wesanderson)

######### New and improved methods
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
#plot(res.mme)
#qqnorm(resid(res.mme))
#qqline(resid(res.mme))
pdf("figures/boxPlot.ddi.sci.ddi.ui.pdf")
#ggline(data,x="Hour",y="AxonLength",color="Group",palette=c("firebrick3","#00AFBB"),add=c("mean_sd"),ylab="Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
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
#res.mme <- lmer(log(AxonLength) ~ Group + Hour + (1|replicate), data = data, REML = F)
res.mme <- lmer(logged ~ Group + Hour + (1|replicate), data = data, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.004626254
#plot(res.mme)
#qqnorm(resid(res.mme))
#qqline(resid(res.mme))
pdf("figures/boxPlot.ddi.sni.ddi.ui.pdf")
#ggline(data,x="Hour",y="logged",color="Group",palette=c("#E7B800","#00AFBB"),add=c("mean_sd"),ylab="Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
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
#res.mme <- lmer(log(AxonLength) ~ Group + Hour + (1|replicate), data = data, REML = F)
res.mme <- lmer(logged ~ Group + Hour + (1|replicate), data = data, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.004626254
#plot(res.mme)
#qqnorm(resid(res.mme))
#qqline(resid(res.mme))
pdf("figures/boxPlot.fa.sci.ddi.ui.pdf")
#ggline(data,x="Hour",y="logged",color="Group",palette=c("chocolate2","#00AFBB"),add=c("mean_sd"),ylab="Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
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
#res.mme <- lmer(log(AxonLength) ~ Group + Hour + (1|replicate), data = data, REML = F)
res.mme <- lmer(logged ~ Group + Hour + (1|replicate), data = data, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.004626254
#plot(res.mme)
#qqnorm(resid(res.mme))
#qqline(resid(res.mme))
pdf("figures/boxPlot.contra.ipsi.sni.pdf")
#ggline(data,x="Hour",y="logged",color="Group",palette=c("forestgreen", "#E7B800"),add=c("mean_sd"),ylab="Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ", signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
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
#ggline(data,x="Hour",y="AxonLength",color="Group",palette=pal,add=c("mean_sd"),ylab="log(Axon Length (um))") + theme(text=element_text(size=16,color="black")) + theme(plot.title = element_text(hjust=0.5))
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



############################################## 
############################################## 
############################################## 
############################################## 
############################################## 
### Old methods
############################################## 
############################################## 
############################################## 
############################################## 
############################################## 

setwd("~/Desktop/Postdoc\ Research/in\ vitro\ folate/data/")

### Untreated SNI v Untreated DDI
sni.ddi <- read.table("untreated.sni.untreated.ddi.txt",header=T)

# Line plot of data
pdf("figures/linePlot.ddi.sni.ddi.ui.pdf")
ggline(sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#E7B800","#00AFBB"),add=c("mean_se","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black"))
dev.off()

# Bar plot of data
pdf("figures/barPlot.ddi.sni.ddi.ui.pdf")
sni.ddi2$Hour <- paste0("H",sni.ddi$Hour)
ggbarplot(sni.ddi2,x="Hour",y="AxonLength",fill="Group",palette=c("#E7B800","#00AFBB"),add=c("mean_se","point"),ylab="Average Axon Length (um)",position=position_dodge(0.9))
dev.off()

# Box plot of data
ggboxplot(sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "#E7B800"),ylab="Average Axon Length (um)")

# Linear model of group and time on axon length
res <- lm(AxonLength ~ Group + Hour, data = sni.ddi)
pGroup <- signif(summary(res)$coef[2,4], digits = 2)

# Mixed effect models
sni.ddi$replicate <- c("A","B","C","D","E","A","B","C","A","B","C","A","B","C","A","B","C","D","A","B","C","A","B","C",
"Z","X","Q","Z","X","Q","Z","X","Q","Z","X","Q","Z","X","Q","Z","X","Q","Z","X","Q")
res.mme <- lmer(AxonLength ~ Group + Hour + (1|replicate), data = sni.ddi, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 1.765173e-08

# Line plot again with calculated P-value
pdf("figures/linePlot.ddi.sni.ddi.ui.pdf")
#ggline(sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#E7B800","#00AFBB"),add=c("mean_se","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",pGroup)) + theme(plot.title = element_text(hjust=0.5))
ggline(sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#E7B800","#00AFBB"),add=c("mean_sd","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()


# Two-way ANOVA of Group and Hours
sni.ddi$Group <- factor(sni.ddi$Group)
sni.ddi$Hour <- factor(sni.ddi$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
res.aov <- aov(AxonLength ~ Group * Hour,data=sni.ddi)
#TukeyHSD(res.aov)
# H5 p.adj = 0.99
# H10 p.adj = 0.336
# H13 p.adj = 0.96
# H17 p.adj = 0.149
# H24 p.adj = 0.0404
# H36 p.adj = 0.49
# H48 p.adj = 0.00067
# Paired T test
pwc <- sni.ddi %>%
group_by(Hour) %>%
pairwise_t_test(
#AxonLength ~ Group, pool.sd = F, p.adjust.methods = "bonferroni")
AxonLength ~ Group, pool.sd = T)
pwc
p.adjust(pwc$p,method="bonferroni")
# H5 p.adj = 0.162
# H10 p.adj = 0.027
# H13 p.adj = 0.273
# H17 p.adj = 0.067
# H24 p.adj = 0.034
# H36 p.adj = 0.051
# H48 p.adj = 0.101
# Wilcoxon test with bonferroni correction
df <- sni.ddi
df <- tbl_df(df)
res <- df %>% group_by(Hour) %>%
do(w = wilcox.test(AxonLength~Group, data=., paired=FALSE)) %>%
summarise(Hour, Wilcox = w$p.value)
res
p.adjust(res$Wilcox,method="bonferroni")





### Untreated SCI v Untreated DDI
sci.ddi <- read.table("untreated.sci.untreated.ddi.txt",header=T)

# Line plot of data
pdf("figures/linePlot.ddi.sci.ddi.ui.pdf")
ggline(sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("firebrick3","#00AFBB" ),add=c("mean_se","point"),ylab="Average Axon Length (um)")
dev.off()

# Bar plot of data
pdf("figures/barPlot.ddi.sci.ddi.ui.pdf")
sci.ddi2$Hour <- paste0("H",sci.ddi$Hour)
ggbarplot(sci.ddi2,x="Hour",y="AxonLength",fill="Group",palette=c("firebrick3","#00AFBB" ),add=c("mean_se","point"),ylab="Average Axon Length (um)",position=position_dodge(0.9))
dev.off()

# Box plot of data
ggboxplot(sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "firebrick3"),ylab="Average Axon Length (um)")

# Linear model of group and time on axon length
res <- lm(AxonLength ~ Group + Hour, data = sci.ddi)
pGroup <- signif(summary(res)$coef[2,4], digits = 2)

# Mixed effect models
sci.ddi$replicate <- c("A", "B", "C", "D", "E", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C",
"Z", "X", "Q", "W", "Z", "X", "Q", "W", "Z", "X", "Q", "Z", "X", "Q", "Z", "X", "Q", "W", "Z", "X", "Q", "Z", "X", "Q")
res.mme <- lmer(AxonLength ~ Group + Hour + (1|replicate), data = sci.ddi, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.00082

# Line plot again with calculated P-value
pdf("figures/linePlot.ddi.sci.ddi.ui.pdf")
#ggline(sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("firebrick3","#00AFBB" ),add=c("mean_se","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",pGroup)) + theme(plot.title = element_text(hjust=0.5))
ggline(sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("chocolate2","#00AFBB"),add=c("mean_sd","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()

# Two-way ANOVA of Group and Hours
sci.ddi$Group <- factor(sci.ddi$Group)
sci.ddi$Hour <- factor(sci.ddi$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
res.aov <- aov(AxonLength ~ Group * Hour,data=sci.ddi)
TukeyHSD(res.aov)
# H5 p.adj = 1
# H10 p.adj = 1
# H13 p.adj = 0.92
# H17 p.adj = 0.33
# H24 p.adj = 0.15
# H36 p.adj = 0.014
# H48 p.adj = 0.000059
pwc <- sci.ddi %>%
group_by(Hour) %>%
pairwise_t_test(
#AxonLength ~ Group, pool.sd = F, p.adjust.methods = "bonferroni")
AxonLength ~ Group, pool.sd = T)
pwc
p.adjust(pwc$p,method="bonferroni")
# Wilcoxon test with bonferroni correction
df <- sci.ddi
df <- tbl_df(df)
res <- df %>% group_by(Hour) %>%
do(w = wilcox.test(AxonLength~Group, data=., paired=FALSE)) %>%
summarise(Hour, Wilcox = w$p.value)
res
p.adjust(res$Wilcox,method="bonferroni")




### FA SCI v Untreated DDI
fa.sci.ddi <- read.table("fa.sci.untreated.ddi.txt",header=T)

# Line plot of data
#pdf("figures/linePlot.fa.sci.ddi.ui.pdf")
#ggline(fa.sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("chocolate2","#00AFBB"),add=c("mean_se","point"),ylab="Average Axon Length (um)")
#dev.off()

# Bar plot of data
pdf("figures/barPlot.fa.sci.ddi.ui.pdf")
fa.sci.ddi2$Hour <- paste0("H",fa.sci.ddi$Hour)
ggbarplot(fa.sci.ddi2,x="Hour",y="AxonLength",fill="Group",palette=c("chocolate2","#00AFBB"),add=c("mean_se","point"),ylab="Average Axon Length (um)",position=position_dodge(0.9))
dev.off()

# Box plot of data
#ggboxplot(fa.sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "chocolate2"),ylab="Average Axon Length (um)")

# Linear model of group and time on axon length
res <- lm(AxonLength ~ Group + Hour, data = fa.sci.ddi)
pGroup <- signif(summary(res)$coef[2,4], digits = 2)

# Mixed effect models
fa.sci.ddi$replicate <- c("A", "B", "C", "D", "E", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C",
"Z", "X", "Q", "W", "R", "Z", "X", "Q", "W", "R", "Z", "X", "Q", "W", "R", "Z", "X", "Q", "W", "R", "T", "Z", "X", "Q", "W", "Z", "X", "Q", "W", "Z", "X", "Q", "W")
res.mme <- lmer(AxonLength ~ Group + Hour + (1|replicate), data = fa.sci.ddi, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.5027892

# Line plot again with calculated P-value
pdf("figures/linePlot.fa.sci.ddi.ui.pdf")
#ggline(fa.sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("chocolate2","#00AFBB"),add=c("mean_se","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",pGroup)) + theme(plot.title = element_text(hjust=0.5))
ggline(fa.sci.ddi,x="Hour",y="AxonLength",color="Group",palette=c("chocolate2","#00AFBB"),add=c("mean_sd","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()

# Two-way ANOVA of Group and Hours
fa.sci.ddi$Group <- factor(fa.sci.ddi$Group)
fa.sci.ddi$Hour <- factor(fa.sci.ddi$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
res.aov <- aov(AxonLength ~ Group * Hour,data=fa.sci.ddi)
TukeyHSD(res.aov)
# H5 p.adj = 1
# H10 p.adj = 1
# H13 p.adj = 1
# H17 p.adj = 1
# H24 p.adj = 0.999
# H36 p.adj = 0.992
# H48 p.adj = 1
pwc <- fa.sci.ddi %>%
group_by(Hour) %>%
pairwise_t_test(
#AxonLength ~ Group, pool.sd = F, p.adjust.methods = "bonferroni")
AxonLength ~ Group, pool.sd = T)
pwc
p.adjust(pwc$p,method="bonferroni")
# Wilcox test
df <- fa.sci.ddi
df <- tbl_df(df)
res <- df %>% group_by(Hour) %>%
do(w = wilcox.test(AxonLength~Group, data=., paired=FALSE)) %>%
summarise(Hour, Wilcox = w$p.value)
res
p.adjust(res$Wilcox,method="bonferroni")




### FA SNI v Untreated DDI
fa.sni.ddi <- read.table("fa.sni.untreated.ddi.txt",header=T)

# Line plot of data
#pdf("figures/linePlot.fa.sni.ddi.ui.pdf")
#ggline(fa.sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "darkgoldenrod3"),add=c("mean_se","point"),ylab="Average Axon Length (um)")
#dev.off()

# Bar plot of data
pdf("figures/barPlot.fa.sni.ddi.ui.pdf")
fa.sni.ddi2$Hour <- paste0("H",fa.sni.ddi$Hour)
ggbarplot(fa.sni.ddi2,x="Hour",y="AxonLength",fill="Group",palette=c("#00AFBB", "darkgoldenrod3"),add=c("mean_se","point"),ylab="Average Axon Length (um)",position=position_dodge(0.9))
dev.off()

# Box plot of data
#ggboxplot(fa.sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "darkgoldenrod3"),ylab="Average Axon Length (um)")

# Linear model of group and time on axon length
res <- lm(AxonLength ~ Group + Hour, data = fa.sni.ddi)
pGroup <- signif(summary(res)$coef[2,4], digits = 2)

# Mixed effect models
fa.sni.ddi$replicate <- c("A", "B", "C", "D", "E", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C",
"Z", "X", "Q", "Z", "X", "Q", "Z", "X", "Z", "X", "Z", "X", "Z", "X", "Q", "Z", "X")
res.mme <- lmer(AxonLength ~ Group + Hour + (1|replicate), data = fa.sni.ddi, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 7.930589e-08

# Line plot again with calculated P-value
pdf("figures/linePlot.fa.sni.ddi.ui.pdf")
#ggline(fa.sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "darkgoldenrod3"),add=c("mean_se","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",pGroup)) + theme(plot.title = element_text(hjust=0.5))
ggline(fa.sni.ddi,x="Hour",y="AxonLength",color="Group",palette=c("chocolate2","#00AFBB"),add=c("mean_sd","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()

# Two-way ANOVA of Group and Hours
fa.sni.ddi$Group <- factor(fa.sni.ddi$Group)
fa.sni.ddi$Hour <- factor(fa.sni.ddi$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
res.aov <- aov(AxonLength ~ Group * Hour,data=fa.sni.ddi)
TukeyHSD(res.aov)
# H5 p.adj = 0.99
# H10 p.adj = 0.053
# H13 p.adj = 0.99
# H17 p.adj = 0.141
# H24 p.adj = 0.00002
# H36 p.adj = 0.22
# H48 p.adj = 0.015
pwc <- fa.sni.ddi %>%
group_by(Hour) %>%
pairwise_t_test(
#AxonLength ~ Group, pool.sd = F, p.adjust.methods = "bonferroni")
AxonLength ~ Group, pool.sd = T)
pwc
p.adjust(pwc$p,method="bonferroni")
# Wilcox test
df <- fa.sni.ddi
df <- tbl_df(df)
res <- df %>% group_by(Hour) %>%
do(w = wilcox.test(AxonLength~Group, data=., paired=FALSE)) %>%
summarise(Hour, Wilcox = w$p.value)
res
p.adjust(res$Wilcox,method="bonferroni")




### FA UI v Untreated UI
fa.ddi.ui <- read.table("fa.ui.untreated.ddi.txt",header=T)

# Line plot of data
#pdf("figures/linePlot.fa.ui.ddi.ui.pdf")
#ggline(fa.ddi.ui,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "magenta"),add=c("mean_se","point"),ylab="Average Axon Length (um)")
#dev.off()

# Bar plot of data
pdf("figures/barPlot.fa.ui.ddi.ui.pdf")
fa.ddi.ui2$Hour <- paste0("H",fa.ddi.ui$Hour)
ggbarplot(fa.ddi.ui2,x="Hour",y="AxonLength",fill="Group",palette=c("#00AFBB", "magenta"),add=c("mean_se","point"),ylab="Average Axon Length (um)",position=position_dodge(0.9))
dev.off()

# Box plot of data
#ggboxplot(fa.ddi.ui,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "magenta"),ylab="Average Axon Length (um)")

# Linear model of group and time on axon length
res <- lm(AxonLength ~ Group + Hour, data = fa.ddi.ui)
pGroup <- signif(summary(res)$coef[2,4], digits = 2)

# Mixed effect models
fa.ddi.ui$replicate <- c("A", "B", "C", "D", "E", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C",
"Z", "X", "Q", "W", "R", "T", "Z", "X", "Q", "W", "R", "Z", "X", "Q", "W", "R", "Z", "X", "Q", "W", "Z", "X", "Q", "Z", "X", "Q", "Z", "X", "Q")
res.mme <- lmer(AxonLength ~ Group + Hour + (1|replicate), data = fa.ddi.ui, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.1479175

# Line plot again with calculated P-value
pdf("figures/linePlot.fa.ui.ddi.ui.pdf")
#ggline(fa.ddi.ui,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "magenta"),add=c("mean_se","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",pGroup)) + theme(plot.title = element_text(hjust=0.5))
ggline(fa.ddi.ui,x="Hour",y="AxonLength",color="Group",palette=c("#00AFBB", "magenta"),add=c("mean_sd","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ", signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()

# Two-way ANOVA of Group and Hours
fa.ddi.ui$Group <- factor(fa.ddi.ui$Group)
fa.ddi.ui$Hour <- factor(fa.ddi.ui$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
res.aov <- aov(AxonLength ~ Group * Hour,data=fa.ddi.ui)
TukeyHSD(res.aov)
# H5 p.adj = 1
# H10 p.adj = 1
# H13 p.adj = 1
# H17 p.adj = 0.999
# H24 p.adj = 1
# H36 p.adj = 0.42
# H48 p.adj = 0.0045
pwc <- fa.ddi.ui %>%
group_by(Hour) %>%
pairwise_t_test(
#AxonLength ~ Group, pool.sd = F, p.adjust.methods = "bonferroni")
AxonLength ~ Group, pool.sd = T)
pwc
p.adjust(pwc$p,method="bonferroni")
# Wilcox test
df <- fa.ddi.ui
df <- tbl_df(df)
res <- df %>% group_by(Hour) %>%
do(w = wilcox.test(AxonLength~Group, data=., paired=FALSE)) %>%
summarise(Hour, Wilcox = w$p.value)
res
p.adjust(res$Wilcox,method="bonferroni")




### Untreated SNI Ipsilateral v Untreated SNI Contralateral
ipsi.contra.sni <- read.table("untreated.sni.ipsi.contra.txt",header=T)

# Line plot of data
#pdf("figures/linePlot.ddi.sni.ipsi.contra.pdf")
#ggline(ipsi.contra.sni,x="Hour",y="AxonLength",color="Group",palette=c("#E7B800", "forestgreen"),add=c("mean_se","point"),ylab="Average Axon Length (um)")
#dev.off()

# Bar plot of data
pdf("figures/barPlot.ddi.sni.ipsi.contra.pdf")
ipsi.contra.sni2$Hour <- paste0("H",ipsi.contra.sni$Hour)
ggbarplot(ipsi.contra.sni2,x="Hour",y="AxonLength",fill="Group",palette=c("#E7B800", "forestgreen"),add=c("mean_se","point"),ylab="Average Axon Length (um)",position=position_dodge(0.9))
dev.off()

# Box plot of data
#ggboxplot(ipsi.contra.sni,x="Hour",y="AxonLength",color="Group",palette=c("#E7B800", "plum3"),ylab="Average Axon Length (um)")

# Linear model of group and time on axon length
res <- lm(AxonLength ~ Group + Hour, data = ipsi.contra.sni)
pGroup <- signif(summary(res)$coef[2,4], digits = 2)

# Mixed effect models
ipsi.contra.sni$replicate <- c("A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", 
"Z", "X", "Z", "X", "Z", "X", "Z", "X", "Z", "X", "Z", "X", "Z", "X")
res.mme <- lmer(AxonLength ~ Group + Hour + (1|replicate), data = ipsi.contra.sni, REML = F)
pGroup.mme <- summary(res.mme)$coefficients[2,5] # 0.8730416

# Line plot again with calculated P-value
pdf("figures/linePlot.ddi.sni.ipsi.contra.pdf")
#ggline(ipsi.contra.sni,x="Hour",y="AxonLength",color="Group",palette=c("#E7B800", "forestgreen"),add=c("mean_se","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ",pGroup)) + theme(plot.title = element_text(hjust=0.5))
ggline(ipsi.contra.sni,x="Hour",y="AxonLength",color="Group",palette=c("#E7B800", "forestgreen"),add=c("mean_sd","point"),ylab="Average Axon Length (um)") + theme(text=element_text(size=16,color="black")) + labs(title=paste0("p = ", signif(pGroup.mme, digits=1))) + theme(plot.title = element_text(hjust=0.5))
dev.off()

# Two-way ANOVA of Group and Hours
ipsi.contra.sni$Group <- factor(ipsi.contra.sni$Group)
ipsi.contra.sni$Hour <- factor(ipsi.contra.sni$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
res.aov <- aov(AxonLength ~ Group * Hour,data=ipsi.contra.sni)
TukeyHSD(res.aov)
# H5 p.adj = 1
# H10 p.adj = 0.99
# H13 p.adj = 0.99
# H17 p.adj = 0.999
# H24 p.adj = 0.98
# H36 p.adj = 0.99
# H48 p.adj = 0.99
pwc <- ipsi.contra.sni %>%
group_by(Hour) %>%
pairwise_t_test(
#AxonLength ~ Group, pool.sd = F, p.adjust.methods = "bonferroni")
AxonLength ~ Group, pool.sd = T)
pwc
p.adjust(pwc$p,method="bonferroni")
# Wilcox test
df <- ipsi.contra.sni
df <- tbl_df(df)
res <- df %>% group_by(Hour) %>%
do(w = wilcox.test(AxonLength~Group, data=., paired=FALSE)) %>%
summarise(Hour, Wilcox = w$p.value)
res
p.adjust(res$Wilcox,method="bonferroni")


### Myelin plated Embryonic vs FA DRG vs DDI DRG
emb.fa.ddi <- read.table("myelin.embryo.ddi.fa.txt",header=T)

# Line plot of data
#library(wesanderson)
#pal <- as.vector(wes_palette("Darjeeling1",3,type="continuous"))
#pdf("figures/linePlot.myelin.emb.fa.ddi.pdf")
#ggline(emb.fa.ddi,x="Hour",y="AxonLength",color="Group",palette=pal,add=c("mean_se","point"),ylab="Average Axon Length (um)")
#dev.off()

# Bar plot of data
pdf("figures/barPlot.myelin.emb.fa.ddi.pdf")
emb.fa.ddi2$Hour <- paste0("H",emb.fa.ddi$Hour)
pal <- as.vector(wes_palette("Darjeeling1",3,type="continuous"))
ggbarplot(emb.fa.ddi2,x="Hour",y="AxonLength",fill="Group",palette=pal,add=c("mean_se","point"),ylab="Average Axon Length (um)",position=position_dodge(0.9))
dev.off()

# Two-way ANOVA of Group and Hours
emb.fa.ddi$Group <- factor(emb.fa.ddi$Group)
emb.fa.ddi$Hour <- factor(emb.fa.ddi$Hour,levels=c(5,10,13,17,24,36,48),label=c("H5","H10","H13","H17","H24","H36","H48"))
res.aov <- aov(AxonLength ~ Group * Hour,data=emb.fa.ddi)
x <- TukeyHSD(res.aov)
x <- x$'Group:Hour'
x <- as.data.frame(x)
x[which(x$'p adj'<0.05),]
# H10 (FA Adult v Emb) = 0.35
# H13 (FA Adult v Emb) = 0.01
# H17 (FA Adult v Emb) = 0.21
# H24 (FA Adult v Emb) = 2.8e-5
# H36 (FA Adult v Emb) = 0.0005
# H48 (FA Adult v Emb) = 0.001
# H10 (Untrt Adult v Emb) = 0.042
# H13 (Untrt Adult v Emb) = 0.001
# H17 (Untrt Adult v Emb) = 0.15
# H24 (Untrt Adult v Emb) = 2.8e-5
# H36 (Untrt Adult v Emb) = 1.2e-6
# H48 (Untrt Adult v Emb) = 0.05
# H10 (Untrt Adult v FA Adult) = 0.999
# H13 (Untrt Adult v FA Adult) = 0.999
# H17 (Untrt Adult v FA Adult) = 1
# H24 (Untrt Adult v FA Adult) = 1
# H36 (Untrt Adult v FA Adult) = 0.96
# H48 (Untrt Adult v FA Adult) = 1.2e-7
pwc <- emb.fa.ddi %>%
group_by(Hour) %>%
pairwise_t_test(
#AxonLength ~ Group, pool.sd = F, p.adjust.methods = "bonferroni")
AxonLength ~ Group, pool.sd = T)
pwc
p.adjust(pwc$p,method="bonferroni")
# Wilcox test
#df <- emb.fa.ddi
#df <- tbl_df(df)
#res <- df %>% group_by(Hour) %>%
#do(w = wilcox.test(AxonLength~Group, data=., paired=FALSE)) %>%
#summarise(Hour, Wilcox = w$p.value)
#res
#p.adjust(res$Wilcox,method="bonferroni")

# Pairwise t test with Bonferroni correction
x <- emb.fa.ddi[which(emb.fa.ddi$Hour==10),]
res10 <- pairwise.t.test(x$AxonLength, x$Group, p.adjust="bonferroni")
x <- emb.fa.ddi[which(emb.fa.ddi$Hour==13),]
res13 <- pairwise.t.test(x$AxonLength, x$Group, p.adjust="bonferroni")
x <- emb.fa.ddi[which(emb.fa.ddi$Hour==17),]
res17 <- pairwise.t.test(x$AxonLength, x$Group, p.adjust="bonferroni")
x <- emb.fa.ddi[which(emb.fa.ddi$Hour==24),]
res24 <- pairwise.t.test(x$AxonLength, x$Group, p.adjust="bonferroni")
x <- emb.fa.ddi[which(emb.fa.ddi$Hour==36),]
res36 <- pairwise.t.test(x$AxonLength, x$Group, p.adjust="bonferroni")
x <- emb.fa.ddi[which(emb.fa.ddi$Hour==48),]
res48 <- pairwise.t.test(x$AxonLength, x$Group, p.adjust="bonferroni")
