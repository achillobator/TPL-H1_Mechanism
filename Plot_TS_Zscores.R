library(ggplot2)
library(ggthemes)
library(nlme)
library(mgcv)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(bigPint)
library(dplyr)
library(plotly)
library(Hmisc)
library(pheatmap)
library(ggrepel)
#install.packages("psych", dependencies=TRUE)
#install.packages("psychTools")
library(psychTools)
install.packages("psych",repos="https://personality-project.org/r/",type="source")

ORFmap <- read.csv(file = "TS_quadruplicate.csv", stringsAsFactors = F)
ORFmap = ORFmap[!duplicated(ORFmap$TSQ.ID),1:7]
#names(ORFmap) <- c("TSQ.ID")

dat.3669.1.ts.ratiozscore <- read.csv(file = "3669-1_T2_TSQIDBavg_Z.csv", stringsAsFactors = F)
trim.3669.1.ts <- dat.3669.1.ts.ratiozscore[ , c("TSQ.ID", "Zscore")]
dat.merge <- merge(ORFmap, trim.3669.1.ts, by = "TSQ.ID", all.y=T)
names(dat.merge)[8] <- c("3669-1.TS")
 
dat.3669.1.NAA.ts.ratiozscore <- read.csv(file = "3669-1_TS_NAA_TSQIDBavg_Z.csv", stringsAsFactors = F)
trim.3669.1.NAA.ts <- dat.3669.1.NAA.ts.ratiozscore[ , c("TSQ.ID", "Zscore")]
dat.merge2 <- merge(dat.merge, trim.3669.1.NAA.ts, by = "TSQ.ID", all.y=T)
names(dat.merge2)[9] <- c("3669-1.NAA.TS")

dat.3669.2.TS.ratiozscore <- read.csv(file = "3669-2_TS_TSQIDBavg_Z.csv", stringsAsFactors = F)
trim.3669.2.TS <- dat.3669.2.TS.ratiozscore[ , c("TSQ.ID", "Zscore")]
dat.merge3 <- merge(dat.merge2, trim.3669.2.TS, by = "TSQ.ID", all.y=T)
names(dat.merge3)[10] <- c("3669-2.TS")

dat.3669.2.NAA.TS.ratiozscore <- read.csv(file = "3669-2_TS_NAA_TSQIDBavg_Z.csv", stringsAsFactors = F)
trim.3669.2.NAA.TS <- dat.3669.2.NAA.TS.ratiozscore[ , c("TSQ.ID", "Zscore")]
dat.merge4 <- merge(dat.merge3, trim.3669.2.NAA.TS, by = "TSQ.ID", all.y=T)
names(dat.merge4)[11] <- c("3669-2.NAA.TS")

dat.3670.1.TS.ratiozscore <- read.csv(file = "3670-1_TS_TSQIDBavg_Z.csv", stringsAsFactors = F)
trim.3670.1.TS <- dat.3670.1.TS.ratiozscore[ , c("TSQ.ID", "Zscore")]
dat.merge5 <- merge(dat.merge4, trim.3670.1.TS, by = "TSQ.ID", all.y=T)
names(dat.merge5)[12] <- c("3670-1.TS")

dat.3670.1.NAA.TS.ratiozscore <- read.csv(file = "3670-1_TS_NAA_TSQIDBavg_Z.csv", stringsAsFactors = F)
trim.3670.1.NAA.TS <- dat.3670.1.NAA.TS.ratiozscore[ , c("TSQ.ID", "Zscore")]
dat.merge6 <- merge(dat.merge5, trim.3670.1.NAA.TS, by = "TSQ.ID", all.y=T)
names(dat.merge6)[13] <- c("3670-1.NAA.TS")

dat.3670.2.TS.ratiozscore <- read.csv(file = "3670-2_TS_TSQIDBavg_Z.csv", stringsAsFactors = F)
trim.3670.2.TS <- dat.3670.2.TS.ratiozscore[ , c("TSQ.ID", "Zscore")]
dat.merge7 <- merge(dat.merge6, trim.3670.2.TS, by = "TSQ.ID", all.y=T)
names(dat.merge7)[14] <- c("3670-2.TS")

### Plot heatmap

colfunc <- colorRampPalette(c("#007AF4","black", "yellow"))
breaks <- c(-10, seq(-5,5, by=0.1), 10)

pheatmap(dat.merge7[,8:14],
         cluster_rows=F, cluster_col=F,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=6,
         fontsize_row=, 
         fontsize_col=6, 
         cellwidth=10, 
         cellheight=, 
         main=)

colfunc <- colorRampPalette(c("#007AF4","black", "yellow"))

dat.merge_dropN188auxin <- dat.merge7[ -c(13) ]
plot(dat.merge_dropN188auxin[,8:13], 
     pch = 10, cex = .1,alpha=0.3,
     main="TS Experimental Correlation",
     lower.panel = NULL)

pairs.panels(dat.merge_dropN188auxin[,8:13],
             pch = 20, cex = 1, alpha = 0.5)

pairs(dat.merge_dropN188auxin[,8:13],
             pch = 20, cex = .1, alpha = 0.3, 
             lower.panel = NULL)

union70.dat <- as.matrix(read.csv(file = "70_union_up.txt", stringsAsFactors = F))
dat.union70 <- subset(dat.merge7, TSQ.ID %in% union70.dat)
rownames(dat.union70) = dat.union70[,1]
pheatmap(dat.union70[,8:14],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=6,
         fontsize_row=, 
         fontsize_col=6, 
         cellwidth=10, 
         cellheight=, 
         main=)


union69.dat <- as.matrix(read.csv(file = "69_union.txt", stringsAsFactors = F))
dat.union69 <- subset(dat.merge7, TSQ.ID %in% union69.dat)
rownames(dat.union69) = dat.union69[,1]
pheatmap(dat.union69[,8:14],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=6,
         fontsize_row=, 
         fontsize_col=6, 
         cellwidth=10, 
         cellheight=, 
         main=)

union69.70.up.dat <- as.matrix(read.csv(file = "all_up_union", stringsAsFactors = F))
dat.union69.70.up <- na.omit(subset(dat.merge7, TSQ.ID %in% union69.70.up.dat))
tsFFs <- as.matrix(read.csv("TS_frequent_fliers.txt", stringsAsFactors = F))
dat.union69.70.up <- subset(dat.union69.70.up, !gene %in% tsFFs)

rownames(dat.union69.70.up) = dat.union69.70.up[,7]
pheatmap(dat.union69.70.up2[,8:14],
         cluster_rows=T, cluster_col=F,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=6,
         fontsize_row=, 
         fontsize_col=6, 
         cellwidth=10, 
         cellheight=,
         cutree_rows = 3,
         cutree_cols = 4, 
         main=)

dat.union69.70.up2 <- dat.union69.70.up[, c(1, 2, 3, 4, 5, 6, 7, 12, 14, 13, 8, 10, 9, 11)]
dat.union69.70.up2_apex <- merge(dat.union69.70.up2, APEX.dat, by = "gene", all.x=TRUE)
dat.union69.70.up2_cyto <- merge(dat.union69.70.up2_apex, TScytoval.H1.dat, by = "TSQ.ID", all.x=TRUE)
dat.union69.70.up2_cyto <- dat.union69.70.up2_cyto[ -c(16) ]
dat.union69.70.up2_cyto <- dat.union69.70.up2_cyto %>% 
  rename("H1-H5_Venus" = "Venus.Amedian",
         "H1-H5_tomato" =  "tdtomato3.Amedian",
         "H1-H5_Ratio" =  "Ratio")
dat.union69.70.up2_cyto2 <- merge(dat.union69.70.up2_cyto, TScytoval.N188.dat, by = "TSQ.ID", all.x=TRUE)
dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2[ -c(19) ]
dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2 %>% 
  rename("N188_Venus" = "Venus.Amedian",
         "N188_tomato" =  "tdtomato3.Amedian",
         "N188_Ratio" =  "Ratio")
dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2[-46,]
dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2[-54,]
dat.union69.70.up2_cyto3 <- dat.union69.70.up2_cyto2[ -c(10) ]
colfunc <- colorRampPalette(c("#007AF4","black", "yellow"))
breaks <- c(-10, seq(-5,5, by=0.2), 10)
rownames(dat.union69.70.up2_cyto3) = dat.union69.70.up2_cyto3[,7]
HMT_A <- pheatmap(dat.union69.70.up2_cyto3[,8:13],
         cluster_rows=T, cluster_col=F,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=6,
         fontsize_row=, 
         fontsize_col=6, 
         cellwidth=10, 
         cellheight=,
         cutree_rows = 4,
         cutree_cols = 4, 
         main=)
HMT_A 

colfunc2 <- colorRampPalette(c("#007AF4","black", "yellow"))


breaks2 = seq(from = 0, to = 250, by = 2.5)
rownames(dat.union69.70.up2_cyto2) = dat.union69.70.up2_cyto2[,7]
HMT_B <- pheatmap(dat.union69.70.up2_cyto2[,16:21],
                  cluster_rows=T, cluster_col=F,  	
                  treeheight_row=, treeheight_col=,
                  border_color='black', 
                  breaks= breaks2, 
                  color=colfunc2(length(breaks)+1), 
                  fontsize=6,
                  fontsize_row=, 
                  fontsize_col=6, 
                  cellwidth=10, 
                  cellheight=,
                  cutree_rows = 3,
                  cutree_cols = 4, 
                  main=)
HMT_B 

library(ComplexHeatmap)
BiocManager::install("ComplexHeatmap")
library(colorRamp2)
col_fun = colorRamp2(c(c(-10,-5, 0,5, 10)), c("#007AF4","#007AF4","black", "yellow","yellow"))
f1 = colorRamp2(seq(min(dat.union69.70.up2_cyto3[,8:13]), max(dat.union69.70.up2_cyto3[,8:13]), length = 3), c("#007AF4","black", "yellow"))
HT1<- Heatmap(dat.union69.70.up2_cyto3[,8:13], name = "Z-score",col = col_fun,
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        row_km = 4,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 10),
        rect_gp = gpar(col = "black", lwd = 1))
HT1

col_fun2 = colorRamp2(c(c(3,3.5,4)), c("red","orange","yellow"))
HT2<- Heatmap(dat.union69.70.up2_cyto3[,14], name = "Log2 APEX",
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col = (col_fun2),
              rect_gp = gpar(col = "black", lwd = 1))
HT2
HT1+HT2

col_fun3 = colorRamp2(c(c(0,100,200)), c("red","orange", "yellow"))
HT3<- Heatmap(dat.union69.70.up2_cyto3[,15:16], name = "Raw FL",
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col = col_fun3,
              row_names_gp = gpar(fontsize = 10),
              rect_gp = gpar(col = "black", lwd = 1))
HT3
HT1+HT2+HT3

col_fun4 = colorRamp2(c(c(0,0.49,0.5,1,2,3)), c("grey","grey","red","orange", "gold","yellow"))
HT4<- Heatmap(dat.union69.70.up2_cyto3[,17], name = "H1-H5 Ratio",
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col = col_fun4,
              rect_gp = gpar(col = "black", lwd = 1))
HT4
HT1+HT2+HT4
HT1+HT4

col_fun5 = colorRamp2(c(c(0,0.04,.041,.1,.3,.5)), c("grey","grey","red","orange","gold", "yellow"))
HT5<- Heatmap(dat.union69.70.up2_cyto3[,20], name = "N188 Ratio",
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              col = col_fun5,
              rect_gp = gpar(col = "black", lwd = 1))
HT5
HT1+HT2+HT3+HT4+HT5
HT1+HT2+HT4+HT5

### apex info
APEX.dat <- read.csv(file = "APEXinfomean.csv", stringsAsFactors = F)
APEX.dat <- APEX.dat[ , c("gene", "log2apex")]
dat.merge.apex <- merge(dat.merge7, APEX.dat, by = "gene", all.x=TRUE)

pheatmap(dat.merge.apex[,8:15],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=6,
         fontsize_row=, 
         fontsize_col=6, 
         cellwidth=10, 
         cellheight=, 
         main=)

dat.merge.apex2 <- merge(dat.merge7, APEX.dat, by = "gene", all.y=FALSE)
dat.merge.apex2_3.5 <- subset(dat.merge.apex2, log2apex > 3.5)

rownames(dat.merge.apex2_3.5) = dat.merge.apex2_3.5[,7]
breaks <- c(-10, seq(-5,5, by=0.1), 10)

pheatmap(dat.merge.apex2_3.5[,8:14],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=4,
         fontsize_row=, 
         fontsize_col=6, 
         cellwidth=10, 
         cellheight=, 
         main=)


#Venn - diagram subset heatmap
union69.70.up.dat <- as.matrix(read.csv(file = "all_up_union", stringsAsFactors = F))
dat.union69.70.up <- na.omit(subset(dat.merge7, TSQ.ID %in% union69.70.up.dat))
tsFFs <- as.matrix(read.csv("TS_frequent_fliers.txt", stringsAsFactors = F))
dat.union69.70.up <- subset(dat.union69.70.up, !gene %in% tsFFs)

specific69.dat <- as.matrix(read.csv(file = "69specificvenn.csv", stringsAsFactors = F))
dat.specific69.up <- na.omit(subset(dat.union69.70.up, TSQ.ID %in% specific69.dat))

rownames(dat.specific69.up) = dat.specific69.up[,7]
pheatmap(dat.specific69.up[,8:14],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=10,
         fontsize_row=, 
         fontsize_col=10, 
         cellwidth=10, 
         cellheight=,
         cutree_rows = 2,
         cutree_cols = 3, 
         main=)

###Adding in cytometry validation for the Union up genes:



###Adding in cytometry validation TAF genes only: 
g_TAFgenes <- c("TSA514", "TSA515", "TSA517", "TSA518", "TSA806", "TSA1494", "TSA1494", "TSA1495", "TSA1495", "TSA1496", "TSA1496", "TSA1497", "TSA1497", "TSA1498", "TSA1498", "TSA2753", "TSA2754", "TSA957", "TSA524", "TSA525", "TSA645", "TSA636", "TSA797", "TSA695", "TSA802", "TSA3112", "TSA942", "TSA509", "TSA512", "TSA513", "TSA608", "TSA638", "TSA730")
tfiid_dat <- subset(dat.merge7, (TSQ.ID %in% g_TAFgenes))
rownames(tfiid_dat) <- tfiid_dat[,7]
pheatmap(tfiid_dat[,8:14],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=10,
         fontsize_row=, 
         fontsize_col=10, 
         cellwidth=10, 
         cellheight=,
         cutree_rows = 2,
         cutree_cols = 3,
         #labels_row=rownames,
         main=)

tfiid_datb <- na.omit(tfiid_dat)
rownames(tfiid_datb) <- tfiid_datb[,7]
TF2Dtest<- Heatmap(tfiid_datb[,8:14], name = "Z-score",col = col_fun,
              cluster_rows = FALSE, 
              cluster_columns = TRUE, 
              row_km = 0,
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 10),
              na_col = "grey90",
              rect_gp = gpar(col = "black", lwd = 1))
TF2Dtest

g_TAFgenes <- c("TSA514", "TSA515", "TSA517", "TSA518", "TSA806", "TSA1494", "TSA1494", "TSA1495", "TSA1495", "TSA1496", "TSA1496", "TSA1497", "TSA1497", "TSA1498", "TSA1498", "TSA2753", "TSA2754", "TSA957", "TSA524", "TSA525", "TSA645", "TSA636", "TSA797", "TSA695", "TSA802", "TSA3112", "TSA942", "TSA509", "TSA512", "TSA513", "TSA608", "TSA638", "TSA730")
tfiid_dat_apex <- subset(dat.merge.apex, (TSQ.ID %in% g_TAFgenes))
write.csv(tfiid_dat_apex, "tfiid_dat_apex.csv")
tfiid_dat_apex2<-read.csv("tfiid_dat_apex.csv")

rownames(tfiid_dat_apex2) <- tfiid_dat_apex2[,8]
pheatmap(tfiid_dat_apex2[,9:14],
         cluster_rows=F, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=10,
         fontsize_row=, 
         fontsize_col=10, 
         cellwidth=10, 
         cellheight=,
         cutree_rows = 2,
         cutree_cols = 4,
         na_col = "grey90",
         #labels_row=rownames,
         main=)

rownames(tfiid_dat_apex2) <- tfiid_dat_apex2[,8]
tfiid_dat_apex3 <- subset(tfiid_dat_apex2, Allele != "taf10-ts34")
tfiid_dat_apex3$order <- factor(tfiid_dat_apex3$order, levels=c(1, 2, 3, 4, 5, 6, 7, 8,9,10,11,13,14,15,16,17,18,19),ordered=TRUE)
tfiid_dat_apex3 <- tfiid_dat_apex3[, c("X","gene","TSQ.ID","plate","r","c","Name","Allele","X3669.1.TS","X3669.2.TS","X3669.1.NAA.TS","X3669.2.NAA.TS","X3670.1.TS","X3670.2.DA","log2apex","order")]

TF2Dtest<- Heatmap(tfiid_dat_apex3[,9:14], name = "Z-score",col = col_fun,
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE, 
                   row_km = 0,
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 10),
                   na_col = "grey90",
                   labCol = FALSE,
                   rect_gp = gpar(col = "black", lwd = 1))
TF2Dtest
pdf("tfiid_dat_apex3.pdf")
Heatmap(tfiid_dat_apex3[,9:14], name = "Z-score",col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        row_km = 0,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 10),
        na_col = "grey90",
        rect_gp = gpar(col = "black", lwd = 1))
dev.off()




TScytoval.H1.dat <- read.csv("H1_5dat_TS.csv")
TScytoval.H1.tfiid.dat<-subset(TScytoval.H1.dat, (TSQ.ID %in% g_TAFgenes))
TScytoval.N188.dat <- read.csv("N188dat_TS.csv")
TScytoval.N188.tfiid.dat<-subset(TScytoval.N188.dat, (TSQ.ID %in% g_TAFgenes))
ts.compiled.dat <- merge(tfiid_dat_apex2, TScytoval.H1.tfiid.dat, by="TSQ.ID", all.x=TRUE)
ts.compiled.dat <- merge(ts.compiled.dat, TScytoval.N188.tfiid.dat, by="TSQ.ID", all.x=TRUE)
write.csv(ts.compiled.dat, "ts.compiled.dat.csv")
ts.compiled.dat2 <- read.csv("ts.compiled.dat.csv")
colnames(ts.compiled.dat2)

breaks <- c(-10, seq(-5,5, by=0.1), 10)
rownames(ts.compiled.dat2) <- ts.compiled.dat2[,9]
pheatmap(ts.compiled.dat2[,10:16],
         cluster_rows=F, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=10,
         fontsize_row=, 
         fontsize_col=10, 
         cellwidth=10, 
         cellheight=,
         cutree_rows = 2,
         cutree_cols = 4,
         #labels_row=rownames,
         main=)

colfunc2 <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))
breaks2 <- c(0, 0.1, .2, .5, 1, 2, 3)
rownames(ts.compiled.dat2) <- ts.compiled.dat2[,9]
pheatmap(ts.compiled.dat2[,18],
         cluster_rows=F, cluster_col=F,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks2, 
         color=colfunc2(length(breaks2)+1), 
         border_color='black', 
         fontsize=10,
         fontsize_row=, 
         fontsize_col=10, 
         cellwidth=10, 
         cellheight=,
         cutree_rows = 2,
         cutree_cols = 4,
         #labels_row=rownames,
         main=)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}



#### use this code for the z-score plotting with floated names 2023.05.08
##Merged dat 6, FFs removed, gene name added, use for plotting 
##3669.1_DAzscores
#plotzscores_dat <- read.csv("dat.merge6_noFFs.csv")
zTS3669.1 <- dat.merge7[ , c("TSQ.ID", "gene", "Allele", "3669-1.TS")]
zTS3669.1 <- na.omit(zTS3669.1)
names(zTS3669.1)[3:4] <- c("Name", "Zscore")
zTS3669.1 = zTS3669.1[order(zTS3669.1$Zscore, decreasing=TRUE),]
#remove FFs
zTS3669.1 = subset(zTS3669.1, !gene %in% removeffs)
removeffs <- c("YAL011W", "YGR002C", "YJL194W", "YOR181W", "YLR274W", "YMR168C", "YOR057W", "YGL130W", "YBR060C", "YPL094C", "YDR334W", "YGR140W", "YCR077C")
number1 <- c(1:(nrow(zTS3669.1)))
zTS3669.1_num <- cbind(number1, zTS3669.1)
increased <- nrow(subset(zTS3669.1_num, Zscore > 2))
decreased <- nrow(subset(zTS3669.1_num, Zscore < -2))
nochange <- nrow(subset(zTS3669.1_num, Zscore < 2 & Zscore > -2))
plotzTS3669.1<- ggplot(zTS3669.1_num, aes(x=number1, y=Zscore, label=Name))+
  geom_point(color = "grey", alpha=0.5)+
  geom_point(data = zTS3669.1_num %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = zTS3669.1_num %>% filter(Zscore <= -2), color = "darkred") +
  #geom_text(data = zTS3669.1_num %>% filter(Zscore >= 2), hjust = 0, nudge_x = 50, size = 2) +
  geom_text_repel(data = zTS3669.1_num %>% filter(Zscore >= 2), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf)+
  geom_text_repel(data = zTS3669.1_num %>% filter(Zscore <= -2.75), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzTS3669.1

##3669.2_DAzscores
zDA3669.2 <- plotzscores_dat[ , c("number", "gene", "name", "DA3669_2")]
names(zDA3669.2)[3:4] <- c("Name", "Zscore")
zDA3669.2 = zDA3669.2[order(zDA3669.2$Zscore, decreasing=TRUE),]
number1 <- c(1:(nrow(zDA3669.2)))
zDA3669.2_num <- cbind(number1, zDA3669.2)
increased <- nrow(subset(zDA3669.2, Zscore > 2))
decreased <- nrow(subset(zDA3669.2, Zscore < -2))
nochange <- nrow(subset(zDA3669.2, Zscore < 2 & Zscore > -2))
plotzDA3669.2<- ggplot(zDA3669.2_num, aes(x=number1, y=Zscore, label=Name))+
  geom_point(color = "grey", alpha=0.5)+
  geom_point(data = zDA3669.2_num %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = zDA3669.2_num %>% filter(Zscore <= -2), color = "darkred") +
  #geom_text(data = zDA3669.2_num %>% filter(Zscore >= 5), hjust = 0, nudge_x = 50, size = 2) +
  geom_text_repel(data = zDA3669.2_num %>% filter(Zscore >= 5), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf)+
  geom_text_repel(data = zDA3669.2_num %>% filter(Zscore <= -4.5), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzDA3669.2

##3669.1.NAA_DAzscores
zDA3669.1.NAA <- plotzscores_dat[ , c("number", "gene", "name", "DA3669_1_NAA")]
names(zDA3669.1.NAA)[3:4] <- c("Name", "Zscore")
zDA3669.1.NAA = zDA3669.1.NAA[order(zDA3669.1.NAA$Zscore, decreasing=TRUE),]

number1 <- c(1:(nrow(zDA3669.1.NAA)))
zDA3669.1.NAA_num <- cbind(number1, zDA3669.1.NAA)
increased <- nrow(subset(zDA3669.1.NAA_num, Zscore > 2))
decreased <- nrow(subset(zDA3669.1.NAA_num, Zscore < -2))
nochange <- nrow(subset(zDA3669.1.NAA_num, Zscore < 2 & Zscore > -2))

plotzDA3669.1.NAA<- ggplot(zDA3669.1.NAA_num, aes(x=number1, y=Zscore, label=Name))+
  geom_point(color = "grey", alpha=0.5)+
  geom_point(data = zDA3669.1.NAA_num %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = zDA3669.1.NAA_num %>% filter(Zscore <= -2), color = "darkred") +
  geom_text_repel(data = zDA3669.1.NAA_num %>% filter(Zscore >= 4), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf)+
  geom_text_repel(data = zDA3669.1.NAA_num %>% filter(Zscore <= -4), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzDA3669.1.NAA

##3669.2.NAA_DAzscores
zDA3669.2.NAA <- plotzscores_dat[ , c("number", "gene", "name", "DA3669_2_NAA")]
names(zDA3669.2.NAA)[3:4] <- c("Name", "Zscore")
zDA3669.2.NAA = zDA3669.2.NAA[order(zDA3669.2.NAA$Zscore, decreasing=TRUE),]

number1 <- c(1:(nrow(zDA3669.2.NAA)))
zDA3669.2.NAA_num <- cbind(number1, zDA3669.2.NAA)
increased <- nrow(subset(zDA3669.2.NAA_num, Zscore > 2))
decreased <- nrow(subset(zDA3669.2.NAA_num, Zscore < -2))
nochange <- nrow(subset(zDA3669.2.NAA_num, Zscore < 2 & Zscore > -2))

plotzDA3669.2.NAA<- ggplot(zDA3669.2.NAA_num, aes(x=number1, y=Zscore, label=Name))+
  geom_point(color = "grey", alpha=0.5)+
  geom_point(data = zDA3669.2.NAA_num %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = zDA3669.2.NAA_num %>% filter(Zscore <= -2), color = "darkred") +
  geom_text_repel(data = zDA3669.2.NAA_num %>% filter(Zscore >= 4.8), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf)+
  geom_text_repel(data = zDA3669.2.NAA_num %>% filter(Zscore <= -4.5), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzDA3669.2.NAA

##3670.1_DAzscores
zDA3670.1 <- plotzscores_dat[ , c("number", "gene", "name", "DA3670_1")]
names(zDA3670.1)[3:4] <- c("Name", "Zscore")
zDA3670.1 = zDA3670.1[order(zDA3670.1$Zscore, decreasing=TRUE),]

number1 <- c(1:(nrow(zDA3670.1)))
zDA3670.1_num <- cbind(number1, zDA3670.1)
increased <- nrow(subset(zDA3670.1_num, Zscore > 2))
decreased <- nrow(subset(zDA3670.1_num, Zscore < -2))
nochange <- nrow(subset(zDA3670.1_num, Zscore < 2 & Zscore > -2))

plotzDA3670.1_num<- ggplot(zDA3670.1_num, aes(x=number1, y=Zscore, label=Name))+
  geom_point(color = "grey", alpha=0.5)+
  geom_point(data = zDA3670.1_num %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = zDA3670.1_num %>% filter(Zscore <= -2), color = "darkred") +
  geom_text_repel(data = zDA3670.1_num %>% filter(Zscore >= 4.2), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf)+
  geom_text_repel(data = zDA3670.1_num %>% filter(Zscore <= -4), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzDA3670.1_num

##3670.2_DAzscores
zDA3670.2 <- plotzscores_dat[ , c("number", "gene", "name", "DA3670_2")]
names(zDA3670.2)[3:4] <- c("Name", "Zscore")
zDA3670.2 = zDA3670.2[order(zDA3670.2$Zscore, decreasing=TRUE),]

number1 <- c(1:(nrow(zDA3670.2)))
zDA3670.2_num <- cbind(number1, zDA3670.2)
increased <- nrow(subset(zDA3670.2_num, Zscore > 2))
decreased <- nrow(subset(zDA3670.2_num, Zscore < -2))
nochange <- nrow(subset(zDA3670.2_num, Zscore < 2 & Zscore > -2))

plotzDA3670.2_num<- ggplot(zDA3670.2_num, aes(x=number1, y=Zscore, label=Name))+
  geom_point(color = "grey", alpha=0.5)+
  geom_point(data = zDA3670.2_num %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = zDA3670.2_num %>% filter(Zscore <= -2), color = "darkred") +
  geom_text_repel(data = zDA3670.2_num %>% filter(Zscore >= 4.5), force = 0.5, nudge_x = 500, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf)+
  geom_text_repel(data = zDA3670.2_num %>% filter(Zscore <= -4), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =2, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzDA3670.2_num
