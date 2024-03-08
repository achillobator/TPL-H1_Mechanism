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
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

ORFmap <- read.csv(file = "ORF.csv", stringsAsFactors = F)
names(ORFmap) <- c("gene")
dat.3669.1.DA.ratiozscore <- read.csv(file = "3669-1_DAbyplate_Z.csv", stringsAsFactors = F)
trim.3669.1.DA <- dat.3669.1.DA.ratiozscore[ , c("gene", "Zscore")]
dat.merge <- merge(ORFmap, trim.3669.1.DA, by = "gene", all.y=T)
names(dat.merge) <- c("gene","3669-1.DA")

dat.3669.1.NAA.DA.ratiozscore <- read.csv(file = "3669-1_NAA_DAbyplate_Z.csv", stringsAsFactors = F)
trim.3669.1.NAA.DA <- dat.3669.1.NAA.DA.ratiozscore[ , c("gene", "Zscore")]
dat.merge2 <- merge(dat.merge, trim.3669.1.NAA.DA, by = "gene", all.y=T)
names(dat.merge2) <- c("gene","3669-1.DA","3669-1.NAA.DA")

dat.3669.2.DA.ratiozscore <- read.csv(file = "3669-2_DAbyplate_Z.csv", stringsAsFactors = F)
trim.3669.2.DA <- dat.3669.2.DA.ratiozscore[ , c("gene", "Zscore")]
dat.merge3 <- merge(dat.merge2, trim.3669.2.DA, by = "gene", all.y=T)
names(dat.merge3) <- c("gene","3669-1.DA","3669-1.NAA.DA","3669-2.DA")

dat.3669.2.NAA.DA.ratiozscore <- read.csv(file = "3669-2_NAA_DAbyplate_Z.csv", stringsAsFactors = F)
trim.3669.2.NAA.DA <- dat.3669.2.NAA.DA.ratiozscore[ , c("gene", "Zscore")]
dat.merge4 <- merge(dat.merge3, trim.3669.2.NAA.DA, by = "gene", all.y=T)
names(dat.merge4) <- c("gene","3669-1.DA","3669-1.NAA.DA","3669-2.DA","3669-2.NAA.DA")

dat.3670.1.DA.ratiozscore <- read.csv(file = "3670-1_DAbyplate_Z.csv", stringsAsFactors = F)
trim.3670.1.DA <- dat.3670.1.DA.ratiozscore[ , c("gene", "Zscore")]
dat.merge5 <- merge(dat.merge4, trim.3670.1.DA, by = "gene", all.y=T)
names(dat.merge5) <- c("gene","3669-1.DA","3669-1.NAA.DA","3669-2.DA","3669-2.NAA.DA","3670-1.DA")

dat.3670.2.DA.ratiozscore <- read.csv(file = "3670-2_DAbyplate_Z.csv", stringsAsFactors = F)
trim.3670.2.DA <- dat.3670.2.DA.ratiozscore[ , c("gene", "Zscore")]
dat.merge6 <- merge(dat.merge5, trim.3670.2.DA, by = "gene", all.y=T)
names(dat.merge6) <- c("gene","3669-1.DA","3669-1.NAA.DA","3669-2.DA","3669-2.NAA.DA","3670-1.DA","3670-2.DA")

### Plot heatmap

colfunc <- colorRampPalette(c("#007AF4","black", "yellow"))
breaks <- c(-10, seq(-5,5, by=0.1), 10)

hm1 <- pheatmap(dat.merge6[,2:7],
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
hm1

colfunc <- colorRampPalette(c("#007AF4","black", "yellow"))

plot(dat.merge6[,2:7], 
     pch = 10, cex = .1)


APEX.dat <- read.csv(file = "APEXinfomean.csv", stringsAsFactors = F)
APEX.dat <- APEX.dat[ , c("gene", "log2apex")]
dat.merge.apex <- merge(dat.merge6, APEX.dat, by = "gene", all.x=TRUE)

hm2 <- pheatmap(dat.merge.apex[,2:8],
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
hm2

dat.merge.apex2 <- merge(dat.merge6, APEX.dat, by = "gene", all.y=FALSE)
dat.merge.apex2_3.9 <- subset(dat.merge.apex2, log2apex > 3.9)

rownames(dat.merge.apex2_3.9) = dat.merge.apex2_3.9[,1]
breaks <- c(-10, seq(-2,2, by=0.1), 10)

hm3 <- pheatmap(dat.merge.apex2_3.9[,2:7],
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
hm3

hm4 <- pheatmap(dat.merge.apex2_3.9[,8])
hm4

library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
hm3 + hm4

####

ff.da.info <- as.matrix(read.csv("DA_frequent_fliers.txt", stringsAsFactors = F))
dat.merge.noFFs <- subset(dat.merge6, !(gene %in% ff.da.info))
dat.merge.noFFs = dat.merge.noFFs[!duplicated(dat.merge.noFFs$gene),1:7]
rownames(dat.merge.noFFs) = dat.merge.noFFs[,1]
breaks <- c(-10, seq(-5,5, by=0.1), 10)

pheatmap(dat.merge.noFFs[,2:7],
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

up.da.pooled.dat <- as.matrix(read.csv("DA_up_hits_noffs.txt", stringsAsFactors = F))
dat.hitsup.noFFs <- subset(dat.merge.noFFs, gene %in% up.da.pooled.dat)

rownames(dat.hitsup.noFFs) = dat.hitsup.noFFs[,1]
breaks <- c(-10, seq(-5,5, by=0.1), 10)


pheatmap(dat.hitsup.noFFs[,2:7],
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
         cutree_rows = 5,
         cutree_cols = 2,
         main=)

dat.hitsup.noFFs.noNA <- na.omit(dat.hitsup.noFFs)
rownames(dat.hitsup.noFFs.noNA) = dat.hitsup.noFFs.noNA[,1]
pheatmap(dat.hitsup.noFFs.noNA[,2:7],
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
         cutree_rows = 10,
         cutree_cols = 2,
         main=)

## more kmeans clustering
dat.clust <- na.omit(dat.merge6)
kmm = kmeans(dat.clust[,2:7],9,nstart = 50,iter.max = 15) #we keep number of iter.max=15 to ensure the algorithm converges and nstart=50 to #ensure that atleat 50 random sets are choosen  
kmm

k.max <- 15
data <- dat.clust
wss <- sapply(1:k.max, 
              function(k){kmeans(data[,2:7], k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
library(factoextra)
fviz_cluster(kmm, data[,2:7],
             #palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

pheatmap(dat.clust[,2:7],
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
         cutree_rows = 7,
         cutree_cols = 2,
         main=)

############## independant plotting of z-scores for the paper

plot(dat.merge6[,2:7], 
     pch = 10, cex = .1)

# Plot Z scores
dat.3669.1.DA.ratiozscore = dat.3669.1.DA.ratiozscore[order(dat.3669.1.DA.ratiozscore$Zscore, decreasing=TRUE),]
number1 <- c(1:(nrow(dat.3669.1.DA.ratiozscore)))
test1 <- cbind(number1, dat.3669.1.DA.ratiozscore)
increased <- nrow(subset(dat.3669.1.DA.ratiozscore, Zscore > 2))
decreased <- nrow(subset(dat.3669.1.DA.ratiozscore, Zscore < -2))
nochange <- nrow(subset(dat.3669.1.DA.ratiozscore, Zscore < 2 & Zscore > -2))

pdf(file = file.path("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/Quant_RSGA_Rsession/Compare_all_by_Zscore.ratio", paste("test1", "Zplot.pdf", sep = "")), width = 6, height = 6)
par(mfrow = c(1, 1))
with(dat.3669.1.DA.ratiozscore, {
  plot(number1,
       Zscore,
       xlab = "",
       ylab = "Venus Abundance (Z-Score)",
       pch = 20,
       col = c(rep("darkgreen", increased), rep("darkgrey", nochange), rep("darkred", decreased)),
       xaxt='n',
       cex.lab = 0.75,
       cex.axis = 0.75,
       mgp = c(2.5,1,0),
       las = 1
       
  )
  abline(h = 0)
  abline(h = 2, lty = "dashed")
  axis(2, at=2,labels=TRUE, mgp = c(2.5,1,0), cex.axis = 0.75, las = 1, tick = FALSE)
  abline(h = (-2), lty = "dashed")
  axis(2, at= (-2),labels=TRUE, mgp = c(2.5,1,0), cex.axis = 0.75, las = 1, tick = FALSE)
  legend("topright", c("Increased Venus (Z > 2)","Decreased Venus (Z < -2)"), text.col = "black", pt.cex = 1, cex = 0.65, bty = "n", pch = 20, col = c("darkgreen", "darkred"))
  
})
dev.off()
ggplot(dat.3669.1.DA.ratiozscore[,"Zscore"])
plotH1H5_DA_zscore<- ggplot(test1, aes(x=number1, y=Zscore, label=gene))+
  geom_point(alpha=0.5)+
  geom_point(data = test1 %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = test1 %>% filter(Zscore <= -2), color = "darkred") +
  geom_text(data = test1 %>% filter(Zscore >= 5), hjust = 0, nudge_x = 50, size = 2) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)") 
  
plotH1H5_DA_zscore
write.csv(dat.merge6, "dat.merge6_forprocessing.csv")

##Merged dat 6, FFs removed, gene name added, use for plotting 
##3669.1_DAzscores
plotzscores_dat <- read.csv("dat.merge6_noFFs.csv")
zDA3669.1 <- plotzscores_dat[ , c("number", "gene", "name", "DA3669_1")]
names(zDA3669.1)[3:4] <- c("Name", "Zscore")
zDA3669.1 = zDA3669.1[order(zDA3669.1$Zscore, decreasing=TRUE),]
number1 <- c(1:(nrow(zDA3669.1)))
zDA3669.1_num <- cbind(number1, zDA3669.1)
increased <- nrow(subset(zDA3669.1, Zscore > 2))
decreased <- nrow(subset(zDA3669.1, Zscore < -2))
nochange <- nrow(subset(zDA3669.1, Zscore < 2 & Zscore > -2))
plotzDA3669.1<- ggplot(zDA3669.1_num, aes(x=number1, y=Zscore, label=Name))+
  geom_point(color = "grey", alpha=0.5)+
  geom_point(data = zDA3669.1_num %>% filter(Zscore >= 2), color = "#FFC003") +
  geom_point(data = zDA3669.1_num %>% filter(Zscore <= -2), color = "darkred") +
  #geom_text(data = zDA3669.1_num %>% filter(Zscore >= 5), hjust = 0, nudge_x = 50, size = 2) +
  geom_text_repel(data = zDA3669.1_num %>% filter(Zscore >= 5), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf)+
  geom_text_repel(data = zDA3669.1_num %>% filter(Zscore <= -4.5), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzDA3669.1

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
  geom_text_repel(data = zDA3669.2_num %>% filter(Zscore >= 5), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf)+
  geom_text_repel(data = zDA3669.2_num %>% filter(Zscore <= -4.5), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf) +
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
  geom_text_repel(data = zDA3670.1_num %>% filter(Zscore >= 4.2), force = 0.5, nudge_x = 300, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf)+
  geom_text_repel(data = zDA3670.1_num %>% filter(Zscore <= -4), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf) +
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
  geom_text_repel(data = zDA3670.2_num %>% filter(Zscore >= 4.5), force = 0.5, nudge_x = 500, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf)+
  geom_text_repel(data = zDA3670.2_num %>% filter(Zscore <= -4), force = 0.5, nudge_x = -500, direction = "y", hjust = 0,segment.size = 0.2, size =3, max.overlaps = Inf) +
  theme_classic()+
  geom_hline(yintercept=2,linetype=2, color="gray")+
  geom_hline(yintercept=-2,linetype=2, color="gray")+
  geom_hline(yintercept=0,linetype=1, color="black")+
  xlab("") +
  ylab("Zscore (Venus/tdTomato)")
plotzDA3670.2_num


##### applying the final code from the TS experiment 
#Applying the pairs function from the pych package for the all by all comparison with correlations
plot(dat.merge6[,2:7], 
     pch = 10, cex = .1,alpha=0.3,
     main="DA Experimental Correlation",
     lower.panel = NULL)

pairs.panels(dat.merge6[,2:7],
             pch = 20, cex = 1, alpha = 0.5)

pairs(dat.merge6[,2:7],
      pch = 20, cex = .1, alpha = 0.3, 
      lower.panel = NULL)

######## Heatmap with cytometry validation data added in for supplemental figure

union70.dat <- as.matrix(read.csv(file = "DA_up_hits_noffs.txt", stringsAsFactors = F))
dat.union70 <- subset(dat.merge6, gene %in% union70.dat)
rownames(dat.union70) = dat.union70[,1]
pheatmap(dat.union70[,2:7],
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

library(tidyverse)
dat.union70 <- distinct(dat.union70, .keep_all = TRUE)

rownames(dat.union70) = dat.union70[,1]
pheatmap(dat.union70[,2:7],
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
         cutree_rows = 3,
         cutree_cols = 2, 
         main=)

#pull in cytometry validation files:
DAcytoval.H1.dat <- read.csv("H1_5dat_DA.csv")
DAcytoval.N188.dat <- read.csv("N188_dat_DA.csv")
write.csv(dat.union70, "dat.union70.csv")
dat.union <- read.csv("dat.union70.csv")
dat.union69.70.up2_apex <- merge(dat.union, APEX.dat, by = "gene", all.x=TRUE)
DAcytoval.H1.dat <- DAcytoval.H1.dat %>% 
  rename("gene" = "Gene")
DAcytoval.N188.dat <- DAcytoval.N188.dat %>% 
  rename("gene" = "Gene")
dat.union69.70.up2_cyto <- merge(dat.union69.70.up2_apex, DAcytoval.H1.dat, by = "gene", all.x=TRUE)
dat.union69.70.up2_cyto <- dat.union69.70.up2_cyto[ -c(10) ]
dat.union69.70.up2_cyto <- dat.union69.70.up2_cyto %>% 
  rename("H1-H5_Venus" = "Venus.Amedian",
         "H1-H5_tomato" =  "tdtomato3.Amedian",
         "H1-H5_Ratio" =  "Ratio")
dat.union69.70.up2_cyto2 <- merge(dat.union69.70.up2_cyto, DAcytoval.N188.dat, by = "gene", all.x=TRUE)
dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2[ -c(13) ]
dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2 %>% 
  rename("N188_Venus" = "Venus.Amedian",
         "N188_tomato" =  "tdtomato3.Amedian",
         "N188_Ratio" =  "Ratio")

#dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2[-46,]
#dat.union69.70.up2_cyto2 <- dat.union69.70.up2_cyto2[-54,]
#dat.union69.70.up2_cyto3 <- dat.union69.70.up2_cyto2[ -c(10) ]

colfunc <- colorRampPalette(c("#007AF4","black", "yellow"))
breaks <- c(-10, seq(-5,5, by=0.1), 10)
rownames(dat.union69.70.up2_cyto2) = dat.union69.70.up2_cyto3[,2]
HMT_A <- pheatmap(dat.union69.70.up2_cyto2[,3:8],
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
                  cutree_rows = 3,
                  cutree_cols = 2, 
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
library(colorRamp2)
#omit NAs
dat.union69.70.up2_cyto3 <- dat.union69.70.up2_cyto2 %>% drop_na(X3669.1.DA)
dat.union69.70.up2_cyto3 <- dat.union69.70.up2_cyto3[, c(1, 2, 3, 5, 4,6,7,8,9,10,11,12,13,14,15)]

col_fun = colorRamp2(c(c(-10,-5, 0,5, 10)), c("#007AF4","#007AF4","black", "yellow","yellow"))
f1 = colorRamp2(seq(min(dat.union69.70.up2_cyto3[,3:8]), max(dat.union69.70.up2_cyto3[,3:8]), length = 3), c("#007AF4","black", "yellow"))
rownames(dat.union69.70.up2_cyto3) = dat.union69.70.up2_cyto3[,2]
HT1<- Heatmap(dat.union69.70.up2_cyto3[,3:8], name = "Z-score",col = col_fun,
              cluster_rows = TRUE, 
              cluster_columns = FALSE, 
              row_km = 6,
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 4),
              column_names_gp = gpar(fontsize = 5),
              rect_gp = gpar(col = "black", lwd = 1))
HT1

col_fun2 = colorRamp2(c(c(3,3.5,4)), c("red","orange","yellow"))
HT2<- Heatmap(dat.union69.70.up2_cyto3[,9], name = "Log2 APEX",col = col_fun2,
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 4),
              column_names_gp = gpar(fontsize = 5),
              rect_gp = gpar(col = "black", lwd = 1))
HT2
HT1+HT2

col_fun3 = colorRamp2(c(c(0,100,200)), c("red","orange", "yellow"))
HT3<- Heatmap(dat.union69.70.up2_cyto3[,10:11], name = "Raw FL (H1-H5)",
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              row_names_side = "right",
              col = col_fun3,
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 5),
              rect_gp = gpar(col = "black", lwd = 1))
HT3
HT1+HT2+HT3

col_fun4 = colorRamp2(c(c(0,0.49,0.5,1,2,3)), c("grey","grey","red","orange", "gold","yellow"))
HT4<- Heatmap(dat.union69.70.up2_cyto3[,12], name = "H1-H5 Ratio",
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              row_names_side = "right",
              col = col_fun4,
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 5),
              rect_gp = gpar(col = "black", lwd = 1))
HT4
HT1+HT2+HT3+HT4

col_fun5 = colorRamp2(c(c(0,0.04,.041,.1,.3,.5)), c("grey","grey","red","orange","gold", "yellow"))
HT5<- Heatmap(dat.union69.70.up2_cyto3[,15], name = "N188 Ratio",
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              row_names_side = "right",
              col = col_fun5,
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 5),
              rect_gp = gpar(col = "black", lwd = 1),
              show_row_names = TRUE)
HT5
HT1+HT2+HT3+HT4+HT5
HT1+HT2+HT4+HT5
HT1
HT1+HT2+HT4+HT5+rowAnnotation(rn = anno_text(rownames(dat.union69.70.up2_cyto3)))


pdf(file = "/Users/alexanderrleydon/Documents/UW/Nemhauser_lab/R-SGA/Compare_all_by_Zscore.ratio/HT1245.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 8) # The height of the plot in inches
HT1+HT2+HT4+HT5+rowAnnotation(rn = anno_text(rownames(dat.union69.70.up2_cyto3)))
# Step 3: Run dev.off() to create the file!
dev.off()


Test1 <- Heatmap(dat.union69.70.up2_cyto3[,3:8], name = "Z-score",col = col_fun,
                 cluster_rows = TRUE, 
                 cluster_columns = FALSE, 
                 row_km = 6,
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 4),
                 column_names_gp = gpar(fontsize = 5),
                 rect_gp = gpar(col = "black", lwd = 1))
Test1+HT3

Test2<- Heatmap(dat.union69.70.up2_cyto3[,9], name = "Log2 APEX",col = col_fun2,
              cluster_rows = FALSE, 
              row_names_side = "right",row_names_gp = gpar(fontsize = 4),
              column_names_gp = gpar(fontsize = 5),
              rect_gp = gpar(col = "black", lwd = 1))
Test2+rowAnnotation(rn = anno_text(rownames(dat.union69.70.up2_cyto3)))
