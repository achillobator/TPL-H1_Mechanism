## Analysis Pipeline for SGA-FLEX Array RNR3 Abundance Screening on Typhoon - 40 hour diploid
## -------------------------------------------------------------------------------------------

## Created by Jason A. Hendry (Originally typhoon_DMAanalysis-1indvSet.R)
## Last updated 2014/11/12
## Adapted for FLEX Collection by Krystal Laframboise
## Last updated 2015/07/16
## Adapted for Auxin circuit R-SGA by Alex Leydon 2022/08/12

###### Chapter 1. Aggregate Size and Intensity Information
# ---------------------------------------------------------

### Size Data
# ------------

## 1. Read colony size files and append rows/rbind.
## 2. Add conditions and replicate vectors.
## 3. Cbind colony sizes with replicate and condition.
## 4. Read 1536 FLEX map and cbind. 
## 5. Name file and write to computer.

#Define a few variables:
spanAnswer <- as.numeric(readline(prompt = "Desired LOESS span: ")) # typically = 1
bioRepAnswer <- as.numeric(readline(prompt = "Biological Replicate: "))
Technical <- as.numeric(readline(prompt = "Technical Replicate: "))
SlowGrow <- as.character(readline(prompt = "Name below 600px file (BioRep_x_Techx_Below600px): "))
#This means name the file that has your below600px colonies

#Chage to the colony size directory for both file path and setwd
sizeFilePath <- dirname("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/3669-1_TS/3669-1_TS_2d/3669-1_TS_colony_size/")
setwd("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/3669-1_TS/3669-1_TS_2d/3669-1_TS_colony_size/")
colonySize <- do.call(rbind, lapply(list.files(pattern = ".dat$"), read.table, skip = 3, colClasses = c(rep("numeric", 4),rep("NULL", 1)), fill = TRUE, header = FALSE))
names(colonySize) <- c("Row", "Column", "Size", "Circularity")
colonySize <- within(colonySize, {
        BioRep <- bioRepAnswer
        PlateNum <- rep(c(1:4), each = 1536)})

## Integrate colony size data with ORF map for FLEX collection.
setwd("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/3669-1_TS/3669-1_TS_2d/")
ORFmap <- read.csv(file = "TS_quadruplicate.csv", stringsAsFactors = F)
annotSize <- cbind(ORFmap[ , c("gene", "TSQ.ID", "plate")], colonySize) # Heirichal order assumption is made here: Plate, Row, Column.
annotSize <- annotSize[ , c( "gene","TSQ.ID","PlateNum", "Row", "Column", "Size", "Circularity","BioRep")] # Clean, by ordering columns

## Name file and write to computer.
# - Below prepares rather elaborate (but clarifying) output file name concatentation
coreName <- substr(        x = paste("BioRep", bioRepAnswer, "aggregateSize.csv", sep = "_"), 
                           start = 1, 
                           stop = nchar(paste("BioRep", bioRepAnswer, "aggregateSize.csv", sep = "_")) - 4)
write.csv(annotSize, file.path(sizeFilePath, coreName), row.names = FALSE)

### Fluorescence Intensity Data
#-------------------------------

# Select GFP and RFP intensity files.
# First select GFP, then RFP from first technical replicate.
## ONLY TAKE THE LAST 5 BLOCKS

cat("Select GFP and then RFP files from first replicate.\n")
fileName <- file("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/3669-1_TS/3669-1_TS_2d/3669-1_TS_2d_data.csv")
pathName <- dirname("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/3669-1_TS/3669-1_TS_2d/")
FPraw <- read.csv(fileName, skip = 32, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1")
#RFPraw <- read.csv(fileName, skip = 32, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1")
GFP <- FPraw[ , c("Block", "Row", "Column", "F488.Median...B488", "F488.Mean...B488")]
RFP <- FPraw[ , c("F532.Median...B532", "F532.Mean...B532")]
names(GFP)[4:5] <- c("GFPmedian", "GFPmean")
names(RFP) <- c("RFPmedian", "RFPmean")

# Combine GFP and RFP information, assumes they have the same order.
finalintensityData <- cbind(GFP, RFP)

# Subset proper blocks
#if (Technical == 1) {
       # intensityData <- intensityData[1:7680,]
#} else {
    #    intensityData <- intensityData[7681:15360,]  
#}

# Add Plate, Replicate Information
#PlateNum <- rep(c(1:8), each = 1536)
#TechRep <- rep(Technical, each = 1536*5)
#BioRep <- bioRepAnswer
#finalintensityData <- cbind(TechRep, intensityData)
#finalintensityData <- cbind(BioRep, finalintensityData)
#finalintensityData <- cbind(PlateNum, finalintensityData)
names(finalintensityData)[1] <- c("PlateNum")
#this just changes block to plate number - prob should do this in a better way

# Aggregate with size information (annotSize) and clean
finalintensityData <- cbind(finalintensityData[ , c("GFPmedian", "GFPmean", "RFPmedian", "RFPmean")], annotSize)
finalintensityData <- finalintensityData[ , c("PlateNum", "Row", "Column", "gene","TSQ.ID", "GFPmean", "GFPmedian", "RFPmean", "RFPmedian", "Size", "Circularity")]

# Save.
write.csv(finalintensityData, file.path(pathName, paste(coreName, "Intensity.csv", sep = "")), row.names = FALSE)


###### Chapter 2. Filtering Low Quality Data
# -------------------------------------------
# 1. Filter border strains using Row/Column information
# 2. Filter By Size
# 3. Filter Undetectable
# 4. Filter variable RFP
# 5. Filter frequent fliers from Henny

# Filter Border, YOR202W and Blanks
intensityFBorder <- subset(finalintensityData, gene != "YMR271C") 
#intensityFBorder <- subset(intensityFBorder, ORF != "blank")
#intensityFBorder <- subset(intensityFBorder, ORF != "YOR202W")

# Filter Size
SizeFilter <- 600 
intensityFBorderSize <- subset(intensityFBorder, Size > 600)
intensityFBorderSizeRemoved <- subset(intensityFBorder, Size < 600)
write.csv(intensityFBorderSizeRemoved, file.path(pathName, paste(coreName, SlowGrow, sep = "")), row.names = FALSE)

# Filter 'Undetectable'
intensityFilt <- subset(intensityFBorderSize, GFPmean > 0 & GFPmedian > 0 & RFPmean > 0 & RFPmedian > 0)

# Write
write.csv(intensityFilt, file.path(pathName, paste(coreName, "IntensityFilter.csv", sep = "")), row.names = FALSE)


###### Chapter 3. Normalization
# ------------------------------------------
# 1. Compute log2 ratio (logRatio) and brightness (Brightness) for filtered data
# 2. Plot ratio vs. size and center such that mean is zero across all sizes
##   - noticed that larger colonies tended to have smaller ratios leading to biases in data

### Input:
#        - intensityFilt

# - Compute log base 2 GFP/RFP ratios and brightness
intensityFilt <- within(intensityFilt, { 
        meanlogRatio <- log2(GFPmean/RFPmean)
        meanBrightness <- log2(GFPmean*RFPmean)
        medlogRatio <- log2(GFPmedian/RFPmedian)
        medBrightness <- log2(GFPmedian*RFPmedian)  
})

# Define LOESS span
span <- spanAnswer

# LOESS normalization

# - apply normalization

intensityNorm <- within(intensityFilt, {
        meanLOESS <- predict(loess(meanlogRatio ~ meanBrightness, span = span))
        medLOESS <- predict(loess(medlogRatio ~ medBrightness, span = span))
        meanlogRatioNorm <- meanlogRatio - meanLOESS
        medlogRatioNorm <- medlogRatio - medLOESS
})


# - order
intensityNormord <- intensityNorm[with(intensityNorm, order(meanlogRatioNorm, decreasing = T)), ]

# LOESS Plot 
# Plots of LOESS Normalization
pdf(file = file.path(pathName, paste(coreName, "IntensityFilterLOESS_Rs.pdf", sep = "")), width = 10, height = 5)
par(mfrow = c(1, 2))
with(intensityNormord, {
        plot(        meanBrightness,
                     meanlogRatio,
                     xlab = "Brightness",
                     ylab = "Mean log(Venus/tdTomato)",
                     main = "Pre-normalization",
                     pch = 16,
                     col = "darkgrey"
        )
        abline(h = 0)
        lines(        sort(meanBrightness), meanLOESS[order(meanBrightness)], lwd = 1.5, col = "firebrick")
        preIQR <- round(IQR(meanlogRatio), 3)
        preMedian <- round(median(meanlogRatio), 3)
        legend(	"topright", c(paste("Median:", preMedian), paste("IQR:", preIQR)), bty = "n")
        legend(	"topleft", "LOESS", lwd = 1.5, col = "firebrick", bty = "n")
        plot(	meanBrightness,
              meanlogRatioNorm,
              xlab = "Brightness",
              ylab = "Mean log(Venus/tdTomato) Normalized",
              main = "Post-normalization",
              pch = 16,
              col = "darkgrey"
        )
        abline(h = 0)		
        postIQR <- round(IQR(meanlogRatioNorm), 3)
        postMedian <- round(median(meanlogRatioNorm), 3)
        legend(	"topright", c(paste("Median:", postMedian), paste("IQR:", postIQR)), bty = "n")		
})

dev.off()

# Name file and save.
write.csv(intensityNormord, file.path(pathName, paste(coreName, "IntensityFilterLOESS_Gal.csv", sep = "")), row.names = FALSE)

###graph by block to make sure it works:
library(ggplot2)
library(ggthemes)
library(nlme)
library(mgcv)
library(dplyr)
library(RColorBrewer)
library(viridis)
#BiocManager::install("bigPint")
library(bigPint)
library(dplyr)
library(plotly)
library(Hmisc)

plot1 <- ggplot(data=intensityNormord, aes(as.factor(PlateNum), x=as.factor(PlateNum), y=medlogRatioNorm,color=as.factor(PlateNum))) + 
  geom_boxplot(show.legend = FALSE)+
  theme_classic(base_family = 'Arial Bold', base_size = 10)+
  ylab('488/532 Ratio (AU)') + 
  xlab('uM NAA') + 
  #coord_cartesian(ylim = c(0, 15))+
  #scale_y_continuous(breaks=c(-1, 0, 1))+
  #theme_minimal()+
  scale_color_viridis(discrete = TRUE, option = "D")
plot1


##
# Calculate Z scores
mean.ratio = aggregate(intensityNormord$meanlogRatioNorm, by=list(intensityNormord$TSQ.ID), mean)
sd.ratio = aggregate(intensityNormord$meanlogRatioNorm, by=list(intensityNormord$TSQ.ID), sd)

mean.sd.ratio = cbind(mean.ratio, sd.ratio[,2])
colnames(mean.sd.ratio) = c("TSQ.ID", "colony.mean.ratio", "colony.sd.ratio")

df.info = intensityNormord[!duplicated(intensityNormord$TSQ.ID),1:5]

mean.ratio.withInfo = merge(df.info, mean.sd.ratio, by.x = "TSQ.ID", by.y = "TSQ.ID", all.x=T) 

Zscore.per.p = as.data.frame(do.call(rbind, lapply(1:4, function(i){
  df = subset(mean.ratio.withInfo, PlateNum == i)
  Zscore = scale(df$colony.mean.ratio)
  df.Zscore = cbind(df[,1], Zscore)
})))

colnames(Zscore.per.p) = c("TSQ.ID", "Zscore")
Zscore.per.p.all.dat = merge(mean.ratio.withInfo, Zscore.per.p, by = "TSQ.ID", all.y=T)
# Write File
write.csv(Zscore.per.p.all.dat, file = file.path(pathName, paste("3669-1_TA_TSQID", "Bavg_Z.csv", sep = "")), row.names = F)

Zscore.per.p.all.dat[,8] = as.numeric(Zscore.per.p.all.dat[,8])
Zscore.per.p.all.dat = Zscore.per.p.all.dat[order(Zscore.per.p.all.dat$Zscore, decreasing=TRUE),]
hits_up = subset(Zscore.per.p.all.dat, Zscore >= 2)
hits_up = hits_up[order(hits_up$Zscore, decreasing=T),]
write.csv(hits_up, file = file.path(pathName, paste(coreName, "Up_hits.csv", sep = "")), row.names = F)

hits_Down = subset(Zscore.per.p.all.dat, Zscore < -2)
hits_Down = hits_Down[order(hits_Down$Zscore, decreasing=T),]
write.csv(hits_Down, file = file.path(pathName, paste(coreName, "Down_hits.csv", sep = "")), row.names = F)

#original zscore - calculated accross all plates, tossed in favor of the by-plate zscore
#Zscore <- (intensityNormord$meanlogRatioNorm - mean(intensityNormord$meanlogRatioNorm))/sd(intensityNormord$meanlogRatioNorm)

# Write File
write.csv(intensityNormord, file = file.path(pathName, paste(coreName, "Bavg_Z.csv", sep = "")), row.names = F)

# Plot Z scores
number <- c(1:(nrow(Zscore.per.p.all.dat)))
Zscore.per.p.all.dat <- cbind(number, Zscore.per.p.all.dat)
increased <- nrow(subset(Zscore.per.p.all.dat, Zscore > 2))
decreased <- nrow(subset(Zscore.per.p.all.dat, Zscore < -2))
nochange <- nrow(subset(Zscore.per.p.all.dat, Zscore < 2 & Zscore > -2))

pdf(file = file.path(pathName, paste(coreName, "Zplot.pdf", sep = "")), width = 6, height = 6)
par(mfrow = c(1, 1))
with(Zscore.per.p.all.dat, {
  plot(number,
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
plot(Zscore.per.p.all.dat[,"Zscore"])


########### Make z-scores for GFP, RFP and colony size and add to the dataset

#GFP
mean.GFP = aggregate(intensityNormord$GFPmean, by=list(intensityNormord$gene), mean)
sd.GFP = aggregate(intensityNormord$GFPmean, by=list(intensityNormord$gene), sd)
mean.sd.GFP = cbind(mean.GFP, sd.GFP[,2])
colnames(mean.sd.GFP) = c("gene", "colony.mean.GFP", "colony.sd.GFP")
df.info.gfp = intensityNormord[!duplicated(intensityNormord$gene),1:5]
mean.GFP.withInfo = merge(df.info.gfp, mean.sd.GFP, by.x = "gene", by.y = "gene", all.x=T) 

Zscore.GFP.per.p = as.data.frame(do.call(rbind, lapply(1:14, function(i){
  df = subset(mean.GFP.withInfo, PlateNum == i)
  Zscore.GFP = scale(df$colony.mean.GFP)
  df.Zscore = cbind(df[,1], Zscore.GFP)
})))

colnames(Zscore.GFP.per.p) = c("gene", "Zscore.GFP")
Zscore.Ratio.GFP.per.p.all.dat = merge(Zscore.per.p.all.dat, Zscore.GFP.per.p, by = "gene", all.y=T)

## write the hits based only on GFP
Zscore.Ratio.GFP.per.p.all.dat[,9] = as.numeric(Zscore.Ratio.GFP.per.p.all.dat[,9])
Zscore.Ratio.GFP.per.p.all.dat = Zscore.Ratio.GFP.per.p.all.dat[order(Zscore.Ratio.GFP.per.p.all.dat$Zscore.GFP, decreasing=TRUE),]
GFP_hits_up = subset(Zscore.Ratio.GFP.per.p.all.dat, Zscore.GFP >= 2)
GFP_hits_up = GFP_hits_up[order(GFP_hits_up$Zscore.GFP, decreasing=T),]
write.csv(GFP_hits_up, file = file.path(pathName, paste(coreName, "GFP_Up_hits.csv", sep = "")), row.names = F)

GFP_hits_Down = subset(Zscore.Ratio.GFP.per.p.all.dat, Zscore.GFP < -2)
GFP_hits_Down = GFP_hits_Down[order(GFP_hits_Down$Zscore.GFP, decreasing=T),]
write.csv(hits_Down, file = file.path(pathName, paste(coreName, "Down_hits.csv", sep = "")), row.names = F)






##RFP
mean.RFP = aggregate(intensityNormord$RFPmean, by=list(intensityNormord$gene), mean)
sd.RFP = aggregate(intensityNormord$RFPmean, by=list(intensityNormord$gene), sd)
mean.sd.RFP = cbind(mean.RFP, sd.RFP[,2])
colnames(mean.sd.RFP) = c("gene", "colony.mean.RFP", "colony.sd.RFP")
df.info.rfp = intensityNormord[!duplicated(intensityNormord$gene),1:5]
mean.RFP.withInfo = merge(df.info.rfp, mean.sd.RFP, by.x = "gene", by.y = "gene", all.x=T) 

Zscore.RFP.per.p = as.data.frame(do.call(rbind, lapply(1:14, function(i){
  df = subset(mean.RFP.withInfo, PlateNum == i)
  Zscore.RFP = scale(df$colony.mean.RFP)
  df.Zscore = cbind(df[,1], Zscore.RFP)
})))

colnames(Zscore.RFP.per.p) = c("gene", "Zscore.RFP")
Zscore.Ratio.GFP.RFP.per.p.all.dat = merge(Zscore.Ratio.GFP.per.p.all.dat, Zscore.RFP.per.p, by = "gene", all.y=T)

##Colony size
mean.Size = aggregate(intensityNormord$Size, by=list(intensityNormord$gene), mean)
sd.Size = aggregate(intensityNormord$Size, by=list(intensityNormord$gene), sd)
mean.sd.Size = cbind(mean.Size, sd.Size[,2])
colnames(mean.sd.Size) = c("gene", "colony.mean.Size", "colony.sd.Size")
df.info.Size = intensityNormord[!duplicated(intensityNormord$gene),1:5]
mean.Size.withInfo = merge(df.info.Size, mean.sd.Size, by.x = "gene", by.y = "gene", all.x=T) 

Zscore.Size.per.p = as.data.frame(do.call(rbind, lapply(1:14, function(i){
  df = subset(mean.Size.withInfo, PlateNum == i)
  Zscore.Size = scale(df$colony.mean.Size)
  df.Zscore = cbind(df[,1], Zscore.Size)
})))

colnames(Zscore.Size.per.p) = c("gene", "Zscore.Size")
Zscore.Ratio.GFP.RFP.Size.per.p.all.dat = merge(Zscore.Ratio.GFP.RFP.per.p.all.dat, Zscore.Size.per.p, by = "gene", all.y=T)
####

head(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat)
class(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,8])

Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,12] = as.numeric(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,12])
Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,9] = as.numeric(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,9])
Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,10] = as.numeric(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,10])
Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,11] = as.numeric(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,11])

plot(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,"Zscore.Size"])
plot(data_subset2)
pheatmap(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,9:12],
         cluster_rows=T, 
         cluster_col=T,
         treeheight_row=,
         treeheight_col=,
         col=rev(brewer.pal(9,"RdBu")),
         border_color='black',
         fontsize=6, 
         cellwidth=10, 
         cellheight=, 
         main=)

colfunc <- colorRampPalette(c("#007AF4","black", "yellow"))
breaks <- seq(min(data_subset2[,2:5]),  max(data_subset2[,2:5], by=0.1))
breaks <- c(-10, seq(-5,5, by=0.1), 10)

pheatmap(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,9:12],
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

##clustering
hclust()
data_subset <- subset(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat, select = c("TSQ.ID", "Zscore","Zscore.GFP","Zscore.RFP","Zscore.Size"))
my_hclust_gene <- hclust(dist(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat[,9:12]), method = "complete")
plot(my_hclust_gene, labels = FALSE, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = "Clusters", xlab = "Gene", ylab = "Height")

write.csv(Zscore.Ratio.GFP.RFP.Size.per.p.all.dat, file = file.path(pathName, paste("3669-1_TS", "_multi-Zscores.csv", sep = "")), row.names = F)

dat_up = merge(hits_up, data_subset, by.x = "TSQ.ID", by.y = "TSQ.ID") 
rownames(dat_up) = dat_up[,1]

pheatmap(dat_up[,9:12],
         cluster_rows=T, 
         cluster_col=T,
         treeheight_row=,
         treeheight_col=,
         col=rev(brewer.pal(9,"RdBu")),
         border_color='black',
         fontsize=4, 
         cellwidth=10, 
         cellheight=, 
         main=)

pheatmap(dat_up[,9:12],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=12,
         #fontsize_row=6, 
         #fontsize_col=6, 
         cellwidth=10, 
         cellheight=, 
         main=)

GFP_dat_up = merge(GFP_hits_up, data_subset, by.x = "gene", by.y = "gene") 

rownames(GFP_dat_up) = GFP_dat_up[,1]
pheatmap(GFP_dat_up[,10:13],
         cluster_rows=T, cluster_col=T,  	
         treeheight_row=, treeheight_col=,
         breaks= breaks, 
         color=colfunc(length(breaks)+1), 
         border_color='black', 
         fontsize=4,
         #fontsize_row=6, 
         #fontsize_col=6, 
         cellwidth=10, 
         cellheight=, 
         main=)
