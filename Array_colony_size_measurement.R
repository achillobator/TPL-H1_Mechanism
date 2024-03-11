source("http://bioconductor.org/biocLite.R")
BiocManager::install("EBImage")
install.packages("githubinstall")
install.packages("devtools")
install_github('omarwagih/gitter')
library(dplyr)
library(tools)
library(devtools)
library(brew)
library(lifecycle)
require(gitter)

# this is used to extract colony size for all array images in a given folder

image.path = paste0("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/3670-2_TS/3670-2_TS /3670-2_TS_colonysize", "/")

setwd(image.path)

setwd = image.path

image.list = dir(image.path)
image.list = image.list[grep("JPG", image.list)]
#image.list = image.list[!grep("dat|gridded", image.list)]


lapply(1:length(image.list), function(i){
	r = paste0(image.path, image.list[1])
	s = paste0(image.path, image.list[i])
	gitter.batch(s, r, plate.format=c(32,48), 'p')
})

# import colony size data and plot

colony.size.files = dir(image.path)
colony.size.files = colony.size.files[grep("dat", colony.size.files)]

colony.size.dat = lapply(1:length(colony.size.files), function(i){
  dat = read.table(paste0(image.path, colony.size.files[i]), head=F, skip=4, sep="\t")
})

names(colony.size.dat) = colony.size.files 
names(colony.size.dat) = 1:length(colony.size.files)

raw.col.size = lapply(colony.size.dat, function(i) i[,3])

pdf("boxplot_colonySize.pdf")
boxplot(raw.col.size[], las =1, notch=T, outline=F)
dev.off()

fig.legend = cbind(image.list, x.tick.label = 1:length(image.list))
write.csv(fig.legend, "fig_legend.csv", row.names=F)

boxplot(raw.col.size[], las =1, notch=T, outline=F, plot=F)

