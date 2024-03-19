#GO analysis graphing for R-SGA upregulated hits. This was used to make the GO analysis from figure 2. 

library(ggplot2)
setwd("~/Documents/UW/nemhauser_lab/Reporter_SGA/Analysis_Toronto/SGA_Analysis/Quant_RSGA_Rsession/GO_analysis")
dat <- read.csv("AllupRSGA_GO.csv", header=TRUE)
dat <- dat[order(Value),]

plot1 <- ggplot(dat, aes(x=CORRECTED_PVALUE, y=TERM))+
  geom_point(aes(size= NUM_LIST_ANNOTATIONS, color = Value, stroke = 1, alpha=0.8)) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80))+
  ylab('GO term (Cellular Component') + 
  xlab('p-Value') + 
  theme_classic(base_family = 'Arial Bold', base_size = 10)+
  theme(text = element_text(size = 5), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(size=.4, color="gray", linetype = "dashed")) +
  scale_color_gradient(low="cyan", high="blue")
plot1 

dat_trim <- read.csv("AllupRSGA_GO_trim.csv", header=TRUE)
dat_trim2 <- subset(dat_trim, NUM_LIST_ANNOTATIONS > 5)

plot2 <- ggplot(dat_trim2, aes(x=logfe, y=reorder(TERM,logfe)))+
  geom_point(aes(size= NUM_LIST_ANNOTATIONS, color = Value, stroke = 1, alpha=0.8)) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80))+
  ylab('GO term (Biological Process)') + 
  xlab('Log10 Fold Enrichment') + 
  theme_classic(base_size = 6)+
  theme(text = element_text(size = 6), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(size=.4, color="gray", linetype = "dashed")) +
  scale_color_gradient(low="cyan", high="blue")
plot2 
ggsave("Trim1_AllUP_RSGA_GO.pdf", width = 5, height = 7)
ggsave("Trim1_AllUP_RSGA_GO.png", width = 5, height = 7)
