#GO analysis graphing for APEX2 proximity labeling. 
#GO analysis was performed using https://go.princeton.edu/
#Raw enrichment was calculated in excel and data was imported into R as follows below.

library(ggplot2)
setwd("~/Documents/UW/nemhauser_lab/APEX/Second_APEX2_run/Go_terms_apex_list")
### Rerunning GO BP to add enrichment
dat828_BPall <- read.csv("828_BP_all.csv")
dat828_BPall <- subset(dat828_BPall, enrichment > 2 )
plot_BPall <- ggplot(dat828_BPall, aes(x=enrichment, y=reorder(term_name, enrichment)))+
  geom_point(aes(size= intersection_size, color = negative_log10_of_adjusted_p_value, stroke = 1, alpha=0.8)) +
  #scale_x_continuous(limits = c(0, 80), breaks = c(0,20,40,60,80))+
  ylab('GO term (Cellular Component)') + 
  xlab('Enrichment)') + 
  theme_classic(base_size = 10)+
  theme(text = element_text(size = 10), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(size=.4, color="gray", linetype = "dashed")) +
  scale_color_gradient(low="cyan", high="blue")
plot_BPall

dat828_BPtrim <- read.csv("BP_trim.csv")
dat828_BPtrim <- subset(dat828_BPtrim, enrichment > 2 )
plot_BPtrim <- ggplot(dat828_BPtrim, aes(x=enrichment, y=reorder(term_name, enrichment)))+
  geom_point(aes(size= intersection_size, color = negative_log10_of_adjusted_p_value, stroke = 1, alpha=0.8)) +
  #scale_x_continuous(limits = c(0, 80), breaks = c(0,20,40,60,80))+
  ylab('GO term (Cellular Component)') + 
  xlab('Enrichment)') + 
  theme_classic(base_size = 5)+
  theme(text = element_text(size = 5), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(size=.4, color="gray", linetype = "dashed")) +
  scale_color_gradient(low="cyan", high="blue")
plot_BPtrim

pdf(file = "/Users/aleydon/Documents/UW/nemhauser_lab/APEX/Second_APEX2_run/Go_terms_apex_list/BP_ENR.pdf",   # The directory you want to save the file in
    width = 4.5, # The width of the plot in inches
    height = 2) # The height of the plot in inches
# Step 2: Create the plot with R code
plot_BPtrim
# Step 3: Run dev.off() to create the file!
dev.off()
