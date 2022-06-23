#Generating barplot of CRSK OTUs that match conserved nubbin OTUs

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(grid)

#Read in data
data <- read.csv(file = "CRSK_OTUs.txt", header=T, sep="\t")

#convert data frame from a "wide" format to a "long" format
melted_data = melt(data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
melted_data$Sample <- factor(melted_data$Sample,levels=unique(melted_data$Sample))

#Making a nice palette with distinctive colours
n <- 10
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Generating plot
barplot = ggplot(melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
    geom_bar(stat = "identity", colour = "black") + 
    theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1), 
    axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), 
    axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(values = col_vector) +
    labs(x = "", y = "Relative Abundance (%)", fill = "OTU")
