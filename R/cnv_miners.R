
library(circlize)


################################################################################
###### CNV

setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-MINERS/CNVkit/cns")
list.files()
X108_88<-read.table(file = "108_88_T.call.cns", sep="\t", head=T)
X108_88<-X108_88[,-c(4,6,7,8,9,10)]
colnames(X108_88)<-c("chr", "start", "end", "value1")

circos.clear()
circos.initializeWithIdeogram()
circos.genomicTrack(X108_88,  track.height = 0.5,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})
title("108_88")

################################################################################################################################################################
# Set the working directory to the directory containing the cns files
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-MINERS/CNVkit/cns")

# Get a list of the cns files
cns_files <- list.files()

# Loop through each cns file
for (cns_file in cns_files) {
  # Read in the cns file
  X <- read.table(file = cns_file, sep = "\t", header = TRUE)
  
  # Remove unnecessary columns
  X <- X[, -c(4, 6:10)]
  
  # Rename columns
  colnames(X) <- c("chr", "start", "end", "value1")
  
  # Create the circos plot
  circos.clear()
  circos.initializeWithIdeogram()
  circos.genomicTrack(X, track.height = 0.5,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  title(gsub("_T.call.cns", "", cns_file))
  
  # Save the plot to a file with the same name as the cns file
  pdf(paste0(cns_file, ".pdf"))
  circos.clear()
  circos.initializeWithIdeogram()
  circos.genomicTrack(X, track.height = 0.5,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                           col = ifelse(value[[1]] > 0, "red", "green"), ...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  title(gsub("_T.call.cns", "", cns_file))
  dev.off()
}




################################################################################################################################################################
# install pdfjam package
# Install and load packages
install.packages(c("pdftools", "magick", "tiff"))
library(pdftools)
library(magick)
library(tiff)

# Merge PDFs
pdf_files <- list.files(pattern = "\\.pdf$")
pdf_combine(pdf_files, output = "merged.pdf")


########################################################################################################################################################################################
# PERCENTAGE OF GAINS OR LOSSES IN EACH CHROMOSOME OVERALL
library(magrittr)
library(ggforestplot)
library(ggplot2)
library(tidyverse)
library(dplyr)


# get list of all .cns files
cns_files <- list.files(pattern = "\\.cns$")


# read in all .cns files and merge into a single dataframe
merged_df <- do.call(rbind, lapply(cns_files, read.table, header = TRUE))

# Remove unnecessary columns
merged_df <- merged_df[, -c(4, 6:10)]
colnames(merged_df) <- c("chr", "start", "end", "value1")
merged_df$chr <- sub("chr", "", merged_df$chr)


head(merged_df, n= 40)

ggplot(data = merged_df, aes(x = factor(chr), y = value1, color = value1 > 0)) +
  theme_minimal()+
  scale_x_discrete(limits = unique(merged_df$chr)) +
  geom_jitter(alpha=0.5) +
  labs(x = "Chromosome", y = "log2")+
  scale_color_manual(values = c("#4ae153", "#e1534a"), labels = c("Loss", "Gain"))




gain_loss_CNV<-ggplot(data = merged_df, aes(x = factor(chr), y = value1, color = value1 > 0)) +
  theme_minimal() +
  scale_x_discrete(limits = unique(merged_df$chr)) +
  geom_jitter(alpha = 0.3) +
  labs(x = "Chromosome", y = "log2") +
  scale_color_manual(values = c("#4ae153", "#e1534a"), labels = c("Loss", "Gain")) +
  # Calculate the mean value of value1 for each chr and value1 category
  stat_summary(fun = mean, aes(group = factor(value1 > 0)), geom = "line", size = 1.25) +
  theme(legend.title = element_blank())+
  ggtitle("Gains and losses Radon-Miners")


filename <- "gain_loss_CNV.png"
filepath <- "C:/Users/bedmar/Desktop/IDIBAPS/RADON-MINERS/CNVkit/cns"

# Save the plot as a PNG
ggsave(file.path(filepath, filename), width = 10, height = 7, dpi = 300)
# Add separate geom_line layers for gains and losses
#geom_line(data = subset(merged_df, value1 > 0), aes(group = chr), color = "#e1534a", size = 1, alpha=0.5) +
#geom_line(data = subset(merged_df, value1 < 0), aes(group = chr), color = "#4ae153", size = 1, alpha=0.5)+


