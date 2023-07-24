setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS")
getwd()

library(maftools)
library(RColorBrewer)
library(stringr) 
library(tibble)
library(mclust)
library(ggplot2)
library(data.table)
library(tidyr)
library(gridExtra)
library(pryr)
library(dplyr)
#HD<-files[-c(2,6,7,8,9,10,11,15,16,20,22,25,26)]
#ID<-files[-c(1,2,3,4,5,12,13,14,15,16,17,18,19,20,21,23,24,27)]
#LD<-files[-c(1,3,4,5,6,7,8,9,10,11,12,13,14,17,18,19,21,22,23,24,25,26,27)]

#laml_HD=merge_mafs(HD)
#laml_ID=merge_mafs(ID)
#laml_LD=merge_mafs(LD)

#save.image(file = "mafs.Rdata")

load("mafs.Rdata")
load("maf.Rdata") # this one contains each miners sample 
load("maf_pathogenic.Rdata")

vc_table <- table(laml@data$Variant_Classification)
vc_pct <- prop.table(vc_table) * 100
as.data.frame(vc_table)


getFields(laml)

#laml_ID@data$Protein_Change<-laml_ID@data$Protein_Change %>% replace_na("")
#laml_LD@data$Protein_Change<-laml_LD@data$Protein_Change %>% replace_na("")



laml@data <- laml@data %>%
  add_column(START = gsub("/.*$","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(END = gsub("^[^/]*/","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(POSITION = laml@data$Protein_position)

laml@data <- laml@data %>%
  add_column(Protein_Change = paste(laml@data$START,laml@data$POSITION,laml@data$END))

laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = "NA", "", laml@data$Protein_Change)

#laml@data<-laml@data[,-c(118,119,120)]

#laml@data <- laml@data[-which(laml@data$Protein_Change == ""), ]
#laml@data$Protein_Change=gsub(pattern = '.*^', replacement = 'p.', laml@data$Protein_Change)


genes<-c("TP53","KRAS","EGFR","STK11","PIK3CA","ERBB4","NTRK3","BRAF","ROS1","ALK","PTEN","MET","ERBB2","RET","NTRK1","NRG1","NTRK2","FGFR2","FGFR4","ERBB3","FGFR1","NRAS","FGFR3","MAPK1","ATM","HPRT")


file=""
for (i in genes){
  filename=paste0("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/variants/tumor/lollipop_plots_drivers/lollipop_", i)
  svg(paste(file,filename,".svg",sep=""), width = 12, height = 4)
  lollipopPlot(
    maf = laml,
    gene = i,
    AACol = 'Protein_Change',
    showMutationRate = TRUE)
  dev.off()
}



lollipopPlot(
  maf = laml_HD,
  gene = "STK11",
  AACol = 'Protein_Change',
  showMutationRate = TRUE)



lollipopPlot(
  maf = laml_ID,
  gene = "EGFR",
  AACol = 'Protein_Change',
  showMutationRate = TRUE)


#PATHWAYS
library(plotly)

OncogenicPathways(maf = laml)

head(laml@data)
samples<-unique(laml@data$Source_MAF)
HD<-samples[-c(2,6,7,8,9,10,11,15,16,20,22,25,26)]
ID<-samples[-c(1,2,3,4,5,12,13,14,15,16,17,18,19,20,21,23,24,27)]
LD<-samples[-c(1,3,4,5,6,7,8,9,10,11,12,13,14,17,18,19,21,22,23,24,25,26,27)]

laml_HD<-laml
laml_HD@data<-laml_HD@data[is.element(laml_HD@data$Source_MAF, HD),]

laml_ID<-laml
laml_ID@data<-laml_ID@data[is.element(laml_ID@data$Source_MAF, ID),]

laml_LD<-laml
laml_LD@data<-laml_LD@data[is.element(laml_LD@data$Source_MAF, LD),]

file=""
#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/TP53_pathogenic",".svg",sep=""), width = 12, height = 4)
library(ggplot2)
library(cowplot)

par(mar = c(5,5,2,5))
help(par)
par(mfcol= c(3,3))
plot(c(1,2))
library(cowplot)
library(irtoys)
library(gridGraphics)

PlotOncogenicPathways(maf = laml_HD, pathways = "TP53", showTumorSampleBarcodes = T); plot1 <- recordPlot()
PlotOncogenicPathways(maf = laml_ID, pathways = "TP53", showTumorSampleBarcodes = T); plot2 <- recordPlot()
PlotOncogenicPathways(maf = laml_LD, pathways = "TP53", showTumorSampleBarcodes = T); plot3 <- recordPlot()

plot_grid(plot1, plot2, plot3, ncol = 3)


PlotOncogenicPathways(maf = laml_HD, pathways = "RTK-RAS", showTumorSampleBarcodes = T); plot4 <- recordPlot()
PlotOncogenicPathways(maf = laml_ID, pathways = "RTK-RAS", showTumorSampleBarcodes = T); plot5 <- recordPlot()
PlotOncogenicPathways(maf = laml_LD, pathways = "RTK-RAS", showTumorSampleBarcodes = T); plot6 <- recordPlot()

plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,
          ncol = 3, nrow = 3)


PlotOncogenicPathways(maf = laml_HD, pathways = "NRF2", showTumorSampleBarcodes = T); plot1 <- recordPlot()
PlotOncogenicPathways(maf = laml_ID, pathways = "NRF2", showTumorSampleBarcodes = T); plot2 <- recordPlot()
PlotOncogenicPathways(maf = laml_LD, pathways = "NRF2", showTumorSampleBarcodes = T); plot3 <- recordPlot()

plot_grid(plot1, plot2, plot3, labels = c("High dose", "Intermediate dose", "Low dose"),
          ncol = 3, label_x = 0, label_y = 0,
          hjust = -0.5, vjust = -0.5)


samples_ordered=c("108_88", "119_90", "166_86", "219_86", "350_88", "372_89", "374_90", "442_87", "443_88", "457_89", "585_88", "616_90", "631_88", "639_90", "269_87", "278_90", "292_87", "300_88", "309_90", "340_87", "603_87", "631_90", "635_89", "118_86", "389_88", "419_87", "520_89")
OncogenicPathways(maf = laml)
#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/RTK-RAS_pathogenic",".svg",sep=""), width = 12, height = 6)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); a <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/WNT_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "WNT", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); b <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/NOTCH_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "NOTCH", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); c <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/Hippo_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "Hippo", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); d <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/PI3K_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "PI3K", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); e <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/Cell_Cycle_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "Cell_Cycle", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); f <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/TGF-Beta_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "TGF-Beta", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); g <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/NRF2_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "TP53", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); h <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/MYC_pathogenic",".svg",sep=""), width = 12, height = 4)
PlotOncogenicPathways(maf = laml, pathways = "MYC", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); i <- recordPlot()
#dev.off()

#svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/pathways/NRF2_pathogenic",".svg",sep=""), width = 12, height = 4)
##################################PlotOncogenicPathways(maf = laml, pathways = "NRF2", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); j <- recordPlot()
#dev.off()

par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(2, 2, 0, 0)) # move plot to the right and up

plot_grid(a, NULL, b,NULL,c, align = "hv",ncol = 1, rel_heights = c(49,0,28,0,27))


plot_grid(a,NULL, align = "hv",ncol = 1, rel_heights = c(1,0))
plot_grid(b,NULL,ncol = 1, rel_heights = c(28,0))


plot_grid(a,b,c,d, ncol = 1, nrow = 3)

#
svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/a",".svg",sep=""), width = 6, height = 10)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); a <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/b",".svg",sep=""), width = 6, height = 6)
PlotOncogenicPathways(maf = laml, pathways = "WNT", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); b <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/c",".svg",sep=""), width = 6, height = 6)
PlotOncogenicPathways(maf = laml, pathways = "NOTCH", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); c <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/d",".svg",sep=""), width = 6, height = 4.25)
PlotOncogenicPathways(maf = laml, pathways = "Hippo", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); d <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/e",".svg",sep=""), width = 6, height = 3.6)
PlotOncogenicPathways(maf = laml, pathways = "PI3K", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); e <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/f",".svg",sep=""), width = 6, height = 2.6)
PlotOncogenicPathways(maf = laml, pathways = "Cell_Cycle", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); f <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/g",".svg",sep=""), width = 6, height = 2.4)
PlotOncogenicPathways(maf = laml, pathways = "TGF-Beta", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); g <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/h",".svg",sep=""), width = 6, height = 2.25)
PlotOncogenicPathways(maf = laml, pathways = "TP53", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); h <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/pathway/i",".svg",sep=""), width = 6, height = 1.75)
PlotOncogenicPathways(maf = laml, pathways = "MYC", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); i <- recordPlot()
dev.off()

#####3as

















#Detecting cancer driver genes based on positional clustering
colnames(laml@data)[118]<-'AAChange'
colnames(laml@data)
getFields(laml)
head(laml@data)
head(laml@data$AAChange, n=40)
head(laml@data$Protein_Change)

laml@data$Protein_Change<-gsub("p.$", "", laml@data$Protein_Change)

laml.sig=oncodrive(maf=laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
?oncodrive()




laml.pfamHD=pfamDomains(maf = laml_HD, AACol = 'Protein_Change', top = 10, labelSize = 1.5)
laml.pfamID=pfamDomains(maf = laml_ID, AACol = 'Protein_Change', top = 10, labelSize = 1.5)
laml.pfamLD=pfamDomains(maf = laml_LD, AACol = 'Protein_Change', top = 10, labelSize = 1.5)

require("ggrepel")
HD<-ggplot(data=laml.pfamHD$domainSummary)+
  geom_point(aes(x=nMuts, y=nGenes, size=nGenes, alpha=0.9))+
  geom_text_repel(data = subset(laml.pfamHD$domainSummary[1:10,]), aes(x=nMuts, y=nGenes, label = DomainLabel))+
  theme_minimal()+
  ggtitle("High dose (N=14)")+
  theme(legend.position = "none")+
  xlab("# mutations")+
  ylab("# genes")

ID<-ggplot(data=laml.pfamID$domainSummary)+
  geom_point(aes(x=nMuts, y=nGenes, size=nGenes, alpha=0.9))+
  geom_text_repel(data = subset(laml.pfamID$domainSummary[1:10,]), aes(x=nMuts, y=nGenes, label = DomainLabel))+
  theme_minimal()+
  ggtitle("Intermediate dose (N=9)")+
  theme(legend.position = "none")+
  xlab("# mutations")+
  ylab("# genes")

LD<-ggplot(data=laml.pfamLD$domainSummary)+
  geom_point(aes(x=nMuts, y=nGenes, size=nGenes, alpha=0.9))+
  geom_text_repel(data = subset(laml.pfamLD$domainSummary[1:10,]), aes(x=nMuts, y=nGenes, label = DomainLabel))+
  theme_minimal()+
  ggtitle("Low dose (N=4)")+
  theme(legend.position = "none")+
  xlab("# mutations")+
  ylab("# genes")

grid.arrange(HD,ID,LD,nrow=1,ncol=3)



library(tidyverse)
library(ggforce)

ggplot(data=laml.pfamHD$domainSummary)+
  geom_point(aes(x=nMuts, y=nGenes, size=nGenes, alpha=0.9))+
  geom_text_repel(data = subset(laml.pfamHD$domainSummary[1:10,]), aes(x=nMuts, y=nGenes, label = DomainLabel))+
  geom_mark_ellipse(aes(x=nMuts, y=nGenes, fill="blue", alpha = 0.0001))+
  theme_minimal()+
  xlim(0,521)+
  ylim(0,155)+
  ggtitle("High dose")+
  theme(legend.position = "none")+
  xlab("# mutations")+
  ylab("# genes")


####################################TEST DATA
laml.mafs = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clins = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

lamls = read.maf(maf = laml.mafs, clinicalData = laml.clins)
getFields(lamls@data)
colnames(lamls@data$Protein_Change)

laml.sig=oncodrive(maf=lamls, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

save.image("maftools.Rdata")


laml.sig <- oncodrive(maf = laml_HD, AACol = 'Protein_Change', minMut = 5)
HD<-ggplot(data=laml.sig)+
  geom_point(aes(x=fract_muts_in_clusters, y=-log10(fdr), color = ifelse(fdr>0.1, 'red', 'blue'), size=2))+
  scale_colour_manual(labels = c("<0.1", ">0.1"), values=c('red', 'blue')) +
  geom_text_repel(data = subset(laml.sig[1:9,]), aes(x=fract_muts_in_clusters, y=-log10(fdr), label = Hugo_Symbol, color="blue"))+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("PlotOncodrive high dose (N=14)")+
  xlab("Fraction of variants within clusters")
  
laml.sig <- oncodrive(maf = laml_ID, AACol = 'Protein_Change', minMut = 5)
ID<-ggplot(data=laml.sig)+
  geom_point(aes(x=fract_muts_in_clusters, y=-log10(fdr), color = ifelse(fdr>0.2, 'red', 'blue'), size=2))+
  scale_colour_manual(labels = c("<0.1", ">0.1"), values=c('red', 'blue')) +
  geom_text_repel(data = subset(laml.sig[1:4,]), aes(x=fract_muts_in_clusters, y=-log10(fdr), label = Hugo_Symbol, color="blue"))+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("PlotOncodrive intermediate dose (N=9)")+
  xlab("Fraction of variants within clusters")

laml.sig <- oncodrive(maf = laml_LD, AACol = 'Protein_Change', minMut = 5)
LD<-ggplot(data=laml.sig)+
  geom_point(aes(x=fract_muts_in_clusters, y=-log10(fdr), color = 'red', size=2))+
  scale_colour_manual(labels = c("<0.1", ">0.1"), values=c('blue')) +
  geom_text_repel(data = subset(laml.sig[1:9,]), aes(x=fract_muts_in_clusters, y=-log10(fdr), label = Hugo_Symbol))+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("PlotOncodrive low dose (N=4)")+
  xlab("Fraction of variants within clusters")

grid.arrange(HD,ID,LD,nrow=1,ncol=3)





################################################################################

#RADON RATS
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea2.2")

load("maf.Rdata")
df<-getGeneSummary(x=laml)
df<-df %>% slice_max(df$total, n = 30)
df<-df[,-11]

df<-melt(df, id.vars = unique(c("Hugo_Symbol")))

head(df)f

ggplot(data=df,aes(x=value, y=reorder(Hugo_Symbol, value), fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.x = element_blank())+
  ylab("count")+
  xlab("M")

###################################################################################################################################################

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1")


load("maf_v2.Rdata")

######NORMALIZED VARIANTS PER CHROMOSOME######
unique(laml@data$Tumor_Sample_Barcode)

subset_low <- laml@data[grep("^D6-", laml@data$Tumor_Sample_Barcode), ]
subset_high <- laml@data[!grep("^D6-", laml@data$Tumor_Sample_Barcode), ]

#cumulative
df<-data.frame(
  Variant_Classification=rownames(as.data.frame(rowSums(table(subset_high$Variant_Classification, subset_high$Tumor_Sample_Barcode)))),
  High=as.data.frame(rowSums(table(subset_high$Variant_Classification, subset_high$Tumor_Sample_Barcode)))$`rowSums(table(subset_high$Variant_Classification, subset_high$Tumor_Sample_Barcode))`/32,
  Low=as.data.frame(rowSums(table(subset_low$Variant_Classification, subset_low$Tumor_Sample_Barcode)))$`rowSums(table(subset_low$Variant_Classification, subset_low$Tumor_Sample_Barcode))`/4
)
df1<-read.table(file = "variant_classification_cumulative", sep = "\t")
df1_prop<-read.table(file = "variant_classification_cumulative_prop", sep = "\t")

#cumulative
df1_long<-melt(df1)
df1_long$group <- factor(df1_long$group, levels = c("High", "Low"))

df1_long_prop<-melt(df1_prop)
df1_long_prop$group <- factor(df1_long_prop$group, levels = c("High", "Low"))

A<-ggplot(data=df1_long, aes(x=group, y=value, fill=variable))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  theme_minimal()+
  scale_fill_manual(values = c("#1f78b4", "#6a3d9a", "#feff99", "#d53e4f", "#33a02b", "#e11a1b", "#000000", "#ff7f00", "#f9d1c5"))+
  xlab("")+
  ylab("Variants relative to n of samples")+
  theme(legend.title = element_blank(), text=element_text(family="Calibri"), legend.position = "none")


B<-ggplot(data=df1_long_prop, aes(x=group, y=value, fill=variable))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  theme_minimal()+
  scale_fill_manual(values = c("#1f78b4", "#6a3d9a", "#feff99", "#d53e4f", "#33a02b", "#e11a1b", "#000000", "#ff7f00", "#f9d1c5"))+
  xlab("")+
  ylab("% Variants")+
  theme(legend.title = element_blank(), text=element_text(family="Calibri"), legend.position = "none")


C<-ggplot(data=df_long_leg, aes(x=group, y=value, fill=variable))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  theme_void()+
  scale_fill_manual(values = c("#1f78b4", "#6a3d9a", "#feff99", "#d53e4f", "#33a02b", "#e11a1b", "#000000", "#ff7f00", "#f9d1c5"))+
  xlab("")+
  ylab("% Variants")+
  theme(legend.title = element_blank(), text=element_text(family="Calibri"))

grid.arrange(A, B, C, ncol=3, layout_matrix = rbind(c(1,1,2,2,3)))

#group
subset1 <- laml@data[grep("^D1-", laml@data$Tumor_Sample_Barcode), ]
subset2 <- laml@data[grep("^D3-", laml@data$Tumor_Sample_Barcode), ]
subset4 <- laml@data[grep("^D6-", laml@data$Tumor_Sample_Barcode), ]
subset3 <- laml@data[grep("^D12-", laml@data$Tumor_Sample_Barcode), ]
subset5 <- laml@data[laml@data$Tumor_Sample_Barcode=="^Pr"]

unique(subset5$Tumor_Sample_Barcode)

df<-data.frame(
  Variant_Classification=rownames(as.data.frame(rowSums(table(subset1$Variant_Classification, subset1$Tumor_Sample_Barcode)))),
  D1=as.data.frame(rowSums(table(subset1$Variant_Classification, subset1$Tumor_Sample_Barcode)))$`rowSums(table(subset1$Variant_Classification, subset1$Tumor_Sample_Barcode))`/12,
  D3=as.data.frame(rowSums(table(subset2$Variant_Classification, subset2$Tumor_Sample_Barcode)))$`rowSums(table(subset2$Variant_Classification, subset2$Tumor_Sample_Barcode))`/7,
  D6=as.data.frame(rowSums(table(subset3$Variant_Classification, subset3$Tumor_Sample_Barcode)))$`rowSums(table(subset3$Variant_Classification, subset3$Tumor_Sample_Barcode))`/4,
  D12=as.data.frame(rowSums(table(subset4$Variant_Classification, subset4$Tumor_Sample_Barcode)))$`rowSums(table(subset4$Variant_Classification, subset4$Tumor_Sample_Barcode))`/9,
  Pr=as.data.frame(rowSums(table(subset5$Variant_Classification, subset5$Tumor_Sample_Barcode)))$`rowSums(table(subset5$Variant_Classification, subset5$Tumor_Sample_Barcode))`/4
)

df<-read.table(file = "variant_classification_group", sep = "\t")
df_prop<-read.table(file = "variant_classification_group_prop", sep = "\t")

df_leg <- data.frame(
  group = c("."),
  Frame_Shift_Del = c(0),
  Frame_Shift_Ins = c(0),
  In_Frame_Del = c(0),
  In_Frame_Ins = c(0),
  Missense_Mutation = c(0),
  Nonsense_Mutation = c(0),
  Nonstop_Mutation = c(0),
  Splice_Site = c(0),
  Translation_Start_Site = c(0)
)

df_long<-melt(df)
df_long$group <- factor(df_long$group, levels = c("D1", "D3", "D6", "D12", "Pr"))

df_long_prop<-melt(df_prop)
df_long_prop$group <- factor(df_long_prop$group, levels = c("D1", "D3", "D6", "D12", "Pr"))

df_long_leg<-melt(df_leg)
df_long_leg$group <- factor(df_long_leg$group, levels = c("."))


D<-ggplot(data=df_long, aes(x=group, y=value, fill=variable))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  theme_minimal()+
  scale_fill_manual(values = c("#1f78b4", "#6a3d9a", "#feff99", "#d53e4f", "#33a02b", "#e11a1b", "#000000", "#ff7f00", "#f9d1c5"))+
  xlab("")+
  ylab("Variants relative to n of samples")+
  theme(legend.title = element_blank(), text=element_text(family="Calibri"), legend.position = "none")


E<-ggplot(data=df_long_prop, aes(x=group, y=value, fill=variable))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  theme_minimal()+
  scale_fill_manual(values = c("#1f78b4", "#6a3d9a", "#feff99", "#d53e4f", "#33a02b", "#e11a1b", "#000000", "#ff7f00", "#f9d1c5"))+
  xlab("")+
  ylab("% Variants")+
  theme(legend.title = element_blank(), text=element_text(family="Calibri"), legend.position = "none")

grid.arrange(D, E, C, ncol=3, layout_matrix = rbind(c(1,1,2,2,3)))








variants<- data.frame(D1=nrow(subset1)/12,
                      D3=nrow(subset2)/7,
                      D6=nrow(subset3)/4,
                      D12=nrow(subset4)/9,
                      Pr=nrow(subset5)/4
                      )

df<-melt(data = variants)

ggplot(data=df, aes(x = variable, y=value))+
  geom_bar(stat="identity")



chrom=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y")
size<-data.frame(
  chromosome = c("chr1", "chr2", "chr4", "chr3", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chrX", "chrY", "chrMT"),
  position = c(282763074, 266435125, 184226339, 177699992, 173707219, 147991367, 145729302, 133307652, 122095297, 112626471, 90463843, 52716770, 114033958, 115493446, 111246239, 90668790, 90843779, 88201929, 62275575, 56205956, 159970021, 3310458, 16313)
)


subset_counts <- list()
subset_sizes <- list()

# Iterate over the subsets
for (i in 1:5) {
  subset <- get(paste0("subset", i))
  tt <- table(subset$Chromosome)[match(chrom, names(table(subset$Chromosome)))]
  ss <- size$position[match(names(tt), gsub("chr", "", size$chromosome))]
    subset_counts[[i]] <- tt
  subset_sizes[[i]] <- ss
}

subset_counts[[1]] <- head(subset_counts[[1]], -1)
subset_counts[[2]] <- head(subset_counts[[2]], -1)
subset_counts[[3]] <- head(subset_counts[[3]], -1)
subset_counts[[4]] <- head(subset_counts[[4]], -1)
subset_counts[[5]] <- head(subset_counts[[5]], -1)

subset_sizes[[1]] <- head(subset_sizes[[1]], -1)
subset_sizes[[2]] <- head(subset_sizes[[2]], -1)
subset_sizes[[3]] <- head(subset_sizes[[3]], -1)
subset_sizes[[4]] <- head(subset_sizes[[4]], -1)
subset_sizes[[5]] <- head(subset_sizes[[5]], -1)

D1<-as.data.frame((subset_counts[[1]]/subset_sizes[[1]])*1000000)
D3<-as.data.frame((subset_counts[[2]]/subset_sizes[[2]])*1000000)
D6<-as.data.frame((subset_counts[[3]]/subset_sizes[[3]])*1000000)
D12<-as.data.frame((subset_counts[[4]]/subset_sizes[[4]])*1000000)
Pr<-as.data.frame((subset_counts[[5]]/subset_sizes[[5]])*1000000)

merged_df <- data.frame(chr=D1$Var1,
                        D1=D1$Freq/12,
                        D3=D3$Freq/7,
                        D6=D6$Freq/4,
                        D12=D12$Freq/9,
                        Pr=Pr$Freq/4
)

long_df=melt(merged_df)

ggplot(data = long_df, aes(x=chr, y=value, fill=variable)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal()+
  scale_fill_viridis_d()+
  xlab("Chromosomes")+
  ylab("")+
  theme(legend.title = element_blank(), axis.ticks = element_blank())






E<-ggplot(data=df, aes(x=Var1, y=Freq, fill = Var2))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  theme_minimal()+
  ggtitle("Variant Classification by Rn Group")+
  scale_fill_manual(values = c("#1f78b4", "#6a3d9a", "#feff99", "#d53e4f", "#33a02b", "#e11a1b", "#000000", "#ff7f00", "#f9d1c5"))+
  xlab("")+
  ylab("variants relative to n of samples")+
  theme(legend.title = element_blank(), text=element_text(family="Calibri"), legend.position = "none")

################################################################################
################################################################################

##############################################



###
oncoplot(maf = laml, top = ngenes,showTumorSampleBarcodes = T, gene_mar = 15)
onco_genes<-c("Cnot6","Il17d","Kat6b","Cd244","Frat2","Olfr1330-ps1","Olr675","Selenoi","Olr1394","Zfp106","AABR07014089.1","AC127887.1","Gli3","Naf1","RT1-A1","Tspan15","AABR07026236.1","Abca14","Baz2a","Echdc3","LOC108348175","Mecom","Obscn","Olr1682","Olr653","Pcdhb12","RGD1562319","Stab2","ENSRNOT00000035992","RT1-CE7")

df<-laml@data[grep(paste(onco_genes, collapse = "|"), laml@data$Hugo_Symbol), ]
unique(df$Hugo_Symbol)


df2<-data.frame(Tumor_Sample_Barcode = laml@data$Tumor_Sample_Barcode,
               Variant_Classification = laml@data$Variant_Classification)
result2 <- table(df2$Tumor_Sample_Barcode, df2$Variant_Classification)
write.table(result2, file = "VARIANT_CLASSIFICATION", sep = "\t")


df3<-data.frame(Tumor_Sample_Barcode = laml@data$Tumor_Sample_Barcode,
                Variant_Type = laml@data$Variant_Type)
result3 <- table(df3$Tumor_Sample_Barcode, df3$Variant_Type)
write.table(result3, file = "VARIANT_TYPE", sep = "\t")




table(laml@data$Variant_Classification)

prop.table(table(laml@data$Variant_Type))*100





result_df <- as.data.frame.matrix(result)

# Calculate percentages for each Hugo_Symbol
result_per <- within(result_df, {
  Total <- rowSums(result_df)  # Calculate the total count for each Hugo_Symbol
  Percentages <- round(100 * (result_df / Total), 2)  # Calculate the percentages
})

# Print the updated table with percentages
print(result_per)

write.table(result_per,file = "top_oncogenes_variant_classification_table.txt", sep = "\t")






drivers<-c("Tp53","Kras","Egfr","Stk11","Pik3Ca","Erbb4","Ntrk3","Braf","Ros1","Alk","Pten","Met","Erbb2","Ret","Ntrk1","Nrg1","Ntrk2","Fgfr2","Fgfr4","Erbb3","Fgfr1","Nras","Fgfr3","Mapk1","Atm","Hprt")

laml_drivers <- laml@data[Hugo_Symbol %in% drivers]


unique(laml_drivers$Source_MAF)
table(laml_drivers$Variant_Type)

unique(laml_drivers$Protein_position)

laml_drivers <- laml_drivers %>%
  add_column(START = gsub("/.*$","",laml_drivers$Amino_acids))

laml_drivers <- laml_drivers %>%
  add_column(END = gsub("^[^/]*/","",laml_drivers$Amino_acids))

laml_drivers <- laml_drivers %>%
  add_column(POSITION = laml_drivers$Protein_position)

laml_drivers <- laml_drivers %>%
  add_column(Protein_Change = paste(laml_drivers$START,laml_drivers$POSITION,laml_drivers$END))

laml_drivers$Protein_Change=sub(pattern = ' ', "", laml_drivers$Protein_Change)
laml_drivers$Protein_Change=sub(pattern = ' ', "", laml_drivers$Protein_Change)
laml_drivers$Protein_Change=sub(pattern = "NA", "", laml_drivers$Protein_Change)

laml_drivers$Protein_Change<-paste0("p.",laml_drivers$Protein_Change)
laml_drivers$Protein_Change[laml_drivers$Protein_Change == "p."] <- ""
#laml_drivers$Protein_Change[laml_drivers$Protein_Change %in% c("-", "?")] <- ""

#laml_drivers$Protein_Change[grepl("[-?*]", laml_drivers$Protein_Change)] <- ""
#laml_drivers <- laml_drivers[complete.cases(laml_drivers$Protein_Change),]

head(laml_drivers)
laml_drivers$Protein_Change

data <- data.frame(Source_MAF = laml_drivers$Source_MAF,
                   Hugo_Symbol = laml_drivers$Hugo_Symbol,
                   Chromosome = laml_drivers$Chromosome,
                   Start_Position = laml_drivers$Start_Position,
                   End_Position = laml_drivers$End_Position,
                   Strand = laml_drivers$Strand,
                   Variant_Classification = laml_drivers$Variant_Classification,
                   Variant_Type = laml_drivers$Variant_Type,
                   Exon_Number = laml_drivers$Exon_Number,
                   Exon_Number = laml_drivers$Amino_acids,
                   Exon_Number = laml_drivers$Protein_position,
                   Protein_Change = laml_drivers$Protein_Change)

data

write.table(data,file = "radon_rats_drivers_protein_changes.txt")




laml@data$Source_MAF
laml@data$Hugo_Symbol <- toupper(laml@data$Hugo_Symbol)

par()
plotProtein(gene = "ALK", refSeqID = "NM_004304")

head(laml@data$Protein_Change)
grep("ALK", laml@data$Hugo_Symbol)
laml@data[c(6248,8598, 14123, 14124, 14125, 17683, 18365, 36047),]

lollipopPlot(maf = laml,
             gene = 'ALK',
             AACol = 'Protein_Change',
             showMutationRate = TRUE)

OncogenicPathways(maf = laml)







samples_ordered<-unique(laml@data$Source_MAF)
samples_ordered <- gsub("\\.maf", "", samples_ordered)

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_Cell_Cycle",".svg",sep=""), width = 10, height = 1.7)
PlotOncogenicPathways(maf = laml, pathways = "Cell_Cycle", showTumorSampleBarcodes = T, removeNonMutated = FALSE,  sampleOrder = samples_ordered); Cell_Cycle <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_Hippo",".svg",sep=""), width = 10, height = 4.25)# good
PlotOncogenicPathways(maf = laml, pathways = "Hippo", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); Hippo <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_MYC",".svg",sep=""), width = 10, height = 2.3)
PlotOncogenicPathways(maf = laml, pathways = "MYC", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); MYC <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_NOTCH",".svg",sep=""), width = 10, height = 6)# good
PlotOncogenicPathways(maf = laml, pathways = "NOTCH", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); NOTCH <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_NRF2",".svg",sep=""), width = 10, height = 1.6)
PlotOncogenicPathways(maf = laml, pathways = "NRF2", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); NRF2 <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_PI3K",".svg",sep=""), width = 10, height = 3.75)# good
PlotOncogenicPathways(maf = laml, pathways = "PI3K", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); PI3K <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_RTK_RAS",".svg",sep=""), width = 10, height = 8.5)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); RTK_RAS <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_TGF_Beta",".svg",sep=""), width = 10, height = 1.6)
PlotOncogenicPathways(maf = laml, pathways = "TGF-Beta", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); TGF_Beta <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_TP53",".svg",sep=""), width = 10, height = 1.425)
PlotOncogenicPathways(maf = laml, pathways = "TP53", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); TP53 <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/ordered_WNT",".svg",sep=""), width = 10, height = 5.75)# good
PlotOncogenicPathways(maf = laml, pathways = "WNT", showTumorSampleBarcodes = T, removeNonMutated = FALSE, sampleOrder = samples_ordered); WNT <- recordPlot()
dev.off()

PlotOncogenicPathways(maf = laml)














svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/Cell_Cycle_2",".svg",sep=""), width = 10, height = 1.7)
PlotOncogenicPathways(maf = laml, pathways = "Cell_Cycle", showTumorSampleBarcodes = T, removeNonMutated = FALSE); Cell_Cycle <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/Hippo_2",".svg",sep=""), width = 10, height = 4.25)# good
PlotOncogenicPathways(maf = laml, pathways = "Hippo", showTumorSampleBarcodes = T, removeNonMutated = FALSE); Hippo <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/MYC_2",".svg",sep=""), width = 10, height = 2.3)
PlotOncogenicPathways(maf = laml, pathways = "MYC", showTumorSampleBarcodes = T, removeNonMutated = FALSE); MYC <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/NOTCH_2",".svg",sep=""), width = 10, height = 6)# good
PlotOncogenicPathways(maf = laml, pathways = "NOTCH", showTumorSampleBarcodes = T, removeNonMutated = FALSE); NOTCH <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/NRF2_2",".svg",sep=""), width = 10, height = 1.6)
PlotOncogenicPathways(maf = laml, pathways = "NRF2", showTumorSampleBarcodes = T, removeNonMutated = FALSE); NRF2 <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/PI3K_2",".svg",sep=""), width = 10, height = 3.75)# good
PlotOncogenicPathways(maf = laml, pathways = "PI3K", showTumorSampleBarcodes = T, removeNonMutated = FALSE); PI3K <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/RTK_RAS_2",".svg",sep=""), width = 10, height = 8.5)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS", showTumorSampleBarcodes = T, removeNonMutated = FALSE); RTK_RAS <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/TGF_Beta_2",".svg",sep=""), width = 10, height = 1.6)
PlotOncogenicPathways(maf = laml, pathways = "TGF-Beta", showTumorSampleBarcodes = T, removeNonMutated = FALSE); TGF_Beta <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/TP53_2",".svg",sep=""), width = 10, height = 1.425)
PlotOncogenicPathways(maf = laml, pathways = "TP53", showTumorSampleBarcodes = T, removeNonMutated = FALSE); TP53 <- recordPlot()
dev.off()

svg(paste(file,"C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1/pathway/WNT_2",".svg",sep=""), width = 10, height = 5.75)# good
PlotOncogenicPathways(maf = laml, pathways = "WNT", showTumorSampleBarcodes = T, removeNonMutated = FALSE); WNT <- recordPlot()
dev.off()

###################################################################################################################################################
#RADON RATS CANCER GENES
library(tidyverse)

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered")

load("maf.Rdata")
unique(laml@data$Source_MAF)
df<-getGeneSummary(x=laml)

df<-df %>% slice_max(df$total, n = 30)
df[,1]
head(df)

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/dose")
####################HIGH
load("highmaf.Rdata")
top50<-c("Gnpnat1","Tf","Kdm1a","Mdc1","Kmt2d","Alb","Msh6","Zfhx3","Ncl","Fasn","Muc1","Serpina1","Chd4","Col4a2","Ager","Eef2","Cyp2a1","Xrcc6","Bag6","Birc6","Smarca4","Plcg1","Apoe","Arid1b","Fat1","Abce1","Tnc","Hif1a","Dicer1","Kmt2c","Ubr5")

frequencies <- as.data.frame(table(laml@data$Hugo_Symbol))
frequencies <- frequencies[order(-frequencies$Freq),]
result <- frequencies$Var1
high_dose_top30_genes<-head(frequencies$Var1, n=30)

df2<-getGeneSummary(x=laml)
df2$Hugo_Symbol

my_data <- df2 %>% filter( 
  Hugo_Symbol %in% top30_union
)

df2 <- subset(my_data, select = -c(total,MutatedSamples,AlteredSamples))

df2<-melt(df2, id.vars = unique(c("Hugo_Symbol")))
df2 <- df2 %>% 
  mutate(value = value / 28)
head(df2)

a<-ggplot(data=df2,aes(x=value, y=reorder(Hugo_Symbol, value), fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  geom_vline(xintercept=median(my_data$total/28), color="red", linetype="dashed")+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.y = element_blank())+
  scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+
  scale_y_discrete(limits = rev(top30_union))+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "#ff7f00", "#f9d2c5"))+
  xlab("#mutations / n")+
  ggtitle(label = "High dose (n=28)",
          subtitle = paste("Median:",round(median(my_data$total/28), digits=2), sep = " "))

####################LOW
load("lowmaf.Rdata")

frequencies <- as.data.frame(table(laml@data$Hugo_Symbol))
frequencies <- frequencies[order(-frequencies$Freq),]
result <- frequencies$Var1
low_dose_top30_genes<-head(frequencies$Var1, n=30)

top30_union<-union(high_dose_top30_genes,low_dose_top30_genes)

df2<-getGeneSummary(x=laml)

df2<-df2 %>% slice_max(df2$total, n = 30)

my_data <- df2 %>% filter( 
  Hugo_Symbol %in% top30_union
)
df2 <- subset(my_data, select = -c(total,MutatedSamples,AlteredSamples))

df2<-melt(df2, id.vars = unique(c("Hugo_Symbol")))
df2 <- df2 %>% 
  mutate(value = value / 8)
head(df2)

b<-ggplot(data=df2,aes(x=value, y=reorder(Hugo_Symbol, value), fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  geom_vline(xintercept=median(my_data$total/8), color="red", linetype="dashed")+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.y = element_blank())+
  scale_x_reverse(limits = c(2.3,0),breaks = seq(0, 2.3, by = 0.2))+
  scale_y_discrete(limits = rev(top30_union))+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "#ff7f00", "#f9d2c5"))+
  xlab("#mutations / n")+
  ggtitle(label = "Low dose (n=8)",
          subtitle = paste("Median:",round(median(my_data$total/8), digits=2), sep = " "))


#################################GROUP
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/groups")

load("D1maf.Rdata")
top50<-c("Gnpnat1","Tf","Kdm1a","Mdc1","Kmt2d","Alb","Msh6","Zfhx3","Ncl","Fasn","Muc1","Serpina1","Chd4","Col4a2","Ager","Eef2","Cyp2a1","Xrcc6","Bag6","Birc6","Smarca4","Plcg1","Apoe","Arid1b","Fat1","Abce1","Tnc","Hif1a","Dicer1","Kmt2c","Ubr5")

df<-getGeneSummary(x=laml)

my_data <- df %>% filter( 
  Hugo_Symbol %in% top50
)

df <- subset(df, select = -c(total,MutatedSamples,AlteredSamples))

df<-melt(df, id.vars = unique(c("Hugo_Symbol")))

df <- df %>% 
  mutate(value = value / 12)

head(df)

c<-ggplot(data=df,aes(x=value, y=reorder(Hugo_Symbol, value), fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  geom_vline(xintercept=median(my_data$total/12), color="red", linetype="dashed")+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.y = element_blank())+
  scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+
  scale_y_discrete(limits = rev(top50))+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "#ff7f00", "#f9d2c5"))+
  xlab("#mutations / n")+
  ggtitle(label = "Group D1 (n=12)",
          subtitle = paste("Median:",round(median(my_data$total/12), digits=2), sep = " "))

load("D12maf.Rdata")
df2<-getGeneSummary(x=laml)

my_data <- df2 %>% filter( 
  Hugo_Symbol %in% top50
)

df2 <- subset(my_data, select = -c(total,MutatedSamples,AlteredSamples))

df2<-melt(df2, id.vars = unique(c("Hugo_Symbol")))
df2 <- df2 %>% 
  mutate(value = value / 9)
head(df2)

d<-ggplot(data=df2,aes(x=value, y=reorder(Hugo_Symbol, value), fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  geom_vline(xintercept=median(my_data$total/9), color="red", linetype="dashed")+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.y = element_blank())+
  scale_x_reverse(limits = c(2.3,0),breaks = seq(0, 2.3, by = 0.2))+
  scale_y_discrete(limits = rev(top50))+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "#ff7f00", "#f9d2c5"))+
  xlab("#mutations / n")+
  ggtitle(label = "Group D12 (n=9)",
          subtitle = paste("Median:", round(median(my_data$total/9), digits=2), sep = " "))


load("D6maf.Rdata")
df2<-getGeneSummary(x=laml)

my_data <- df2 %>% filter( 
  Hugo_Symbol %in% top50
)
df2 <- subset(my_data, select = -c(total,MutatedSamples,AlteredSamples))

df2<-melt(df2, id.vars = unique(c("Hugo_Symbol")))
df2 <- df2 %>% 
  mutate(value = value / 4)
head(df2)

e<-ggplot(data=df2,aes(x=value, y=reorder(Hugo_Symbol, value), fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  geom_vline(xintercept=median(my_data$total/4), color="red", linetype="dashed")+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.y = element_blank())+
  scale_x_reverse(limits = c(2.3,0),breaks = seq(0, 2.3, by = 0.2))+
  scale_y_discrete(limits = rev(top50))+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "#ff7f00", "#f9d2c5"))+
  xlab("#mutations / n")+
  ggtitle(label = "Group D6 (n=4)",
          subtitle = paste("Median:", round(median(my_data$total/4),digits=2), sep = " "))

load("D3maf.Rdata")
df2<-getGeneSummary(x=laml)

my_data <- df2 %>% filter( 
  Hugo_Symbol %in% top50
)
df2 <- subset(my_data, select = -c(total,MutatedSamples,AlteredSamples))

df2<-melt(df2, id.vars = unique(c("Hugo_Symbol")))
df2 <- df2 %>% 
  mutate(value = value / 7)
head(df2)

f<-ggplot(data=df2,aes(x=value, y=reorder(Hugo_Symbol, value), fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  geom_vline(xintercept=median(my_data$total/7), color="red", linetype="dashed")+
  theme_minimal()+
  theme(legend.title=element_blank(),panel.grid.major.y = element_blank())+
  scale_x_reverse(limits = c(2.3,0),breaks = seq(0, 2.3, by = 0.2))+
  scale_y_discrete(limits = rev(top50))+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "#ff7f00", "#f9d2c5"))+
  xlab("#mutations / n")+
  ggtitle(label = "Group D3 (n=7)",
          subtitle = paste("Median:", round(median(my_data$total/7),digits=2), sep = " "))
f


#SAVE PLOTS
ggsave(file="C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/dose/lowvshigh.svg", plot=grid.arrange(b+theme(legend.position = "none")+ylab(element_blank()),a+theme(legend.position = "none")+ylab(element_blank()),ncol=2), width=10, height=6)


ggsave(file="C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/D12vsD1.svg", plot=grid.arrange(d+theme(legend.position = "none")+ylab(element_blank()),c+theme(legend.position = "none")+ylab(element_blank()),ncol=2), width=10, height=6)
ggsave(file="C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/D3vsD1.svg", plot=grid.arrange(e+theme(legend.position = "none")+ylab(element_blank()),c+theme(legend.position = "none")+ylab(element_blank()),ncol=2), width=10, height=6)
ggsave(file="C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/D6vsD1.svg", plot=grid.arrange(f+theme(legend.position = "none")+ylab(element_blank()),c+theme(legend.position = "none")+ylab(element_blank()),ncol=2), width=10, height=6)
ggsave(file="C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/D3vsD12.svg", plot=grid.arrange(e+theme(legend.position = "none")+ylab(element_blank()),d + scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+theme(legend.position = "none")+ylab(element_blank()),ncol=2), width=10, height=6)
ggsave(file="C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/D6vsD12.svg", plot=grid.arrange(f+theme(legend.position = "none")+ylab(element_blank()),d + scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+theme(legend.position = "none")+ylab(element_blank()),ncol=2), width=10, height=6)
ggsave(file="C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cancer_genes_filtered/plots/D6vsD3.svg", plot=grid.arrange(f+theme(legend.position = "none")+ylab(element_blank()),e + scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+theme(legend.position = "none")+ylab(element_blank()),ncol=2), width=10, height=6)



grid.arrange(b+theme(legend.position = "none")+ylab(element_blank()),a+theme(legend.position = "none")+ylab(element_blank()),ncol=2)#low VS high
grid.arrange(b+theme(legend.position = "none")+ylab(element_blank()),a+ylab(element_blank()),ncol=2)#low VS high

#D1

grid.arrange(d+theme(legend.position = "none")+ylab(element_blank()),c+theme(legend.position = "none")+ylab(element_blank()),ncol=2)#D12 VS D1
grid.arrange(e+theme(legend.position = "none")+ylab(element_blank()),c+theme(legend.position = "none")+ylab(element_blank()),ncol=2)#D3 VS D1
grid.arrange(f+theme(legend.position = "none")+ylab(element_blank()),c+theme(legend.position = "none")+ylab(element_blank()),ncol=2)#D6 VS D1

grid.arrange(e+theme(legend.position = "none")+ylab(element_blank()),d + scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+theme(legend.position = "none")+ylab(element_blank()),ncol=2)#D3 VS D12
grid.arrange(f+theme(legend.position = "none")+ylab(element_blank()),d + scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+theme(legend.position = "none")+ylab(element_blank()),ncol=2)#D6 VS D12

grid.arrange(f+theme(legend.position = "none")+ylab(element_blank()),e + scale_x_continuous(limits = c(0, 2.3),breaks = seq(0, 2.3, by = 0.2))+theme(legend.position = "none")+ylab(element_blank()),ncol=2)#D6 VS D3









################################################################################

#VARIANTS PLOT PROPORTIONS

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea2.1")
load("maf.Rdata")


ngenes=30
file=""
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top=ngenes,showBarcodes =T)

df <- subset(laml@data, select = c(Source_MAF,Variant_Classification))
df <- count(df, Source_MAF, Variant_Classification)

df$n[df$Source_MAF  == "D1.maf"] <- round(df$n[df$Source_MAF  == "D1.maf"]/12, digits= 2)
df$n[df$Source_MAF  == "D3.maf"] <- round(df$n[df$Source_MAF  == "D3.maf"]/7, digits= 2)
df$n[df$Source_MAF  == "D6.maf"] <- round(df$n[df$Source_MAF  == "D6.maf"]/4, digits= 2)
df$n[df$Source_MAF  == "D12.maf"] <- round(df$n[df$Source_MAF  == "D12.maf"]/9, digits= 2)
df$n[df$Source_MAF  == "Pr.maf"] <- round(df$n[df$Source_MAF  == "Pr.maf"]/4, digits= 2)

nrow(laml@data)
df$n[df$Source_MAF  == "D1.maf"] <- round(df$n[df$Source_MAF  == "D1.maf"]/28427, digits= 4)
df$n[df$Source_MAF  == "D3.maf"] <- round(df$n[df$Source_MAF  == "D3.maf"]/28427, digits= 4)
df$n[df$Source_MAF  == "D6.maf"] <- round(df$n[df$Source_MAF  == "D6.maf"]/28427, digits= 4)
df$n[df$Source_MAF  == "D12.maf"] <- round(df$n[df$Source_MAF  == "D12.maf"]/28427, digits= 4)
df$n[df$Source_MAF  == "Pr.maf"] <- round(df$n[df$Source_MAF  == "Pr.maf"]/28427, digits= 4)

metrics<-aggregate(df$n, by=list(Category=df$Source_MAF), FUN=sum)

ggplot(data=df,aes(x= factor(Source_MAF, level = c('D1.maf', 'D3.maf', 'D6.maf', 'D12.maf', 'Pr.maf')), y=n, fill=Variant_Classification))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "black", "#ff7f00", "#f9d2c5"))+
  xlab(element_blank())+
  scale_x_discrete(labels=c('D1','D3','D6','D12','Pr'))+
  ggtitle(label = "",
          subtitle = paste("Median:", round(median(metrics$x),digits=2), sep = " "))


#high vs low

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea2.2")
load("maf.Rdata")

df <- subset(laml@data, select = c(Source_MAF,Variant_Classification))
df<-count(df, Source_MAF, Variant_Classification)

nrow(laml@data)
df$n[df$Source_MAF  == "high.maf"] <- round(df$n[df$Source_MAF  == "high.maf"]/28, digits= 2)
df$n[df$Source_MAF  == "low.maf"] <- round(df$n[df$Source_MAF  == "low.maf"]/8, digits= 2)

df$n[df$Source_MAF  == "high.maf"] <- round(df$n[df$Source_MAF  == "high.maf"]/23158, digits= 4)
df$n[df$Source_MAF  == "low.maf"] <- round(df$n[df$Source_MAF  == "low.maf"]/23158, digits= 4)
metrics<-aggregate(df$n, by=list(Category=df$Source_MAF), FUN=sum)
df

ggplot(data=df,aes(x=factor(Source_MAF, level = c('high.maf', 'low.maf')), y=n, fill=Variant_Classification))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "black", "#ff7f00", "#f9d2c5"))+
  xlab(element_blank())+
  scale_x_discrete(labels=c('high', 'low'))+
  ggtitle(label = "",
          subtitle = paste("Median:", round(median(metrics$x),digits=2), sep = " "))


ggplot(data=df,aes(x= factor(Source_MAF, level = c('high.maf', 'low.maf')), y=n, fill=Variant_Classification))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "black", "#ff7f00", "#f9d2c5"))+
  xlab(element_blank())+
  scale_x_discrete(labels=c('D1','D3','D6','D12','Pr'))+
  ggtitle(label = "",
          subtitle = paste("Median:", round(median(metrics$x),digits=2), sep = " "))


#long vs short
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea3")
load("maf.Rdata")

df <- subset(laml@data, select = c(Source_MAF,Variant_Classification))
df<-count(df, Source_MAF, Variant_Classification)

df$n[df$Source_MAF  == "long.maf"] <- round(df$n[df$Source_MAF  == "long.maf"]/9, digits= 2)
df$n[df$Source_MAF  == "short.maf"] <- round(df$n[df$Source_MAF  == "short.maf"]/23, digits= 2)
df$n[df$Source_MAF  == "D6.maf"] <- round(df$n[df$Source_MAF  == "D6.maf"]/4, digits= 2)

metrics<-aggregate(df$n, by=list(Category=df$Source_MAF), FUN=sum)

ggplot(data=df,aes(x= factor(Source_MAF, level = c('long.maf', 'short.maf', 'D6.maf')), y=n, fill=Variant_Classification))+
  geom_bar(stat="identity", position = "stack")+
  geom_hline(yintercept = median(metrics$x), color="red", linetype="dashed")+
  theme_minimal()+
  scale_fill_manual(values=c("#1f78b4", "#6a3d9a", "#ffff99", "#d53e4f", "#33a02c", "#e31a1c", "black", "#ff7f00", "#f9d2c5"))+
  xlab(element_blank())+
  scale_x_discrete(labels=c('long', 'short', 'D6'))+
  ggtitle(label = "",
          subtitle = paste("Median:", round(median(metrics$x),digits=2), sep = " "))




  




df<-read.csv(file="radon-rats_variant classification.CSV",sep=",",header=T)
df$sample <- factor(df$sample, levels = unique(df$sample))
df$group <- factor(df$group, levels = unique(df$group))

DEL<-ggplot(data=df, aes(x=sample, y=DEL, fill=group))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.title = element_blank(), text=element_text(family="Calibri"))+
  ylab("Variants")+
  ggtitle("DEL")+ 
  scale_fill_viridis_d()+
  xlab(element_blank())

DNP<-ggplot(data=df, aes(x=sample, y=DNP, fill=group))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.title = element_blank(), text=element_text(family="Calibri"))+
  ylab("")+
  ggtitle("DNP")+
  scale_fill_viridis_d()+
  xlab(element_blank())

INS<-ggplot(data=df, aes(x=sample, y=INS, fill=group))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.title = element_blank(), text=element_text(family="Calibri"))+
  ylab("")+
  ggtitle("INS")+
  scale_fill_viridis_d()+
  xlab(element_blank())

ONP<-ggplot(data=df, aes(x=sample, y=ONP, fill=group))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.title = element_blank(), text=element_text(family="Calibri"))+
  ylab("Variants")+
  ggtitle("ONP")+
  scale_fill_viridis_d()+
  xlab(element_blank())

SNP<-ggplot(data=df, aes(x=sample, y=SNP, fill=group))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.title = element_blank(), text=element_text(family="Calibri"))+
  ylab("")+
  ggtitle("SNP")+
  scale_fill_viridis_d()+
  xlab(element_blank())

TNP<-ggplot(data=df, aes(x=sample, y=TNP, fill=group))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank(), legend.title = element_blank(), text=element_text(family="Calibri"))+
  ylab("")+
  ggtitle("TNP")+
  scale_fill_viridis_d()+
  xlab(element_blank())

DEL <- DEL + theme(legend.position = "none")
DNP <- DNP + theme(legend.position = "none")
INS <- INS + theme(legend.position = "none")
ONP <- ONP + theme(legend.position = "none")
SNP <- SNP + theme(legend.position = "none")
TNP <- TNP + theme(legend.position = "none")


grid.arrange(DEL, DNP, INS, ONP, SNP, TNP, ncol=3)


################################################################################
# Define the sample sizes and DEL values for each group
# Create a dataframe with the given values
df <- data.frame(
  group = c("D1", "D3", "D6", "D12", "Pr"),
  DEL = c(2870, 1493, 850, 1719, 433),
  DNP = c(496, 432, 219, 430, 96),
  INS = c(512, 578, 279, 545, 155),
  ONP = c(2, 1, 0, 2, 0),
  SNP = c(8676, 5698, 3229, 6137, 1564),
  TNP = c(24, 20, 9, 23, 6)
)

y

# Perform ANOVA test
model <- lm(DEL + DNP + INS + ONP + SNP + TNP  ~ group, data = df, weights = n)
summary(model)


df<-data.frame(
  group = laml@data$Source_MAF,
  type = laml@data$Variant_Type
)


df

model <- lm(type  ~ group, data = df)
anova(model)















setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea2.1")
load("maf.Rdata")
by(laml@data$Variant_Type, laml@data$Source_MAF, function(x) table(x))
nrow(laml@data)


setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea2.2")
load("maf.Rdata")
by(laml@data$Variant_Type, laml@data$Source_MAF, function(x) table(x))





# Create a contingency table of counts
counts <- table(laml@data$Variant_Type, laml@data$Source_MAF)

chisq.test(table)





# create a data frame with the counts and group sizes 
counts <- data.frame(
  DEL = c(2482, 1435, 1232, 793, 378),
  DNP = c(369, 275, 281, 190, 76),
  INS = c(280, 285, 322, 225, 111),
  ONP = c(2, 2, 1, 0, 0),
  SNP = c(6991, 4406, 4064, 2905, 1261),
  TNP = c(16, 15, 17, 8, 5),
  group = factor(rep(c("D1", "D12", "D3", "D6", "Pr"), c(12, 9, 7, 4, 4)))
)

# check the normality assumption
library(ggplot2)
ggplot(counts, aes(x = DEL)) + geom_density()  # check one group at a time

# check the equal variance assumption
bartlett.test(counts[, -6], counts$group)

# conduct ANOVA
fit <- aov(DEL ~ group, data = counts)
summary(fit)





































#####################################################################################################################



setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/drivers")

files=dir(pattern=".maf$", full.names = T)
laml=merge_mafs(files)

laml@data
unique(laml@data$Tumor_Sample_Barcode)
table(laml@data$Variant_Type)

unique(laml@data$Protein_position)

laml@data <- laml@data %>%
  add_column(START = gsub("/.*$","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(END = gsub("^[^/]*/","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(POSITION = laml@data$Protein_position)

laml@data <- laml@data %>%
  add_column(Protein_Change = paste(laml@data$START,laml@data$POSITION,laml@data$END))

laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = "NA", "", laml@data$Protein_Change)

laml@data$Protein_Change<-paste0("p.",laml@data$Protein_Change)
laml@data$Protein_Change[laml@data$Protein_Change == "p."] <- ""
laml@data$Protein_Change[laml@data$Protein_Change %in% c("-", "?")] <- ""

laml@data$Protein_Change[grepl("[-?*]", laml@data$Protein_Change)] <- ""
laml@data <- laml@data[complete.cases(laml@data$Protein_Change),]


laml@data


laml@data <- laml@data %>%
  add_column(CDS_Change = paste0("c.",laml@data$CDS_position,laml@data$Reference_Allele,">",laml@data$Tumor_Seq_Allele2))

laml@data


df<-laml@data(
  
)

df <- data.frame(Tumor_Sample_Barcode = laml@data$Tumor_Sample_Barcode,
                 Hugo_Symbol = laml@data$Hugo_Symbol,
                 Chromosome = laml@data$Chromosome,
                 Start_Position = laml@data$Start_Position,
                 End_Position = laml@data$End_Position,
                 Strand = laml@data$Strand,
                 Consequence = laml@data$Consequence,
                 Protein_Change = laml@data$Protein_Change,
                 Variant_Type = laml@data$Variant_Type,
                 CDS_Change = laml@data$CDS_Change)

df


write.table(df, file = "table.txt", sep = "\t", row.names = FALSE)


#####################################################################################################################

files=dir(pattern=".maf$", full.names = T)
laml=merge_mafs(files)

laml@data
unique(laml@data$Tumor_Sample_Barcode)
table(laml@data$Variant_Type)

unique(laml@data$Protein_position)

laml@data <- laml@data %>%
  add_column(START = gsub("/.*$","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(END = gsub("^[^/]*/","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(POSITION = laml@data$Protein_position)

laml@data <- laml@data %>%
  add_column(Protein_Change = paste(laml@data$START,laml@data$POSITION,laml@data$END))

laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = "NA", "", laml@data$Protein_Change)

laml@data$Protein_Change<-paste0("p.",laml@data$Protein_Change)
laml@data$Protein_Change[laml@data$Protein_Change == "p."] <- ""
laml@data$Protein_Change[laml@data$Protein_Change %in% c("-", "?")] <- ""

laml@data$Protein_Change[grepl("[-?*]", laml@data$Protein_Change)] <- ""
laml@data <- laml@data[complete.cases(laml@data$Protein_Change),]


laml@data


laml@data <- laml@data %>%
  add_column(CDS_Change = paste0("c.",laml@data$CDS_position,laml@data$Reference_Allele,">",laml@data$Tumor_Seq_Allele2))

laml@data


df<-laml@data(
  
)

df <- data.frame(Tumor_Sample_Barcode = laml@data$Tumor_Sample_Barcode,
                 Hugo_Symbol = laml@data$Hugo_Symbol,
                 Chromosome = laml@data$Chromosome,
                 Start_Position = laml@data$Start_Position,
                 End_Position = laml@data$End_Position,
                 Strand = laml@data$Strand,
                 Consequence = laml@data$Consequence,
                 Protein_Change = laml@data$Protein_Change,
                 Variant_Type = laml@data$Variant_Type,
                 CDS_Change = laml@data$CDS_Change)

df

write.table(df, file = "table_pathogenic.txt", sep = "\t", row.names = FALSE)












#####################################################################################################################



setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/drivers_orthologs")

files=dir(pattern=".maf$", full.names = T)
#Silent variants: 467 

laml=merge_mafs(files)

unique(laml@data$IMPACT)

nrow(laml@data)
unique(laml@data$Tumor_Sample_Barcode)
table(laml@data$Variant_Type)

unique(laml@data$Protein_position)

laml@data <- laml@data %>%
  add_column(START = gsub("/.*$","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(END = gsub("^[^/]*/","",laml@data$Amino_acids))

laml@data <- laml@data %>%
  add_column(POSITION = laml@data$Protein_position)

laml@data <- laml@data %>%
  add_column(Protein_Change = paste(laml@data$START,laml@data$POSITION,laml@data$END))

laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = ' ', "", laml@data$Protein_Change)
laml@data$Protein_Change=sub(pattern = "NA", "", laml@data$Protein_Change)

laml@data$Protein_Change<-paste0("p.",laml@data$Protein_Change)
laml@data$Protein_Change[laml@data$Protein_Change == "p."] <- ""
laml@data$Protein_Change[laml@data$Protein_Change %in% c("-", "?")] <- ""

laml@data$Protein_Change[grepl("[-?*]", laml@data$Protein_Change)] <- ""


laml@data <- laml@data %>%
  add_column(AA_mutation = paste0("p.",laml@data$Protein_position,laml@data$START,">",laml@data$END, " (", laml@data$Variant_Type, " - ", laml@data$Consequence, ", position ",  laml@data$Protein_position, ", ", laml@data$START, " ➞ ", laml@data$END,")"))

laml@data <- laml@data %>%
  add_column(CDS_mutation = paste0("c.",laml@data$CDS_position,laml@data$Reference_Allele,">",laml@data$Tumor_Seq_Allele2, " (", laml@data$Variant_Type, " - ", laml@data$Consequence, ", position ",  laml@data$CDS_position, ", ", laml@data$Reference_Allele, " ➞ ", laml@data$Tumor_Seq_Allele2, ")"))


df <- data.frame(Tumor_Sample_Barcode = laml@data$Tumor_Sample_Barcode,
                 Hugo_Symbol = laml@data$Hugo_Symbol,
                 Chromosome = laml@data$Chromosome,
                 AA_mutation = laml@data$AA_mutation,
                 CDS_mutation = laml@data$CDS_mutation)

df


write.table(df, file = "driver_mutations_non_silent_radon-rats.txt", sep = "\t", row.names = FALSE)


#####################################################################################################################



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenVisR")



install.packages("devtools")
install.packages(c("rmarkdown", "knitr", "testthat"))
devtools::install_github("griffithlab/GenVisR")


library(GenVisR)


lolliplot_mutationObs








