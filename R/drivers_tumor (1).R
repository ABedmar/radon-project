setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/drivers/tumor/total")
getwd()
library(maftools)
library(RColorBrewer)
library(stringr)
library(tibble)
library(mclust)
library(ggplot2)
library(data.table)
library(tidyr)

load("maf.Rdata")

laml@variant.type.summary


data<- data.frame(laml@data[,c(2,17)])
write.table(data, file = "drivers")

data2<-melt(data, id.vars = unique(c("Tumor_Sample_Barcode", "Hugo_Symbol")))

data<-data.frame(data)
duplicated[(data$Hugo_Symbol), ]

data<-laml@gene.summary
laml@gene.summary<-laml@gene.summary[order(-laml@gene.summary$total),]
laml@gene.summary<-laml@gene.summary[1:50,]

ggplot(laml@data, aes(x=laml@data$Tumor_Sample_Barcode, y=which(table(laml@data$Hugo_Symbol)>1)))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



####LOAD DATA RATS
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/drivers_orthologs/")
drivers<-read.csv(file = "Radon_rats_DRIVERS-group.CSV", sep = ",", header = TRUE)
drivers<-read.csv(file = "Radon_rats_DRIVERS-dose.CSV", sep = ",", header = TRUE)
drivers<-read.csv(file = "Radon_rats_DRIVERS-duration.CSV", sep = ",", header = TRUE)

drivers <- melt(drivers, id.vars="gene")

ggplot(data=drivers)+
  geom_col(aes(x=reorder(gene, -value),value,fill=variable), position = "dodge")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("# mutations")



drivers<-read.csv(file = "Radon_rats_DRIVERS-group.CSV", sep = ",", header = TRUE)
drivers<-drivers[-c(2,4,5,6,7,8,9,11,13,18,19,20,21,22),]
drivers<-drivers[22,]
ggplot(data=drivers, aes(fill=variable, y=value, x=gene))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("")+
  ylab("# mutations")
  
mean_D1=mean(drivers$D1.n.12.)
mean_D3=mean(drivers$D3.n.7.)
mean_D12=mean(drivers$D12.n.9.)
mean_D6=mean(drivers$D6.n.4.)
mean_Pr=mean(drivers$Pr.n.4.)


drivers <- melt(drivers, id.vars="gene")

head(drivers)

ggplot(data=drivers)+
  geom_col(aes(x=reorder(gene, -value),value,fill=variable), position = "dodge")+
  geom_hline(yintercept = mean_D1, color = "#f8766d", linetype="dashed")+
  geom_hline(yintercept = mean_D3, color = "#00bf7d", linetype="dashed")+
  geom_hline(yintercept = mean_D12, color = "#a3a500", linetype="dashed")+
  geom_hline(yintercept = mean_D6, color = "#00b0f6", linetype="dashed")+
  geom_hline(yintercept = mean_Pr, color = "#e76bf3", linetype="dashed")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("# mutations")+
  ggtitle("Group: D1 vs. D12 vs. D3 vs. D6 vs. Pr")


drivers<-read.csv(file = "Radon_rats_DRIVERS-dose.CSV", sep = ",", header = TRUE)
drivers<-drivers[-c(22),]

mean_high=mean(drivers$high.n.28.)
mean_low=mean(drivers$low.n.8.)

drivers <- melt(drivers, id.vars="gene")

head(drivers)

ggplot(data=drivers)+
  geom_col(aes(x=reorder(gene, -value),value,fill=variable), position = "dodge")+
  geom_hline(yintercept = mean_high, color = "#f8766d", linetype="dashed")+
  geom_hline(yintercept = mean_low, color = "#00bfc4", linetype="dashed")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("# mutations")+
  ggtitle("Dose: high vs. low")






drivers<-read.csv(file = "Radon_rats_DRIVERS-duration.CSV", sep = ",", header = TRUE)
drivers<-drivers[-c(22),]

mean_long=mean(drivers$long.n.9.)
mean_short=mean(drivers$short.n.23.)
mean_D6=mean(drivers$D6.n.4.)

drivers <- melt(drivers, id.vars="gene")

head(drivers)

ggplot(data=drivers)+
  geom_col(aes(x=reorder(gene, -value),value,fill=variable), position = "dodge")+
  geom_hline(yintercept = mean_long, color = "#f8766d", linetype="dashed")+
  geom_hline(yintercept = mean_short, color = "#00ba38", linetype="dashed")+
  geom_hline(yintercept = mean_D6, color = "#619cff", linetype="dashed")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("# mutations")+
  ggtitle("Duration: long vs. short vs. D6")


################################################################################



#HIGH VS LOW
drivers<-read.csv(file = "high_vs_low.CSV", sep = ",", header = TRUE)

mean_drivers<-drivers[-21,-c(2,3,4)]
drivers<-drivers[-21,-c(4,5)]

mean_high=mean(drivers$high_dose)
mean_low=mean(drivers$low_dose)
  
drivers <- melt(drivers, id.vars="Hugo_Symbol")

head(drivers)


ggplot(data=drivers)+
  geom_col(aes(x=reorder(Hugo_Symbol, -value),value,fill=variable), position = "dodge")+
  geom_hline(aes(yintercept= mean_high, linetype = "dose group mean"), colour= '#440154FF') +
  scale_linetype_manual(name ="", values = c('dashed'))+
  geom_hline(yintercept = mean_low, color = "#FDE725FF", linetype="dashed")+
  geom_line(data=mean_drivers, aes(x=Hugo_Symbol, y=mean, group=1), color="#20b2aa")+
  geom_point(data=mean_drivers, aes(x=Hugo_Symbol, y=mean), color="#20b2aa")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  scale_fill_viridis_d()+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("Mutations")

ggplot(data=drivers)+
  geom_col(aes(x=Hugo_Symbol,value,fill=variable), position = "dodge")+
  geom_hline(aes(yintercept= mean_high, linetype = "dose group mean"), colour= '#440154FF') +
  scale_linetype_manual(name ="", values = c('dashed'))+
  geom_hline(yintercept = mean_low, color = "#FDE725FF", linetype="dashed")+
  geom_line(data=mean_drivers, aes(x=Hugo_Symbol, y=mean, group=1), color="#20b2aa")+
  geom_point(data=mean_drivers, aes(x=Hugo_Symbol, y=mean), color="#20b2aa")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  scale_fill_viridis_d()+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("Mutations")



#HIGH VS INTERMEDIATE VS LOW
drivers2<-read.csv(file = "high_vs_intermediate_vs_low.CSV", sep = ",", header = TRUE)


mean_drivers2<-drivers2[-21,-c(2,3,4,5)]
drivers2<-drivers2[-21,-c(5,6)]

mean_high=mean(drivers2$high_dose)
mean_intermediate=mean(drivers2$intermediate_dose)
mean_low=mean(drivers2$low_dose)

drivers2 <- melt(drivers2, id.vars="Hugo_Symbol")

head(drivers2)

ggplot(data=drivers2)+
  geom_col(aes(x=reorder(Hugo_Symbol, -value),value,fill=variable), position = "dodge")+
  geom_hline(aes(yintercept= mean_high, linetype = "dose group mean"), colour= '#440154FF') +
  scale_linetype_manual(name ="asdsad", values = c('dashed'))+
  geom_hline(yintercept = mean_intermediate, color = "#1F968BFF", linetype="dashed")+
  geom_hline(yintercept = mean_low, color = "#FDE725FF", linetype="dashed")+
  geom_line(data=mean_drivers2, aes(x=Hugo_Symbol, y=mean, group=1), color="#20b2aa")+
  geom_point(data=mean_drivers2, aes(x=Hugo_Symbol, y=mean), color="#20b2aa")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  scale_fill_viridis_d()+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("Mutations")







# SCATTER PLOT
drivers<-read.csv(file = "high_vs_low.CSV", sep = ",", header = TRUE)

mean_drivers<-drivers[-21,-c(2,3,4)]
drivers<-drivers[-21,-c(4,5)]

mean_high=mean(drivers$high_dose)
mean_low=mean(drivers$low_dose)

drivers <- melt(drivers, id.vars="Hugo_Symbol")

head(drivers)


ggplot(data=drivers)+
  geom_col(aes(x=reorder(Hugo_Symbol, -value),value,fill=variable), position = "dodge")+
  geom_hline(aes(yintercept= mean_high, linetype = "dose group mean"), colour= '#440154FF') +
  scale_linetype_manual(name ="", values = c('dashed'))+
  geom_hline(yintercept = mean_low, color = "#FDE725FF", linetype="dashed")+
  geom_line(data=mean_drivers, aes(x=Hugo_Symbol, y=mean, group=1), color="#20b2aa")+
  geom_point(data=mean_drivers, aes(x=Hugo_Symbol, y=mean), color="#20b2aa")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  scale_fill_viridis_d()+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("Mutations")





############################################## NORMALIZATION
library("limma")
library("edgeR")

counts <- read.csv(file = "drivers2.csv", sep = ",", header = T)

rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[-nrow(counts),-ncol(counts)]

counts <- data.matrix(counts)

dge <- DGEList(counts=counts)
design <- model.matrix(~ 0+factor(colnames(counts)))
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)
topTable(fit, coef=ncol(design))
write.table(v$E,file = "counts_normalized.txt", sep = "\t") 


v$E

##Plot unnormalized vs normalized data
par(mfrow=c(1,2))
boxplot(counts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(counts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")








####### TTN
TTN<-read.csv(file = "TTN.CSV", sep = ",", header = TRUE)

ggplot(data=TTN)+
  geom_col(aes(x=reorder(Sample, -TTN_count), y=TTN_count, fill=Radon.dose))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "bottom")+
  ggtitle("Pathogenic TTN mutations (samples)")+
  xlab(element_blank())+
  ylab("# TTN mutations")

means<-aggregate(TTN$TTN_count, list(TTN$Radon.dose), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(TTN$Sample, by=list(TTN$Radon.dose), FUN=length)

ggplot(data=TTN)+
  geom_boxplot(aes(x=reorder(Radon.dose, -TTN_count), y=TTN_count))+
  geom_text(data = means, aes(x=Group.1, label =  paste("Mean:", x), y = 6))+
  geom_text(data = N, aes(x=Group.1, label =  paste("N =", x), y = 6.3))+
  geom_jitter(aes(x=Radon.dose, y=TTN_count, color=Smoking.status))+
  theme_minimal()+ 
  ggtitle("Pathogenic TTN mutations (radon dose)")+
  xlab(element_blank())+
  ylab("# TTN mutations")


ggplot(data=TTN)+
  geom_col(aes(x=reorder(Sample, -TTN_count), y=TTN_count, fill=Smoking.status))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "bottom")+
  ggtitle("Pathogenic TTN mutations (samples)")+
  xlab(element_blank())+
  ylab("# TTN mutations")



means<-aggregate(TTN$TTN_count, list(TTN$Smoking.status), FUN=mean)
means[2]<-round(means$x, digits = 2)
N <- aggregate(TTN$Sample, by=list(TTN$Smoking.status), FUN=length)

ggplot(data=TTN)+
  geom_boxplot(aes(x=reorder(Smoking.status, -TTN_count), y=TTN_count, color=Smoking.status))+
  geom_text(data = means, aes(x=Group.1, label =  paste("Mean:", x), y = 6))+
  geom_text(data = N, aes(x=Group.1, label =  paste("N =", x), y = 6.3))+
  geom_jitter(aes(x=Smoking.status, y=TTN_count, colour=Radon.dose, shape=Molecular.smoking))+
  theme_minimal()+ 
  ggtitle("Pathogenic TTN mutations (Smoking status)")+
  xlab(element_blank())+
  ylab("# TTN mutations")



#MINERS
#####################################################################DRIVERS ANALYSIS COMO EN RATS
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/")
load("maf.Rdata")
load("maf_pathogenic.Rdata")

############################################# PREPARE DATA FOR ALL ANALYSIS

drivers<-c("Tp53","Kras","Egfr","Stk11","Pik3Ca","Erbb4","Ntrk3","Braf","Ros1","Alk","Pten","Met","Erbb2","Ret","Ntrk1","Nrg1","Ntrk2","Fgfr2","Fgfr4","Erbb3","Fgfr1","Nras","Fgfr3","Mapk1","Atm","Hprt")
drivers<-toupper(drivers)
dim(laml@data)


laml_drivers <- laml@data[laml@data$Hugo_Symbol %in% drivers]

braf<-laml_drivers[laml_drivers$Hugo_Symbol %in% "BRAF"]

df <- data.frame(Hugo_Symbol = laml_drivers$Hugo_Symbol,
                 Tumor_Sample_Barcode = laml_drivers$Tumor_Sample_Barcode)
# Create a new data frame with unique Hugo Symbols
unique_genes <- data.frame(Hugo_Symbol = unique(df$Hugo_Symbol))
# Rename columns in df
colnames(df) <- c("Hugo_Symbol", "Tumor_Sample_Barcode")
# Melt the data frame by Hugo Symbol
melted_df <- melt(df, id.vars = "Hugo_Symbol")

write.table(table(melted_df$Hugo_Symbol, melted_df$value), file = "Radon_miners_DRIVERS_pathogenic.txt", sep="\t")





setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/drivers")
################################################################
#high dose
data<-read.table(file = "high_pathogenic.tsv", header = T, sep = "\t")

head(data)
# Set the order of drivers and reverse it
driver_order <- rev(c("TP53","EGFR","ATM","ALK","ERBB4","RET","BRAF","FGFR1","MET","NTRK1","PIK3CA","PTEN","FGFR3","ROS1","ERBB3","KRAS","NRAS","NRG1","NTRK2"))

# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "brown1") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("High dose")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
heatmap

# Add percentages as annotations per row
heatmap_with_labels_high <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
heatmap_with_labels_high




################################################################
#inter dose
data<-read.table(file = "inter_pathogenic.tsv", header = T, sep = "\t")


# Set the order of drivers and reverse it
driver_order <- rev(c("TP53","EGFR","ATM","ALK","ERBB4","RET","BRAF","FGFR1","MET","NTRK1","PIK3CA","PTEN","FGFR3","ROS1","ERBB3","KRAS","NRAS","NRG1","NTRK2"))

# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "goldenrod1") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("Intermediate dose")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
heatmap

# Add percentages as annotations per row
heatmap_with_labels_inter <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
heatmap_with_labels_inter



################################################################
#low dose
data<-read.table(file = "low_pathogenic.tsv", header = T, sep = "\t")


# Set the order of drivers and reverse it
driver_order <- rev(c("TP53","EGFR","ATM","ALK","ERBB4","RET","BRAF","FGFR1","MET","NTRK1","PIK3CA","PTEN","FGFR3","ROS1","ERBB3","KRAS","NRAS","NRG1","NTRK2"))

# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "darkolivegreen3") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("Low dose")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
heatmap

# Add percentages as annotations per row
heatmap_with_labels_low <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
heatmap_with_labels_low

combined_plot <- heatmap_with_labels_high / heatmap_with_labels_inter / heatmap_with_labels_low + plot_layout(widths = c(15, 10, 5), ncol = 3)

combined_plot





#####################################################################################################################



setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/")
load("maf_pathogenic.Rdata")

laml@data <- laml@data[laml@data$Hugo_Symbol %in% drivers]

nrow(laml@data)
unique(laml@data$Tumor_Sample_Barcode)
unique(laml@data$Hugo_Symbol)

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


write.table(df, file = "pathogenic_drivers_protein_change.txt", sep = "\t", row.names = FALSE)


#####################################################################################################################










setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/cancer_predisposition_pathogenic/")

load("cancer_predisposition_N_maf.Rdata")

df <- data.frame(Hugo_Symbol = laml@data$Hugo_Symbol,
                 Tumor_Sample_Barcode = laml@data$Source_MAF)

df$Tumor_Sample_Barcode <- gsub("\\.maf", "", df$Tumor_Sample_Barcode)

# Create a new data frame with unique Hugo Symbols
unique_genes <- data.frame(Hugo_Symbol = unique(df$Hugo_Symbol))
# Rename columns in df
colnames(df) <- c("Hugo_Symbol", "Tumor_Sample_Barcode")
# Melt the data frame by Hugo Symbol
melted_df <- melt(df, id.vars = "Hugo_Symbol")

write.table(table(melted_df$Hugo_Symbol, melted_df$value), file = "Radon_miners_cancer_predisposition_pathogenic.txt", sep="\t")

################################################################
#cancer predisposition genes
data<-read.table(file = "Radon_miners_cancer_predisposition_pathogenic.tsv", header = T, sep = "\t")


# Set the order of drivers and reverse it
driver_order <- rev(c("TP53","APC","ATM","BRCA2","NF1","EGFR","BRAF","PTCH1","RB1","CDH1","JMJD1C","PTEN","COL7A1","KIT","SETBP1","ATR","AXIN2","BRIP1","CBL","DICER1","ERCC3","FANCA","FANCI","MET","MSH2","PDGFRA","SMAD4","VHL","WT1","RET","SMARCA4","ALK","AP002956.1","CTR9","CYLD","FAM71B","MARCH9","PMS1","RAF1","RHBDF2","STAT3","TGFBR1","UROD","WRN","ABCB11","B3GALT2","BRCA1","CDC73","CDKN1B","CFD","CSDE1","DOCK8","DROSHA","ELANE","ERCC2","ERCC4","ETFRF1","FAM20A","FANCE","FANCM","FH","FLCN","GPC3","HMBS","ITK","LRRC14","MSH6","MUTYH","NRAS","PALB2","PAX5","POLE","POLH","RAD51C","SDHA","SMARCB1","SOS1","STK11","TERT","TMEM127","TOE1","TSC1","TSHR","TSPAN31","WAS","XPA"))

# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("Cancer predisposition pathogenic genes in normal tissue")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
heatmap

# Add percentages as annotations per row
heatmap_with_labels_low <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
heatmap_with_labels_low







################################################################
load("cancer_predisposition_N_maf.Rdata")
load("cancer_predisposition_maf.Rdata")

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


df <- data.frame(Normal_Sample_Barcode = laml@data$Source_MAF,
                 Hugo_Symbol = laml@data$Hugo_Symbol,
                 Chromosome = laml@data$Chromosome,
                 Start_Position = laml@data$Start_Position,
                 End_Position = laml@data$End_Position,
                 Strand = laml@data$Strand,
                 Consequence = laml@data$Consequence,
                 Protein_Change = laml@data$Protein_Change,
                 Variant_Type = laml@data$Variant_Type,
                 CDS_Change = laml@data$CDS_Change)

df$Normal_Sample_Barcode <- gsub("\\.maf", "", df$Normal_Sample_Barcode)

df


write.table(df, file = "AA_changes_cancer_predisposition_VERY_pathogenic_genes_in_normal_tissue.txt", sep = "\t", row.names = FALSE)






  ###
















drivers<-read.csv(file = "Radon_rats_DRIVERS.csv", sep = ",", header = T)
head(drivers)

df<-drivers[,c(1,3)]

df_summary <- df %>% 
  group_by(Hugo_Symbol, group) %>% 
  summarize(count = n()) %>% 
  ungroup()

head(df_summary)
# create a data frame with all possible combinations of Hugo_Symbol and group
all_combinations <- expand.grid(unique(df$Hugo_Symbol), unique(df$group))
colnames(all_combinations) <- c("Hugo_Symbol", "group")

# join the all_combinations data frame with df_summary to get the count for each combination
df_summary_new <- left_join(all_combinations, df_summary, by = c("Hugo_Symbol", "group"))
# replace NA values with 0
df_summary_new$count[is.na(df_summary_new$count)] <- 0
head(df_summary_new)


# Calculate mean count for each group
group_means <- df_summary_new %>%
  group_by(group) %>%
  summarize(mean_count = mean(count)) %>%
  ungroup()

df_summary_new$group <- factor(df_summary_new$group, levels = unique(df$group))

A<-ggplot(df_summary_new, aes(x = reorder(Hugo_Symbol, -count), y = count, fill = factor(group, levels = c("D1", "D3", "D6", "D12", "Pr")))) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  theme_minimal() +
  ylab("# Mutations")+
  geom_hline(yintercept = group_means$mean_count[1], color = "#440154", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[2], color = "#5dc863", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[3], color = "#3b528b", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[4], color = "#21908c", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[5], color = "#fde725", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ggtitle("Gene drivers by radon group")























