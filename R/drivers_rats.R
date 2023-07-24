setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/cumulative_dose/")

getwd()
library(maftools)
library(RColorBrewer)
library(stringr)
library(tibble)
library(mclust)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(stringr)


######################COMMULATIVE

drivers<-c("AKT1","ALK","APC","ARHGAP35","ARID14","ARID1A","ARID2","ATM","ATRX","BRAF","CDKN2A","CHD2","CR2","CTNNB1","CUL3","DHX30","DICER1","EGFR","EIF1AX","EP300","ERBB2","ERBB3","ERBB4","FAT1","FBXW7","FGFR1","FGFR2","FGFR3","FGFR4","GNAS","HLA-A","HPRT","HRAS","HUWE1","IGF2","IGF2BP3","KDM6A","KEAP1","KLF5","KMT2D","KRAS","LRP1B","LTK","MAP2K1","MAPK1","MEN1","MET","MGA","NF1","NFE2L2","NKX2-1","NOTCH1","NRAS","NRG1","NTRK1","NTRK2","NTRK3","PDGFRA","PIK3CA","PKHD1L1","PPARG","PTEN","RASA1","RB1","RBM10","RET","RIT1","RM10","ROS1","SETD2","SMAD4","SMARCA4","STK11","TENM3","TERT","TG","TP53","TSC1","TSC2","TSHR","U2AF1","UBA1")
load("maf_comulative.Rdata")

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea2.1//")
load("maf.Rdata")
head(laml@data$Source_MAF)
dim(laml@data)

laml_drivers <- laml@data[laml@data$Hugo_Symbol %in% drivers]
colnames(laml@data)
dim(laml@data)
dim(a)
table(a$Hugo_Symbol)
merge_mafs()

dim(laml_drivers$Hugo_Symbol)
table(laml_drivers$Hugo_Symbol)
df<-table(laml@data$Hugo_Symbol)

laml_drivers <- df[drivers]

drivers <- data.frame(Hugo_Symbol = laml_drivers$Hugo_Symbol,
                      variable = laml_drivers$Source_MAF)

drivers<-drivers %>%
  group_by(Hugo_Symbol, variable) %>%
  summarise(value = n()) %>%
  ungroup()

drivers$variable<-gsub(".maf", "", drivers$variable)

head(drivers)
drivers<-complete(drivers, nesting(Hugo_Symbol), nesting(variable = factor(variable, levels = c("high_comultaive", "low_comulative"))), fill = list(value = 0))

means<-drivers %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value))

gene_means <- drivers %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(mean_value = mean(value))


#NORMALIZE BY NUMBER OF SAMPLES
drivers_modified <- drivers %>%
  mutate(value = ifelse(variable == "high_comultaive", value / 32, value / 4))
means<-drivers_modified %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value))
gene_means <- drivers_modified %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(mean_value = mean(value))


ggplot(data=drivers_modified)+
  geom_col(aes(x=reorder(Hugo_Symbol, -value),value,fill=variable), position = "dodge", width = 0.6)+
  geom_hline(yintercept = means$mean_value[1], color = "#440154FF", linetype="dashed")+
  geom_hline(yintercept = means$mean_value[2], color = "#FDE725FF", linetype="dashed")+
  geom_line(data=gene_means, aes(x=Hugo_Symbol, y=mean_value, group=1), color="#20b2aa")+
  geom_point(data=gene_means, aes(x=Hugo_Symbol, y=mean_value), color="#20b2aa")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  scale_fill_viridis_d(name="Mutation category", labels=c("High cumulative", "Low cumulative"))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("# Mutations")
################################################################################
################################################################################

######################PAEC


setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/PAEC_dose/")


load("maf_PAEC.Rdata")
laml_drivers <- laml@data[Hugo_Symbol %in% drivers]


drivers <- data.frame(Hugo_Symbol = laml_drivers$Hugo_Symbol,
                      variable = laml_drivers$Source_MAF)

drivers<-drivers %>%
  group_by(Hugo_Symbol, variable) %>%
  summarise(value = n()) %>%
  ungroup()

drivers$variable<-gsub(".maf", "", drivers$variable)

head(drivers)
drivers<-complete(drivers, nesting(Hugo_Symbol), nesting(variable = factor(variable, levels = c("high_PAEC", "low_PAEC"))), fill = list(value = 0))

means<-drivers %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value))

gene_means <- drivers %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(mean_value = mean(value))

ggplot(data=drivers)+
  geom_col(aes(x=reorder(Hugo_Symbol, -value),value,fill=variable), position = "dodge")+
  geom_hline(yintercept = means$mean_value[1], color = "#440154FF", linetype="dashed")+
  geom_hline(yintercept = means$mean_value[2], color = "#FDE725FF", linetype="dashed")+
  geom_line(data=gene_means, aes(x=Hugo_Symbol, y=mean_value, group=1), color="#20b2aa")+
  geom_point(data=gene_means, aes(x=Hugo_Symbol, y=mean_value), color="#20b2aa")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  scale_fill_viridis_d()+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ylab("# Mutations")










################################################################################ GENERATE DRIVERS TABLE
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/idea_1")

drivers<-c("AKT1","ALK","APC","ARHGAP35","ARID14","ARID1A","ARID2","ATM","ATRX","BRAF","CDKN2A","CHD2","CR2","CTNNB1","CUL3","DHX30","DICER1","EGFR","EIF1AX","EP300","ERBB2","ERBB3","ERBB4","FAT1","FBXW7","FGFR1","FGFR2","FGFR3","FGFR4","GNAS","HLA-A","HPRT","HRAS","HUWE1","IGF2","IGF2BP3","KDM6A","KEAP1","KLF5","KMT2D","KRAS","LRP1B","LTK","MAP2K1","MAPK1","MEN1","MET","MGA","NF1","NFE2L2","NKX2-1","NOTCH1","NRAS","NRG1","NTRK1","NTRK2","NTRK3","PDGFRA","PIK3CA","PKHD1L1","PPARG","PTEN","RASA1","RB1","RBM10","RET","RIT1","RM10","ROS1","SETD2","SMAD4","SMARCA4","STK11","TENM3","TERT","TG","TP53","TSC1","TSC2","TSHR","U2AF1","UBA1")
drivers <- str_to_title(tolower(drivers))

load("maf_v2.Rdata")
head(laml@data)
laml_drivers <- laml@data[Hugo_Symbol %in% drivers]

dim(laml_drivers)

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

head(laml_drivers)
laml_drivers$Protein_Change

simplified_drivers<-data.frame(Source_MAF = laml_drivers$Source_MAF,
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
write.table(simplified_drivers, file = "C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/cumulative_dose/simplified_EXTENDED_drivers.txt", sep = "\t")

df <- data.frame(Hugo_Symbol = laml_drivers$Hugo_Symbol,
                 Tumor_Sample_Barcode = laml_drivers$Tumor_Sample_Barcode)
# Create a new data frame with unique Hugo Symbols
unique_genes <- data.frame(Hugo_Symbol = unique(df$Hugo_Symbol))
# Rename columns in df
colnames(df) <- c("Hugo_Symbol", "Tumor_Sample_Barcode")
# Melt the data frame by Hugo Symbol
melted_df <- melt(df, id.vars = "Hugo_Symbol")

write.table(table(melted_df$Hugo_Symbol, melted_df$value), file = "Radon_Rats_EXTENDED_drivers.txt", sep="\t")



################################################################
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/ideas/drivers/")
#D1
data<-read.table(file = "D1_(2).tsv", header = T, sep = "\t")
data<-data[,-1]
head(data)
# Set the order of drivers and reverse it
driver_order <- rev(c("Nf1","Nras","Met","Ctnnb1","Erbb3","Erbb4","Fgfr2","Keap1","Akt1","Fgfr3","Fgfr4","Mapk1","Pten","Alk","Braf","Egfr","Erbb2","Fgfr1","Ntrk3","Ret","Ros1","Stk11","Hras","Kras","Map2k1","Nrg1","Ntrk1","Ntrk2","Pik3ca","Tp53",".",",","Huwe1","Mga","Pkhd1l1","Fat1","Kdm6a","Kmt2d","Atrx","Dicer1","Smarca4","Apc","Arid1a","Chd2","Gnas","Notch1","Smad4","Uba1","Cul3","Ep300","Ltk","Rasa1","Rbm10","Setd2","Tert","Arid2","Atm","Cr2","Igf2bp3","Nfe2l2","Tenm3","Tsc2","U2af1","Arhgap35","Fbxw7","Igf2","Men1","Nkx2-1","Tg","Tsc1","Cdkn2a","Hprt","Lrp1b","Pdgfra","Rb1","Rit1"))
# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))
sample_order<- c("D1.132","D1.214","D1.195","D1.162","D1.104","D1.207","D1.36","D1.81","D1.134","D1.98","D1.171","D1.189")
data_long$variable <- factor(data_long$variable, levels = sample_order)

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "#440154") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("D1")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
#heatmap

# Add percentages as annotations per row
D1_heatmap_with_labels <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
d1<-D1_heatmap_with_labels+theme(axis.text.x = element_text(margin = margin(t = 13.5)),
                                 plot.title = element_text(hjust = 0.5),
                                 text=element_text(family="Calibri")) +
  geom_text(aes(label = ifelse(value>0, value, "")),
            colour = ifelse(data_long$value > 4, "white", "black"), 
            size=3,
            vjust=0.4)

d1
################################################################
################################################################
#D3
data<-read.table(file = "D3_(2).tsv", header = T, sep = "\t")
data<-data[,-1]

head(data)

# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))
sample_order<- c("D3.232","D3.92","D3.168","D3.212","D3.203","D3.233","D3.71")
data_long$variable <- factor(data_long$variable, levels = sample_order)

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "#3b528b") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("D3")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
#heatmap

# Add percentages as annotations per row
D3_heatmap_with_labels <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
d3<-D3_heatmap_with_labels+theme(axis.text.x = element_text(margin = margin(t = 13.5)),
                                 plot.title = element_text(hjust = 0.5),
                                 text=element_text(family="Calibri")) +
  geom_text(aes(label = ifelse(value>0, value, "")),
            colour = ifelse(data_long$value > 1, "white", "black"), 
            size=3,
            vjust=0.4)
d3
################################################################
################################################################
#D6
data<-read.table(file = "D6_(2).tsv", header = T, sep = "\t")
data<-data[,-1]

head(data)
# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))
sample_order<- c("D6.132","D6.204","D6.152","D6.179")
data_long$variable <- factor(data_long$variable, levels = sample_order)

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "#21908c") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("D6")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
#heatmap

# Add percentages as annotations per row
D6_heatmap_with_labels <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
d6<-D6_heatmap_with_labels+theme(axis.text.x = element_text(margin = margin(t = 13.5)),
                                 plot.title = element_text(hjust = 0.5),
                                 text=element_text(family="Calibri")) +
  geom_text(aes(label = ifelse(value>0, value, "")),
            colour = ifelse(data_long$value > 2, "white", "black"), 
            size=3,
            vjust=0.4)
d6
################################################################
################################################################
#D12
data<-read.table(file = "D12_(2).tsv", header = T, sep = "\t")
data<-data[,-1]

head(data)

# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))
sample_order<- c("D12.145","D12.184","D12.178","D12.51","D12.231","D12.229","D12.185.1","D12.185.2","D12.49")
data_long$variable <- factor(data_long$variable, levels = sample_order)

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "#5dc863") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("D12")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
#heatmap

# Add percentages as annotations per row
D12_heatmap_with_labels <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
d12<-D12_heatmap_with_labels+theme(plot.title = element_text(hjust = 0.5),
                                   text=element_text(family="Calibri")) +
  geom_text(aes(label = ifelse(value>0, value, "")),
            colour = ifelse(data_long$value > 4, "white", "black"), 
            size=3,
            vjust=0.4)
d12
################################################################
################################################################
#Pr
data<-read.table(file = "Pr_(2).tsv", header = T, sep = "\t")
data<-data[,-1]

head(data)

# Convert 'driver' column to factor with reversed order
data$driver <- factor(data$driver, levels = driver_order)

# Reshape the data to long format
data_long <- melt(data, id.vars = c("driver", "percentage_of_samples_with_driver"))
sample_order<- c("Pr140","Pr109","Pr214","Pr146")
data_long$variable <- factor(data_long$variable, levels = sample_order)

# Calculate the unique percentages for each driver
unique_percentages <- aggregate(percentage_of_samples_with_driver ~ driver, data_long, unique)

# Create the heatmap
heatmap <- ggplot(data_long, aes(x = variable, y = driver)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "#fde725") +
  theme_minimal() +
  xlab(label = element_blank())+
  ylab(label = element_blank())+
  ggtitle("Pr")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")

# Display the heatmap
#heatmap

# Add percentages as annotations per row
Pr_heatmap_with_labels <- heatmap +
  geom_text(data = unique_percentages, aes(x = Inf, label = paste0(percentage_of_samples_with_driver, "%")),
            color = "black", size = 2.5, hjust = 1)

# Display the heatmap with labels
pr<-Pr_heatmap_with_labels+theme(axis.text.x = element_text(margin = margin(t = 19.5)),
                                 plot.title = element_text(hjust = 0.5),
                                 text=element_text(family="Calibri")) +
  geom_text(aes(label = ifelse(value>0, value, "")),
            colour = ifelse(data_long$value > 4, "white", "black"), 
            size=3,
            vjust=0.4)

pr
################################################################
# Combine the plots

combined_plot <- grid.arrange(d1, d3, d6, d12, pr,
                              ncol = 5, widths = c(12, 8, 6, 10, 6))

#Display combined plot
print(combined_plot)

help(grid.arrange)





aaaaaaaaaa












































































































































































































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
    
################################################################################
#CUMULATIVE
df<-drivers[,c(1,4)]

df_summary <- df %>% 
  group_by(Hugo_Symbol, CUMULATIVE) %>% 
  summarize(count = n()) %>% 
  ungroup()

head(df_summary)
all_combinations <- expand.grid(unique(df$Hugo_Symbol), unique(df$CUMULATIVE))
colnames(all_combinations) <- c("Hugo_Symbol", "CUMULATIVE")

df_summary_new <- left_join(all_combinations, df_summary, by = c("Hugo_Symbol", "CUMULATIVE"))
df_summary_new$count[is.na(df_summary_new$count)] <- 0
head(df_summary_new)


group_means <- df_summary_new %>%
  group_by(CUMULATIVE) %>%
  summarize(mean_count = mean(count)) %>%
  ungroup()


B<-ggplot(df_summary_new, aes(x = reorder(Hugo_Symbol, -count), y = count, fill = CUMULATIVE)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  theme_minimal() +
  ylab("# Mutations")+
  geom_hline(yintercept = group_means$mean_count[1], color = "#440154", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[2], color = "#fde725", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ggtitle("Gene drivers by cumulative dose")
################################################################################
#PAEC
df<-drivers[,c(1,5)]

df_summary <- df %>% 
  group_by(Hugo_Symbol, PAEC) %>% 
  summarize(count = n()) %>% 
  ungroup()

head(df_summary)
all_combinations <- expand.grid(unique(df$Hugo_Symbol), unique(df$PAEC))
colnames(all_combinations) <- c("Hugo_Symbol", "PAEC")

df_summary_new <- left_join(all_combinations, df_summary, by = c("Hugo_Symbol", "PAEC"))
df_summary_new$count[is.na(df_summary_new$count)] <- 0
head(df_summary_new)


group_means <- df_summary_new %>%
  group_by(PAEC) %>%
  summarize(mean_count = mean(count)) %>%
  ungroup()


C<-ggplot(df_summary_new, aes(x = reorder(Hugo_Symbol, -count), y = count, fill = PAEC)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  theme_minimal() +
  ylab("# Mutations")+
  geom_hline(yintercept = group_means$mean_count[1], color = "#440154", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[2], color = "#fde725", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ggtitle("Gene drivers by PEAC dose")
################################################################################
#HISTOLOGY
df<-drivers[,c(1,6)]

df_summary <- df %>% 
  group_by(Hugo_Symbol, histology) %>% 
  summarize(count = n()) %>% 
  ungroup()

head(df_summary)
all_combinations <- expand.grid(unique(df$Hugo_Symbol), unique(df$histology))
colnames(all_combinations) <- c("Hugo_Symbol", "histology")

df_summary_new <- left_join(all_combinations, df_summary, by = c("Hugo_Symbol", "histology"))
df_summary_new$count[is.na(df_summary_new$count)] <- 0
head(df_summary_new)


group_means <- df_summary_new %>%
  group_by(histology) %>%
  summarize(mean_count = mean(count)) %>%
  ungroup()


D<-ggplot(df_summary_new, aes(x = reorder(Hugo_Symbol, -count), y = count, fill = histology)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  theme_minimal() +
  ylab("# Mutations")+
  geom_hline(yintercept = group_means$mean_count[1], color = "#440154", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[2], color = "#21908c", linetype="dashed")+
  geom_hline(yintercept = group_means$mean_count[3], color = "#fde725", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.9, hjust=1))+
  theme(legend.title = element_blank())+
  xlab(element_blank())+
  ggtitle("Gene drivers by histology")


grid.arrange(A,D,B,C,ncol=2)




