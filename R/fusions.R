setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/fusions")

library(ggplot2)
library(dplyr)
library(circlize)
library(gridExtra)
library(plotly)
library(reshape2)


########################## ########################## ########################## FUSIONS BARPLOT MEAN
df <- data.frame(group = c("D1", "D3", "D6", "D12", "Pr"), 
                 mean_fusions = c(12.44444444, 12.5, 13.25, 29, 1))
group_order <- c("D1", "D3", "D6", "D12", "Pr")
df$group <- factor(df$group, levels = group_order)

mean(df$mean_fusions)
# Create a barplot
ggplot(df, aes(x = group, y = mean_fusions, fill=group)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = mean(mean_fusions)), linetype = "dashed", color = "red") +
  xlab("Group") +
  ylab("Mean Fusions per Group") +
  ggtitle("Mean Fusions per Group 'rn6'") +
  theme_minimal()+
  labs(caption = "mean = 13.64")+
  theme(legend.position = "none")+
  scale_fill_viridis_d()


########################## ########################## ##########################
########################## ########################## ########################## FUSIONS BARPLOT MEAN
df <- data.frame(group = c("D1", "D3", "D6", "D12", "Pr"), 
                 mean_fusions = c(14, 14.83333333, 16.75, 83.875, 5.5))
group_order <- c("D1", "D3", "D6", "D12", "Pr")
df$group <- factor(df$group, levels = group_order)

mean(df$mean_fusions)
# Create a barplot
ggplot(df, aes(x = group, y = mean_fusions, fill=group)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = mean(mean_fusions)), linetype = "dashed", color = "red") +
  xlab("Group") +
  ylab("Mean Fusions per Group") +
  ggtitle("Mean Fusions per Group 'hg38'") +
  theme_minimal()+
  labs(caption = "mean = 26.99")+
  theme(legend.position = "none")+
  scale_fill_viridis_d()


########################## ########################## ##########################




data<-read.csv(file="fusions_radon-rats.CSV",sep=",",header=T)

data<-melt(data,id.vars = "Samples")
ggplot(data=data,aes(x=Samples,y=value, fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#201923","#ffffff","#fcff5d","#7dfc00","#0ec434","#228c68","#8ad8e8","#235b54","#29bdab","#3998f5","#37294f","#277da7","#3750db","#f22020","#991919","#ffcba5","#e68f66","#c56133","#96341c","#632819","#ffc413","#f47a22","#2f2aa0","#b732cc","#772b9d","#f07cab","#d30b94","#edeff3","#c3a5b4","#946aa2","#5d4c86"))+
  ggtitle("Predicted gene fusion effects")


########################### SNAKEY
links<-read.csv(file="fusion_partners_sample_count.CSV",sep=",",header=T)
  
labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

links<-links[-c(69:253),]

plot_ly(type = "sankey",
        node = list(label = nodes$label),
        link = list(source = match(links$from, nodes$label) - 1,
                    target = match(links$to, nodes$label) - 1,
                    value = links$value))

########################### CHORD DIAGRAM
circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")

#circos.genomicLabels(links, labels.column = 1, side = "outside")




################################################################################
fusions<-read.table(file="all_fusions3.TXT",sep="\t",header=F)
colnames(fusions) <- c("chr1", "po1", "gene1", "chr2", "po2", "gene2");
cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$po1,
                   end = fusions$po1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$po2,
                   end = fusions$po2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)
bed1 <- bed1 %>% 
  group_by(value1) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(value1 = ifelse(n < 25, NA, value1)) %>% 
  select(-n)

# Group by "value1" and add a "rank" column based on the order of rows within the group
bed1_ranked <- bed1 %>%
  group_by(value1) %>%
  mutate(rank = row_number())

# Convert "value1" to NA for rows with rank > 1 within each group
bed1<-bed1_ranked %>%
  mutate(value1 = ifelse(rank > 1, NA, value1)) %>%
  select(-rank)

head(as.data.frame(bed1), n=100)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "23", cols)
cols <- gsub("Y", "24", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.initializeWithIdeogram()
circos.genomicLabels(bed1, labels.column = 4, side = "inside")

circos.genomicLink(bed1, bed2, col = cols, border = NA)
################################################################################

#ALK FUSIONS
# Install and load the Gviz package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Gviz")

positions <- c(29831320, 29831322, 29926080, 29909462, 29831260) # specify positions

library(GenomicRanges)
library(Gviz)

# Create a GRanges object with the positions
gr <- GRanges(seqnames = rep("2", length(positions)),
              ranges = IRanges(start = positions, end = positions))

# Create the data track
dat <- DataTrack(gr, name = "Markers")

# Create the genome axis track
gat <- GenomeAxisTrack()

# Plot the tracks
plotTracks(tracks = list(gat, dat),
           chromosome = "2",
           from = 29000000,
           to = 31000000)







########################### PARTNERS
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/fusions")

library(plotly)
links<-read.table(file="FUSION_PARTNERS_RATS_FINAL.TXT",sep="\t",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

links<-links[-c(30:305),]


########################### CHORD DIAGRAM
circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")

#circos.genomicLabels(links, labels.column = 1, side = "outside")

################################################################################
links<-read.table(file="FUSION_PARTNERS_RATS_FINAL_gene_id.TXT",sep="\t",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

links<-links[-c(15:317),]


########################### CHORD DIAGRAM
circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")

################################################################################
fusions<-read.table(file="all_fusions_rat_genome.TXT",sep="\t",header=F)
colnames(fusions) <- c("chr1", "po1", "gene1", "chr2", "po2", "gene2");
cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$po1,
                   end = fusions$po1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$po2,
                   end = fusions$po2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)
bed1 <- bed1 %>% 
  group_by(value1) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(value1 = ifelse(n < 25, NA, value1)) %>% 
  select(-n)

# Group by "value1" and add a "rank" column based on the order of rows within the group
bed1_ranked <- bed1 %>%
  group_by(value1) %>%
  mutate(rank = row_number())

# Convert "value1" to NA for rows with rank > 1 within each group
bed1<-bed1_ranked %>%
  mutate(value1 = ifelse(rank > 1, NA, value1)) %>%
  select(-rank)

head(as.data.frame(bed1), n=100)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "21", cols)
cols <- gsub("Y", "22", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.clear()
circos.initializeWithIdeogram(species = "rn6")
#circos.genomicLabels(bed1, labels.column = 4, side = "inside")
circos.genomicLink(bed1, bed2, col = cols, border = NA)
################################################################################


setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/fusions/groups")

#BY GROUPS


#D1
fusions<-read.csv(file="D1.CSV",sep=",",header=T)
cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$pos1,
                   end = fusions$pos1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$pos2,
                   end = fusions$pos2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "21", cols)
cols <- gsub("Y", "22", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.clear()
circos.initializeWithIdeogram(species = "rn6")
#circos.genomicLabels(bed1, labels.column = 4, side = "inside")
circos.genomicLink(bed1, bed2, col = cols, border = NA)
title("D1")
################################################################################

#D3

fusions<-read.csv(file="D3.CSV",sep=",",header=T)
cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$pos1,
                   end = fusions$pos1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$pos2,
                   end = fusions$pos2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "21", cols)
cols <- gsub("Y", "22", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.clear()
circos.initializeWithIdeogram(species = "rn6")
#circos.genomicLabels(bed1, labels.column = 4, side = "inside")
circos.genomicLink(bed1, bed2, col = cols, border = NA)
title("D3")


################################################################################

#D6

fusions<-read.csv(file="D6.CSV",sep=",",header=T)
cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$pos1,
                   end = fusions$pos1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$pos2,
                   end = fusions$pos2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "21", cols)
cols <- gsub("Y", "22", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.clear()
circos.initializeWithIdeogram(species = "rn6")
#circos.genomicLabels(bed1, labels.column = 4, side = "inside")
circos.genomicLink(bed1, bed2, col = cols, border = NA)

title("D6")

################################################################################

#D12

fusions<-read.csv(file="D12.CSV",sep=",",header=T)
cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$pos1,
                   end = fusions$pos1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$pos2,
                   end = fusions$pos2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "21", cols)
cols <- gsub("Y", "22", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.clear()
circos.initializeWithIdeogram(species = "rn6")
#circos.genomicLabels(bed1, labels.column = 4, side = "inside")
circos.genomicLink(bed1, bed2, col = cols, border = NA)

title("D12")


################################################################################

#Pr

fusions<-read.csv(file="Pr.CSV",sep=",",header=T)
cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$pos1,
                   end = fusions$pos1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$pos2,
                   end = fusions$pos2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "21", cols)
cols <- gsub("Y", "22", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.clear()
circos.initializeWithIdeogram(species = "rn6")
#circos.genomicLabels(bed1, labels.column = 4, side = "inside")
circos.genomicLink(bed1, bed2, col = cols, border = NA)
title("Pr")




################################################################################
links<-read.table(file="D1.txt.TXT",sep="\t",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
title("D1")
################################################################################
links<-read.table(file="D3.txt.TXT",sep="\t",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
title("D3")
################################################################################
links<-read.table(file="D6.txt.TXT",sep="\t",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
title("D6")
################################################################################
links<-read.table(file="D12.txt.TXT",sep="\t",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
title("D12")
################################################################################
links<-read.table(file="Pr.txt.TXT",sep="\t",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

par(cex = 1, mar = c(0, 0, 0, 0))
circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(links$xmax)[i],by=5), 
            minor.ticks=1, labels.away.percentage = 0.15)
title("Pr")





########################################################################################################################################################################################################
#RAT FUSIONS

library(maftools)
library(data.table)
library(rstatix)

# load the data frame
df <- data.frame(
  group = c("CTRL", "CTRL", "D1", "D1", "D1", "D1", "D1", "D1", "D1", "D1",
            "D1", "D3", "D3", "D3", "D3", "D3", "D3", "D6", "D6", "D6", "D6",
            "D12", "D12", "D12", "D12", "D12", "D12", "D12", "D12", "D12", "Pr", "Pr", "Pr", "Pr"),
  sample = c("44514-ARN", "44541-ARN", "RnD1-36", "RnD1-81", "RnD1-104", "RnD1-132",
             "RnD1-134", "RnD1-162", "RnD1-195", "RnD1-207", "RnD1-214", "RnD3-168",
             "RnD3-203", "RnD3-212", "RnD3-232", "RnD3-233", "RnD3-92", "RnD6-132",
             "RnD6-152", "RnD6-179", "RnD6-204", "RnD12-49", "RnD12-51", "RnD12-145",
             "RnD12-178", "RnD12-184", "RnD12-185-1", "RnD12-185-2", "RnD12-229", "RnD12-231",
             "RnPr109", "RnPr140", "RnPr146", "RnPr214"),
  n_fusions = c(6, 144, 19, 22, 16, 5, 19, 2, 12, 14, 3, 9, 18, 20, 9, 6, 13,
                6, 22, 20, 5, 44, 135, 1, 6, 3, 9, 23, 6, 33, 0, 3, 1, 0)
)

pairwise_t_test(data=df,n_fusions ~ group, p.adjust.method = "bonferroni")

# check the structure of the data frame
str(df)

# perform the one-way ANOVA test
fit <- aov(n_fusions ~ group, data = df)
summary(fit)

# check the assumptions of the ANOVA test
plot(fit)

# perform the Tukey HSD post-hoc test
TukeyHSD(fit)





















setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-MINERS/genefuse/")

################################################################################ FUSIONS POSITIONS RADON-MINERS
fusions<-read.table(file="fusions_positions.TXT",sep="\t",header=T)

#REMOVE FUSIONS WITH LESS THAN 10 READS
fusions<-fusions[-c(3,9,12,13,17,20,36,46,49,52,53,67,68,69,71,72),]

cols<-fusions$chr1
fusions$chr1<-paste0("chr",fusions$chr1) 
fusions$chr2<-paste0("chr",fusions$chr2)

head(fusions)
# create bed1 dataframe
bed1 <- data.frame(chr = fusions$chr1,
                   start = fusions$pos1,
                   end = fusions$pos1,
                   value1 = fusions$gene1)

# create bed2 dataframe
bed2 <- data.frame(chr = fusions$chr2,
                   start = fusions$pos2,
                   end = fusions$pos2,
                   value1 = fusions$gene2)

colnames(bed1) <- c("chr", "start", "end", "value1");
colnames(bed2) <- c("chr", "start", "end", "value1");

head(bed1, n=30)
bed1 <- bed1 %>% 
  group_by(value1) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(value1 = ifelse(n < 25, NA, value1)) %>% 
  select(-n)

# Group by "value1" and add a "rank" column based on the order of rows within the group
bed1_ranked <- bed1 %>%
  group_by(value1) %>%
  mutate(rank = row_number())

# Convert "value1" to NA for rows with rank > 1 within each group
bed1<-bed1_ranked %>%
  mutate(value1 = ifelse(rank > 1, NA, value1)) %>%
  select(-rank)

head(as.data.frame(bed1), n=100)

# Replace "X" with "21" and "Y" with "22"
cols <- gsub("X", "23", cols)
cols <- gsub("Y", "24", cols)

# Convert to numeric format
cols <- as.numeric(cols)

circos.clear()
circos.initializeWithIdeogram(species = "hg18")
#circos.genomicLabels(bed1, labels.column = 4, side = "inside")
circos.genomicLink(bed1, bed2, col = cols, border = NA)




################################################################################ PARTNERS RADON-MINERS
links<-read.csv(file="fusions_partners.CSV",sep=",",header=T)
links<-read.csv(file="fusions_partners10reads.CSV",sep=",",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")

########################################################
########################################################
links<-read.csv(file="only_fusions_partners.CSV",sep=",",header=T)
links<-read.csv(file="only_fusions_partners10reads.CSV",sep=",",header=T)

labels<- c(union(unique(links$from),unique(links$to)))
nodes <- data.frame(label = labels)

circos.clear()
chordDiagram(links, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")




############PREDICTED EFFECT
setwd("C:/Users/bedmar/Desktop/IDIBAPS/RADON-RATS/fusions/")


pe_df<-read.csv(file="predicted_effect.csv",sep=",",header=T)

pe<-round(prop.table(table(pe_df$Predicted_effect)) * 100, 2)
pe<-as.data.frame(pe)
colnames(pe) <- c("Predicted_effect", "count")

library(scales)

ggplot(data=pe, aes(x="", y=count, fill=reorder(Predicted_effect, -count)))+
  geom_bar(stat = "identity", color="white")+
  coord_polar(theta = "y")+
  geom_text(data=pe,
            aes(label=percent(count/100)),
            position = position_stack(vjust = 0.5), color="white", size=3.5)+
  theme_void()+
  scale_fill_viridis_d(name="Predicted fusion effect")













####################################################RADON RATS FUSIONS INTER INTRA CHROMOSOMAL
# Load required libraries
library(ggplot2)
library(reshape2)

# Create the dataframe with the data
df <- data.frame(
  group = c("D1", "D3", "D6", "D12", "Pr"),
  inter_relative = c(10.56, 10.33, 10.25, 25.44, 0.50),
  intra_relative = c(1.89, 2.17, 3.00, 3.56, 0.50)
)
df_long<-melt(df,id.vars = "group")

# Calculate mean inter and intra values overall
mean_inter <- mean(df$inter_relative)
mean_intra <- mean(df$intra_relative)

df_long$group <- factor(df_long$group, levels = unique(df_long$group))
# Create the barplot
ggplot(data = df_long, aes(y = value, x = group, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = mean_inter, linetype = "dashed", color = "#338b85") +
  geom_hline(yintercept = mean_intra, linetype = "dashed", color = "#9ce0db") +
  scale_fill_manual(values = c("#338b85", "#9ce0db"), name = NULL, labels = c("Inter chromosomal", "Intra chromosomal")) +
  theme_minimal() +
  labs(x = "Rn group", y = "Relative to n of samples", fill = NULL)+
  theme(legend.position = "top")
  
+l