setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/snp_density/")

library("handyFunctions")
library(ggplot2)
library(tools)

help(ShowSNPDensityPlot)

dir_path <- "C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/snp_density/"
pattern <- "*snpden"

# Get file names
file_names <- list.files(dir_path, pattern = pattern)

#RADON-RATS
for (file_name in file_names) {
  density_data <- read.table(paste0(dir_path, file_name), header = TRUE)
  file_title <- file_path_sans_ext(basename(file_name))
  
  plot <- ShowSNPDensityPlot(density_data, binSize = 1e6) + ggtitle(file_title)
  
  png(file = paste0(dir_path, gsub(pattern = pattern, replacement = ".png", x = file_name)),
      width = 850, height = 1080, res = 150)
  print(plot)
  dev.off()

}

file_names
data <- read.table(paste0(dir_path, "RnD12-135.snpden"), header = TRUE)
data <- data[data$CHROM == 15, ]
max(data$BIN_START)
unique(data$CHROM)
data <- data[data$CHROM == "X", ]


ggplot()+
  geom_line(data = data, aes(x = BIN_START, y = SNP_COUNT, colour=..y..), alpha = 0.4) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+
  scale_colour_gradientn(colors=c("gray", "green", "yellow","red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red", "red2","red2", "red2", "red2", "red2", "red4", "red4", "red4","red4", "red4", "red4", "red4", "red4", "red4", "red4", "red4"),
                         limits=c(0, 188))+ 
  scale_x_continuous(breaks = seq(0, 112000000000, by = 5000000))




chr10_snp <- ggplot()

# Iterate through the files and add each plot to the combined plot
for (file_name in file_names) {
  data <- read.table(paste0(dir_path, file_name), header = TRUE)
  file_title <- file_path_sans_ext(basename(file_name))
  
  data <- data[data$CHROM == 10, ]
  # Add the plot for the current file to the combined plot
  chr10_snp <- chr10_snp +
    geom_line(data = data, aes(x = BIN_START, y = SNP_COUNT, colour=..y..), alpha = 0.4) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+
    scale_colour_gradientn(colors=c("gray", "green", "yellow","red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red", "red2","red2", "red2", "red2", "red2", "red4", "red4", "red4","red4", "red4", "red4", "red4", "red4", "red4", "red4", "red4"),
                           limits=c(0, 188))+ 
    scale_x_continuous(breaks = seq(0, 112000000000, by = 4000000))
}


# Display the combined plot
chr10_snp <- chr10_snp +
  ggtitle("SNP density on chromosome 10")

chr10_snp




chr12_snp <- ggplot()

# Iterate through the files and add each plot to the combined plot
for (file_name in file_names) {
  data <- read.table(paste0(dir_path, file_name), header = TRUE)
  file_title <- file_path_sans_ext(basename(file_name))
  
  data <- data[data$CHROM == 12, ]
  # Add the plot for the current file to the combined plot
  chr12_snp <- chr12_snp +
    geom_line(data = data, aes(x = BIN_START, y = SNP_COUNT, colour=..y..), alpha = 0.4) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+
    scale_colour_gradientn(colors=c("gray", "green", "yellow","red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red", "red2","red2", "red2", "red2", "red2", "red4", "red4", "red4","red4", "red4", "red4", "red4", "red4", "red4", "red4", "red4"),
                           limits=c(0, 175))+ 
    scale_x_continuous(breaks = seq(0, 112000000000, by = 2000000))
}

# Display the combined plot
chr12_snp <- chr12_snp +
  ggtitle("SNP density on chromosome 12")

chr12_snp




chr15_snp <- ggplot()

# Iterate through the files and add each plot to the combined plot
for (file_name in file_names) {
  data <- read.table(paste0(dir_path, file_name), header = TRUE)
  file_title <- file_path_sans_ext(basename(file_name))
  
  data <- data[data$CHROM == 15, ]
  # Add the plot for the current file to the combined plot
  chr15_snp <- chr15_snp +
    geom_line(data = data, aes(x = BIN_START, y = SNP_COUNT, colour=..y..), alpha = 0.4) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+
    scale_colour_gradientn(colors=c("gray", "green", "yellow","red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red", "red2","red2", "red2", "red2", "red2", "red4", "red4", "red4","red4", "red4", "red4", "red4", "red4", "red4", "red4", "red4"),
                           limits=c(0, 80))+ 
    scale_x_continuous(breaks = seq(0, 112000000000, by = 4000000))
}

# Display the combined plot
chr15_snp <- chr15_snp +
  ggtitle("SNP density on chromosome 15")

chr15_snp





chr20_snp <- ggplot()

# Iterate through the files and add each plot to the combined plot
for (file_name in file_names) {
  data <- read.table(paste0(dir_path, file_name), header = TRUE)
  file_title <- file_path_sans_ext(basename(file_name))
  
  data <- data[data$CHROM == 20, ]
  # Add the plot for the current file to the combined plot
  chr20_snp <- chr20_snp +
    geom_line(data = data, aes(x = BIN_START, y = SNP_COUNT, colour=..y..), alpha = 0.4) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+
    scale_colour_gradientn(colors=c("gray", "green", "yellow","red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red", "red2","red2", "red2", "red2", "red2", "red4", "red4", "red4","red4", "red4", "red4", "red4", "red4", "red4", "red4", "red4"),
                           limits=c(0, 260))+ 
    scale_x_continuous(breaks = seq(0, 112000000000, by = 2000000))
}

# Display the combined plot
chr20_snp <- chr20_snp +
  ggtitle("SNP density on chromosome 20")

chr20_snp






chrX_snp <- ggplot()

# Iterate through the files and add each plot to the combined plot
for (file_name in file_names) {
  data <- read.table(paste0(dir_path, file_name), header = TRUE)
  file_title <- file_path_sans_ext(basename(file_name))
  
  data <- data[data$CHROM == "X", ]
  # Add the plot for the current file to the combined plot
  chrX_snp <- chrX_snp +
    geom_line(data = data, aes(x = BIN_START, y = SNP_COUNT, colour=..y..), alpha = 0.4) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")+
    scale_colour_gradientn(colors=c("gray", "green", "yellow","red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red", "red2","red2", "red2", "red2", "red2", "red4", "red4", "red4","red4", "red4", "red4", "red4", "red4", "red4", "red4", "red4"),
                           limits=c(0, 260))+ 
    scale_x_continuous(breaks = seq(0, 112000000000, by = 2000000))
}

# Display the combined plot
chrX_snp <- chrX_snp +
  ggtitle("SNP density on chromosome X")

chrX_snp





################################################################################

setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/snp_density/")
dir_path <- "C:/Users/Bedmar/Desktop/IDIBAPS/RADON-MINERS/snp_density/"
pattern <- "*snpden"

# Get file names
file_names <- list.files(dir_path, pattern = pattern)

#RADON-MINERS
for (file_name in file_names) {
  density_data <- read.table(paste0(dir_path, file_name), header = TRUE)
  file_title <- file_path_sans_ext(basename(file_name))
  density_data$CHROM <- gsub("chr", "", density_data$CHROM)
  
  plot <- ShowSNPDensityPlot(density_data, binSize = 1e6, chromSet = c(1:22)) + ggtitle(file_title)
  
  png(file = paste0(dir_path, gsub(pattern = pattern, replacement = ".png", x = file_name)),
      width = 850, height = 1080, res = 150)
  print(plot)
  dev.off()
}




################################################################################MEAN PER GROUP





################################################################################
library(tidyr)
library(ggplot2)
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/snp_density/")
df<-read.csv(file = "snp_density_CHR_MEAN_PER_GROUP.CSV", sep=",", head=T)
df<-df[-c(21,23),]
# Assuming your data is stored in a data frame called df

df[, -1] <- lapply(df[, -1], function(x) as.numeric(gsub(",", ".", x)))
df$Mean <- rowMeans(df[, -1], na.rm = TRUE)
df

# Convert the data from wide to long format
df_long <- pivot_longer(df, cols = c(D1, D3, D6, D12, Pr), names_to = "group", values_to = "mean")

# Convert the mean column from character to numeric
df_long$mean <- as.numeric(gsub(",", ".", df_long$mean))


group_order <- c("D1", "D3", "D6", "D12", "Pr")
df_long$group <- factor(df_long$group, levels = group_order)

chr_order <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X")
df_long$chr <- factor(df_long$chr, levels = chr_order)


head(df_long)
# Create the plot
ggplot(df_long, aes(x = chr, y = mean, group = group, color = group)) +
  geom_line() +
  geom_point() +
  geom_line(aes(y = Mean), linetype = "dashed", color = "black") +  # Add dashed line for Mean
  labs(x = "Chromosome", y = "Mean SNP density count") +
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(), legend.title = element_blank())+
  scale_colour_viridis_d()



################################################################################
#ANOVA

df<-read.csv(file = "snp_density_CHR_MEAN_PER_GROUP.CSV", sep=",", head=T)
df<-df[-c(21,22,23),]

df <- data.frame(lapply(df, function(x) as.numeric(gsub(",", ".", x))))

# Perform one-way ANOVA
head(df)
anova_result <- aov(D1 ~ chr + D3 + D6 + D12 + Pr, data = df)
summary(anova_result)


library(dplyr)

# Reshape the data
melted_df <- melt(df, id.vars = "chr", variable.name = "group", value.name = "snp_density")

# Perform the Kruskal-Wallis test for each chromosome
result <- melted_df %>%
  group_by(chr) %>%
  summarize(p_value = kruskal.test(snp_density ~ group)$p.value)

# Print the results
print(result)










results <- lapply(df[-1], function(group) {
  aov_result <- aov(group ~ chr, data = df)
  summary(aov_result)
})

results





pattern <- "^RnD1-.*snpden$"
D1 <- list.files(dir_path, pattern = pattern)
D1_files <- c(D1)
pattern <- "^RnD3-.*snpden$"
D3 <- list.files(dir_path, pattern = pattern)
D3_files <- c(D3)
pattern <- "^RnD6-.*snpden$"
D6 <- list.files(dir_path, pattern = pattern)
D6_files <- c(D6)
pattern <- "^RnD12-.*snpden$"
D12 <- list.files(dir_path, pattern = pattern)
D12_files <- c(D12)
pattern <- "^RnPr.*snpden$"
Pr <- list.files(dir_path, pattern = pattern)
Pr_files <- c(Pr)

head(Pr_files)


# Function to read SNP density data file and extract relevant columns
read_snp_density <- function(file) {
  data <- read.table(file, header = TRUE)
  data <- data[, c("CHROM", "SNP_COUNT")]
  return(data)
}

# Combine the SNP density data for all groups into a single data frame
df <- rbind(
  do.call(rbind, lapply(D1_files, read_snp_density)),
  do.call(rbind, lapply(D3_files, read_snp_density)),
  do.call(rbind, lapply(D6_files, read_snp_density)),
  do.call(rbind, lapply(D12_files, read_snp_density)),
  do.call(rbind, lapply(Pr_files, read_snp_density))
)

# Perform the Kruskal-Wallis test
kruskal_result <- kruskal.test(SNP_COUNT ~ CHROM, data = df)

# Print the p-value
print(kruskal_result$p.value)






























# Reshape the data to long format
df_long <- df %>% pivot_longer(-chr, names_to = "Group", values_to = "Mean_SNP")

df_long$Group <- factor(df_long$Group, levels=c("D1", "D3", "D6", "D12", "Pr"))

ggplot(df_long, aes(x = Group, y = Mean_SNP, fill = Group)) +
  geom_boxplot(aes(alpha=0.5)) +
  geom_jitter()+
  labs(x = "Group", y = "Mean SNP") +
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none")
  


file_names <- list.files(dir_path, pattern = pattern)

D1<-read.table(file = "D1.snpden", sep="\t", head=T)
D3<-read.table(file = "D3.snpden", sep="\t", head=T)
D6<-read.table(file = "D6.snpden", sep="\t", head=T)
D12<-read.table(file = "D12.snpden", sep="\t", head=T)
Pr<-read.table(file = "Pr.snpden", sep="\t", head=T)



D1$CHROM <- factor(D1$CHROM, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y","MT"))
D3$CHROM <- factor(D3$CHROM, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y","MT"))
D6$CHROM <- factor(D6$CHROM, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y","MT"))
D12$CHROM <- factor(D12$CHROM, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y","MT"))
Pr$CHROM <- factor(Pr$CHROM, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y","MT"))


ggplot()+
  geom_violin(data=D1, aes(x=CHROM ,y = SNP_COUNT, color="red"))+
  geom_violin(data=D3, aes(x=CHROM ,y = SNP_COUNT, color="blue"))+
  geom_violin(data=D6, aes(x=CHROM ,y = SNP_COUNT, color="green"))+
  geom_violin(data=D12, aes(x=CHROM ,y = SNP_COUNT, color="yellow"))+
  geom_violin(data=Pr, aes(x=CHROM ,y = SNP_COUNT, color="black"))+
  geom_violin(data=Pr, aes(x=CHROM ,y = SNP_COUNT, color="black"))+
  scale_color_brewer(palette = "Dark2", labels = c("D1", "D3", "D6", "D12", "Pr"))+
  theme_minimal()
                
              
              
              