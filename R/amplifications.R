setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/cnv/amplificaitons/")
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(tibble)



gene_drivers<-c("Atrx","Nf1","Nras","Met","Ctnnb1","Erbb3","Erbb4","Fgfr2","Keap1","Akt1","Atm","Fgfr3","Fgfr4","Mapk1","Pten","Alk","Braf","Egfr","Erbb2","Fgfr1","Ntrk3","Ret","Ros1","Stk11","Hras","Kras","Map2k1","Nrg1","Ntrk1","Ntrk2","Pik3ca","Tp53")

# Create a list to store data from each file
file_names <- list.files(pattern = "*.txt")

for (file in file_names) {
  # Read the lines of the file
  lines <- readLines(file)
  # Initialize an empty list to store the data
  data_list <- list()
  # Determine the maximum number of columns in the table
  max_cols <- 0
  for (line in lines) {
    elements <- strsplit(line, "\t")[[1]]
    if (length(elements) > max_cols) {
      max_cols <- length(elements)
    }
  }
  # Process each line and pad with NA values if necessary
  for (line in lines) {
    elements <- strsplit(line, "\t")[[1]]
    padded_elements <- c(elements, rep(NA, max_cols - length(elements)))
    data_list <- c(data_list, list(padded_elements))
  }
  # Create a data frame from the list
  data <- as.data.frame(do.call(rbind, data_list))
  data<-data[-1,8:ncol(data)]
  data_melt<-melt(data, id=("V8"))
  
  df<-data_melt[,-2]
  
  df_matrix <- as.matrix(df)
  
  name <- gsub("_amp.txt$", "", basename(file))
  write.table(df_matrix, file = paste0(name, ".tsv"), sep = "\t")
}


file_names <- list.files(pattern = "*.tsv")
data_frames <- list()

# Loop through each file, read it, and store it in the data_frames list
for (file in file_names) {
  data <- read.table(file, header = TRUE, sep = "\t")
  data$Filename <- basename(file)  # Add a new column for the filename
  data_frames[[file]] <- data
}

# Combine all data frames into a single data frame
combined_data <- bind_rows(data_frames)

rownames(combined_data) <- rep(1:nrow(combined_data), 1)
combined_data$V8 <- gsub("^SVTYPE=", "", combined_data$V8)
combined_data$Filename <- gsub("\\.tsv", "", combined_data$Filename)
colnames(combined_data) <- c("type", "driver", "sample")
head(combined_data)


combined_data<-as.data.frame(combined_data)

result_df <- pivot_wider(
  data = combined_data,
  names_from = sample,
  values_from = type,
  values_fn = list(type = toString)
)


result_df<-as.data.frame(result_df)
result_df$driver <- replace(result_df$driver, is.na(result_df$driver), "NA")

rownames(result_df) <- result_df$driver

result_df <- result_df[, -1]

result_df<-result_df[-6,]
head(result_df)

write.table(result_df, file = "amplifications_per_sample.tsv", sep = "\t")
