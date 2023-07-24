setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS")

library(STRINGdb)
library(DESeq2)
library("limma")
library("edgeR")
library(dplyr)
library(tibble)
library(biomaRt)
library(RColorBrewer)


setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS")

counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)
counts<- counts[,-1:-2]
head(counts)
counts <- data.matrix(counts)
dge <- DGEList(counts=counts)

#################design cumulative
cumulative_dose <- factor(colnames(counts))
# Create a vector of the same length as cumulative_dose, initialized to "high_cumulative"
dose_type <- rep("high_cumulative", length(cumulative_dose))
# Set the dose type to "low_comultative" for levels starting with "RnD6"
dose_type[grep("^RnD6", cumulative_dose)] <- "low_comultative"
dose_type[grep("^X445", cumulative_dose)] <- "CONTROL"
cumulative_dose<-factor(dose_type)
design_cumulative <- model.matrix(~ 0 + cumulative_dose)

colnames(design_cumulative) <- c("high", "low")



keep <- filterByExpr(dge, design_cumulative)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

fit <- lmFit(dge$counts, design_cumulative)

contrast.matrix <- makeContrasts(low-high,
                                 levels=design_cumulative)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
genes_cumulative<-topTable(fit2, coef=1, adjust="BH", number = nrow(fit2))



genes_cumulative$diffexpressed <- "NO"
genes_cumulative$diffexpressed[genes_cumulative$logFC > 0.6 & genes_cumulative$P.Value < 0.05] <- "UP"
genes_cumulative$diffexpressed[genes_cumulative$logFC < -0.6 & genes_cumulative$P.Value < 0.05] <- "DOWN"

genes_cumulative$delabel <- NA
non_missing_genes <- !is.na(rownames(genes_cumulative))
genes_cumulative$delabel[genes_cumulative$diffexpressed != "NO" & non_missing_genes] <- rownames(genes_cumulative)[genes_cumulative$diffexpressed != "NO" & non_missing_genes]


ggplot(data=genes_cumulative, aes(x=logFC, y=-log10(P.Value), col = genes_cumulative$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("Cumulative high vs. low")+
  xlim(-1500, 550)



genes_cumulative<-subset(genes_cumulative, P.Value < 0.05)
dim(genes_cumulative)

write.table(genes)

#################design D1+D12 vs. D6
counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)
counts<- counts[,-c(1,2,24,25,26,27,28,29,30,35,36,37,38)]
d1d12_d6_dose <- factor(colnames(counts))
counts <- data.matrix(counts)
dge <- DGEList(counts=counts)

# Create a vector of the same length as d1d12_d6_dose, initialized to "high"
dose_type <- rep("high", length(d1d12_d6_dose))
# Set the dose type to "low" for levels starting with "RnD6"
dose_type[grep("^RnD6", d1d12_d6_dose)] <- "low"

d1d12_d6_dose<-factor(dose_type)
design_d1d12_d6_dose <- model.matrix(~ 0 + d1d12_d6_dose)

colnames(design_d1d12_d6_dose) <- c("high", "low")



keep <- filterByExpr(dge, design_d1d12_d6_dose)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

fit <- lmFit(dge$counts, design_d1d12_d6_dose)

contrast.matrix <- makeContrasts(low-high,
                                 levels=design_d1d12_d6_dose)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
genes_d1d12_d6_dose<-topTable(fit2, coef=1, adjust="BH", number = nrow(fit2))
genes_d1d12_d6_dose<-subset(genes_d1d12_d6_dose, P.Value < 0.05)
dim(genes_d1d12_d6_dose)


#################design PAEC
counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)
counts<- counts[,-1:-2]
head(counts)
PAEC_dose <- factor(colnames(counts))
dose_type <- rep("high_PAEC", length(PAEC_dose))
dose_type[grep("^RnD3", PAEC_dose)] <- "low_PAEC"
dose_type[grep("^RnD6", PAEC_dose)] <- "low_PAEC"
dose_type[grep("^RnD12", PAEC_dose)] <- "low_PAEC"
dose_type[grep("^X445", PAEC_dose)] <- "CONTROL"

PAEC_dose<-factor(dose_type)
design_PAEC <- model.matrix(~ 0 + PAEC_dose)


colnames(design_PAEC) <- c("high", "low")

keep <- filterByExpr(dge, design_PAEC)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

fit <- lmFit(dge$counts, design_PAEC)

contrast.matrix <- makeContrasts(low-high,
                                 levels=design_PAEC)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
genes_PAEC<-topTable(fit2, coef=1, adjust="BH", number = nrow(fit2))



genes_PAEC$diffexpressed <- "NO"
genes_PAEC$diffexpressed[genes_PAEC$logFC > 0.6 & genes_PAEC$P.Value < 0.05] <- "UP"
genes_PAEC$diffexpressed[genes_PAEC$logFC < -0.6 & genes_PAEC$P.Value < 0.05] <- "DOWN"

genes_PAEC$delabel <- NA
non_missing_genes <- !is.na(rownames(genes_PAEC))
genes_PAEC$delabel[genes_PAEC$diffexpressed != "NO" & non_missing_genes] <- rownames(genes_PAEC)[genes_PAEC$diffexpressed != "NO" & non_missing_genes]


ggplot(data=genes_PAEC, aes(x=logFC, y=-log10(P.Value), col = genes_PAEC$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("PAEC high(D1+Pr) vs. low(D3+D6+D12)")+
  xlim(-3000, 2000)

ggplot(data=genes_PAEC, aes(x=logFC, y=-log10(P.Value), col = genes_PAEC$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("PAEC high(D1+Pr) vs. low(D3+D6+D12)")+
  xlim(-10, 10)




genes_PAEC<-subset(genes_PAEC, P.Value < 0.05)
dim(genes_PAEC)


diff <- anti_join(genes_cumulative, genes_PAEC)

################################################################################
colnames(genes_cumulative) <- c("logFC", "AveExpr", "t", "pvalue", "adj.P.Val", "B")
genes_cumulative <- rownames_to_column(genes_cumulative, var = "gene")

colnames(genes_PAEC) <- c("logFC", "AveExpr", "t", "pvalue", "adj.P.Val", "B")
genes_PAEC <- rownames_to_column(genes_PAEC, var = "gene")

colnames(genes_d1d12_d6_dose) <- c("logFC", "AveExpr", "t", "pvalue", "adj.P.Val", "B")
genes_d1d12_d6_dose <- rownames_to_column(genes_d1d12_d6_dose, var = "gene")



genes_cumulative_pval_0_05<-genes_cumulative[,c(5,2,1)]
genes_PAEC_pval_0_05<-genes_PAEC[,c(5,2,1)]
genes_d1d12_d6_dose_pval_0_05<-genes_d1d12_d6_dose[,c(5,2,1)]

################################################################################
# set up Ensembl Mart database
ensembl <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
# define function to convert transcript IDs to gene IDs
convert_transcript_to_gene <- function(transcript_ids){
  # get gene IDs from Ensembl Mart
  gene_ids <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), 
                    filters = "ensembl_transcript_id",
                    values = transcript_ids,
                    mart = ensembl)
  # return gene IDs as a named vector
  setNames(gene_ids$ensembl_gene_id, gene_ids$ensembl_transcript_id)
}

# create a named vector of gene IDs for each transcript ID in the dataframe
gene_ids_PAEC <- convert_transcript_to_gene(genes_PAEC_pval_0_05$gene)

# add gene IDs as a column to the dataframe
genes_PAEC_pval_0_05$ensembl_gene_id <- gene_ids_PAEC[genes_PAEC_pval_0_05$gene]

# remove transcript ID column
genes_PAEC_pval_0_05$gene <- NULL



# create a named vector of gene IDs for each transcript ID in the dataframe
gene_ids_cumulative <- convert_transcript_to_gene(genes_cumulative_pval_0_05$gene)
# add gene IDs as a column to the dataframe
genes_cumulative_pval_0_05$ensembl_gene_id <- gene_ids_cumulative[genes_cumulative_pval_0_05$gene]
# remove transcript ID column
genes_cumulative_pval_0_05$gene <- NULL


# now with d1d12_d6
gene_ids_d1d12 <- convert_transcript_to_gene(genes_d1d12_d6_dose_pval_0_05$gene)
genes_d1d12_d6_dose_pval_0_05$ensembl_gene_id <- gene_ids_d1d12[genes_d1d12_d6_dose_pval_0_05$gene]
genes_d1d12_d6_dose_pval_0_05$gene <- NULL


colnames(genes_cumulative_pval_0_05)<-c("pvalue", "logFC", "gene")
colnames(genes_PAEC_pval_0_05)<-c("pvalue", "logFC", "gene")
colnames(genes_d1d12_d6_dose_pval_0_05)<-c("pvalue", "logFC", "gene")

head(genes_cumulative_pval_0_05)
head(genes_PAEC_pval_0_05)
head(genes_d1d12_d6_dose_pval_0_05)


write.table(genes_cumulative_pval_0_05, file = "genes_cumulative_pval_0_05.tsv", sep = "\t")

write.table(genes_PAEC_pval_0_05, file = "genes_PAEC_pval_0_05.tsv", sep = "\t")

# protein protein interaction
####################################################
####################################################

string_db<-STRINGdb$new(version)



#Map the gene ID’s with STRING database
string_db <- STRINGdb$new(version="11.5", species=10116 , score_threshold=400, network_type="full", input_directory="")
####################################################
####################################################CUMULTAIVE PPI
####################################################
mapped <- string_db$map(genes_cumulative_pval_0_05, "gene", removeUnmappedRows = TRUE)

hits <- mapped$STRING_id

#Generate clusters on the PPI network
clustersList <- string_db$get_clusters(hits, algorithm="fastgreedy")

string_db$plot_network(clustersList[[1]])

#par(mfrow=c(4,3))
#for(i in seq(1:12)){
#  string_db$plot_network(clustersList[[i]])
#}


#PAYLOAD MECHANISM
mapped_pval05 <- string_db$add_diff_exp_color(subset(mapped, pvalue<0.05),
                                              logFcColStr="logFC")


#ADD COLOR GRADIENT TO INDICATE OVER OR UNDER REPRESENTATION
# Define function to convert RGB to HEX
rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)

# Define the color gradient for negative values
neg_colors <- c("#ffffff", "#ffbaba", "#ff7a7a", "#ff4545", "#ff0000")

# Define the color gradient for positive values
pos_colors <- c("#ffffff", "#b4ffb3", "#82ff80", "#4aff47", "#04ff00")

# Define a function to get the corresponding color for a logFC value
get_color <- function(logFC) {
  if (logFC >= 0) {
    # Get the index of the corresponding color in the positive gradient
    color_index <- floor((logFC / max(mapped_pval05$logFC)) * (length(pos_colors) - 1)) + 1
    # Get the corresponding color
    color <- pos_colors[color_index]
  } else {
    # Get the index of the corresponding color in the negative gradient
    color_index <- floor((logFC / min(mapped_pval05$logFC)) * (length(neg_colors) - 1)) + 1
    # Get the corresponding color
    color <- neg_colors[color_index]
  }
  # Convert the color to RGB and then to hex
  rgb <- col2rgb(color)
  hex <- rgb2hex(rgb[1], rgb[2], rgb[3])
  return(hex)
}


# Apply the function to the logFC column and assign the resulting colors to the color column
mapped_pval05$color <- sapply(mapped_pval05$logFC, get_color)


payload_id <- string_db$post_payload(mapped_pval05$STRING_id,
                                     colors=mapped_pval05$color)

#string_db$plot_network(clustersList[[1]], payload_id=payload_id)

par(mfrow=c(2, 1))
for(i in seq(1:2)){
  string_db$plot_network(clustersList[[i]], payload_id=payload_id)
}
####################################################
####################################################PAEC PPI
####################################################
mapped <- string_db$map(genes_PAEC_pval_0_05, "gene", removeUnmappedRows = TRUE)

hits <- mapped$STRING_id

#Generate clusters on the PPI network
clustersList <- string_db$get_clusters(hits, algorithm="fastgreedy")
string_db$plot_network(clustersList[[1]])

#par(mfrow=c(4,3))
#for(i in seq(1:12)){
#  string_db$plot_network(clustersList[[i]])
#}


#PAYLOAD MECHANISM
mapped_pval05 <- string_db$add_diff_exp_color(subset(mapped, pvalue<0.05),
                                              logFcColStr="logFC")


# Apply the function to the logFC column and assign the resulting colors to the color column
mapped_pval05$color <- sapply(mapped_pval05$logFC, get_color)

##########


payload_id <- string_db$post_payload(mapped_pval05$STRING_id,
                                     colors=mapped_pval05$color)

par(mfrow=c(1,1))
string_db$plot_network(clustersList[[2]], payload_id=payload_id)

par(mfrow=c(3,5))
for(i in seq(1:15)){
  string_db$plot_network(clustersList[[i]], payload_id=payload_id)
}




####################################################
####################################################
####################################################




####################################################
####################################################d1d12_d6 PPI
####################################################
mapped <- string_db$map(genes_d1d12_d6_dose_pval_0_05, "gene", removeUnmappedRows = TRUE)

hits <- mapped$STRING_id


#Generate clusters on the PPI network
clustersList <- string_db$get_clusters(hits, algorithm="fastgreedy")
string_db$plot_network(clustersList[[1]])

#par(mfrow=c(4,3))
#for(i in seq(1:12)){
#  string_db$plot_network(clustersList[[i]])
#}


#PAYLOAD MECHANISM
mapped_pval05 <- string_db$add_diff_exp_color(subset(mapped, pvalue<0.05),
                                              logFcColStr="logFC")


# Apply the function to the logFC column and assign the resulting colors to the color column
mapped_pval05$color <- sapply(mapped_pval05$logFC, get_color)

##########


payload_id <- string_db$post_payload(mapped_pval05$STRING_id,
                                     colors=mapped_pval05$color)

par(mfrow=c(1,1))
string_db$plot_network(clustersList[[2]], payload_id=payload_id)

par(mfrow=c(1,2))
for(i in seq(1:2)){
  string_db$plot_network(clustersList[[i]], payload_id=payload_id)
}




####################################################
####################################################
####################################################





#GO Biological Process
BP <- string_db$get_enrichment(clustersList[[1]], category='Process', methodMT='fdr')
BP$genes <- paste0(BP$number_of_genes,'/', BP$number_of_genes_in_background)
BP <- BP[,c('category','term','description','genes', 'p_value', 'fdr')]
rownames(BP) <- NULL # "reset" row indexes
head(BP)
#Cancer hallmarks
GO_cancer <- c('GO:0030307', # Positive regulation of cell growth
               'GO:0008284', # Positive regulation of cell proliferation
               'GO:0045787', # Positivie regulation of cell cycle
               'GO:0045786', # negative regulation of cell cycle
               'GO:0008283', # Cell Proliferation
               'GO:0007049', # Cell Cycle
               'GO:0051301', # Cell Division
               'GO:0051781', # Positive Regulation of Cell Division
               'GO:0009968', # Negative regulation of signal transduction
               'GO:0030308', # Negative regulation of cell growth
               'GO:0008285', # Negative regulation of cell proliferation
               'GO:0043069', # Negative regualtion of programmed cell death
               'GO:0012501', # programmed cell death
               'GO:0043067', # regulation of programmed cell death
               'GO:0001302', # Replicative cell aging           
               'GO:0090398', # Cellular senescence          
               'GO:0000723', # telomere maintenance         
               'GO:0032204', # regulation of telomere maintenance           
               'GO:2000772', # Regulation of cellular senescence            
               'GO:0032200', # telomere Organization
               'GO:0045766', # Positive regulation of angiogenesis
               'GO:0001570', # Vasculogenesis           
               'GO:0001525', # angiogenesis 
               'GO:0007162', # Negative regulation of cell adhesion         
               'GO:0001837', # Epithelial to mesenchymal transition         
               'GO:0016477', # Cell migration           
               'GO:0007155', # Cell adhesion            
               'GO:0034330', # Cell junction Organization           
               'GO:0030030', # Cell projection Organization
               'GO:0030155', # Regulation of Cell Adhesion
               'GO:0045005', # Maintenance of fidelity involved in DNA-dependent DNA replication            
               'GO:0006281', # DNA repair           
               'GO:0006282', # regulation of DNA repair         
               'GO:0031570', # DNA integrity Checkpoint
               'GO:0002367', # Cytokine production involved in immune response          
               'GO:0050727', # regulation of inflammatory response          
               'GO:0006096', # Glycolysis
               'GO:0071456', # Cellular response to hypoxia
               'GO:0002837', # Regulation of immune response to tumor cell          
               'GO:0002418', # Immune response to tumor cells           
               'GO:0050776', # Regulation of immune response            
               'GO:0006955') # Immune Response          


BP[BP$term %in% GO_cancer,]




#KEGG
KEGG <- string_db$get_enrichment(clustersList[[1]], category='KEGG', methodMT='fdr')
KEGG$genes <- paste0(KEGG$number_of_genes,'/', KEGG$number_of_genes_in_background)
KEGG <- KEGG[,c('category','term','description','genes', 'p_value', 'fdr')]
rownames(KEGG) <- NULL # "reset" row indexes
head(KEGG)
#Cancer-related pathways
KEGG_cancer <- c('hsa04010', # MAPK signaling pathway
                 'hsa04020', # Calcium signaling pathway
                 'hsa04024', # cAMP signaling pathway
                 'hsa04060', # Cytokine-cytokine receptor interaction
                 'hsa04066', # HIF-1 signaling pathway
                 'hsa04110', # Cell cycle
                 'hsa04115', # p53 signaling pathway
                 'hsa04150', # mTOR signaling pathway
                 'hsa04151', # PI3K-Akt signaling pathway
                 'hsa04210', # Apoptosis
                 'hsa04310', # Wnt signaling pathway
                 'hsa04330', # Notch signaling pathway
                 'hsa04340', # Hedgehog signaling pathway
                 'hsa04350', # TGF-beta signaling pathway
                 'hsa04370', # VEGF signaling pathway
                 'hsa04510', # Focal adhesion
                 'hsa04512', # ECM-receptor interaction
                 'hsa04520', # Adherens junction
                 'hsa04630', # JAK-STAT signaling pathway
                 'hsa04915') # Estrogen signaling pathway


KEGG[KEGG$term %in% KEGG_cancer,]













counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)
head(counts)
counts <- data.matrix(counts)
dge <- DGEList(counts=counts)

#################design cumulative
cumulative_dose <- factor(colnames(counts))
# Create a vector of the same length as cumulative_dose, initialized to "high_cumulative"
dose_type <- rep("TUMOR", length(cumulative_dose))
# Set the dose type to "low_comultative" for levels starting with "RnD6"
dose_type[grep("^X445", cumulative_dose)] <- "CONTROL"
cumulative_dose<-factor(dose_type)
design_cumulative <- model.matrix(~ 0 + cumulative_dose)

colnames(design_cumulative) <- c("CONTROL", "TUMOR")


keep <- filterByExpr(dge, design_cumulative)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

fit <- lmFit(dge$counts, design_cumulative)

contrast.matrix <- makeContrasts(TUMOR-CONTROL,
                                 levels=design_cumulative)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
genes_cumulative<-topTable(fit2, coef=1, adjust="BH", number = nrow(fit2))



genes_cumulative$diffexpressed <- "NO"
genes_cumulative$diffexpressed[genes_cumulative$logFC > 0.6 & genes_cumulative$P.Value < 0.05] <- "UP"
genes_cumulative$diffexpressed[genes_cumulative$logFC < -0.6 & genes_cumulative$P.Value < 0.05] <- "DOWN"

genes_cumulative$delabel <- NA
non_missing_genes <- !is.na(rownames(genes_cumulative))
genes_cumulative$delabel[genes_cumulative$diffexpressed != "NO" & non_missing_genes] <- rownames(genes_cumulative)[genes_cumulative$diffexpressed != "NO" & non_missing_genes]


ggplot(data=genes_cumulative, aes(x=logFC, y=-log10(P.Value), col = genes_cumulative$diffexpressed, label=delabel)) +
  geom_point()+
  geom_text(size = 2.5, vjust=-0.8)+
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype= "dashed") +
  theme_minimal()+
  xlab("log2(Fold Change)")+
  ylab("-log10(pval)")+
  labs(color='Expression regulation')+
  ggtitle("CONTROL vs. TUMOR")+
  xlim(-100, 100)




