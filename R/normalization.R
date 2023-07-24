setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/deconvolution")
getwd()

library(ggplot2)
library("limma")
library("edgeR")
install.packages("devtools")
library("devtools")
library("ggradar")
library(reshape2)

counts <- read.csv(file = "matri_orth_uniq.CSV", sep = ",", header = T)
counts<- read.table(file="matrix.TXT", sep="\t", header=T)
counts<-counts[,-1]
head(counts)

rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- data.matrix(counts)

dge <- DGEList(counts=counts)
design2 <- model.matrix(~ 0 + factor(colnames(counts)))

keep <- filterByExpr(dge, design2)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design2, plot=TRUE)

min(v$E)
v$E<-v$E-min(v$E)

write.table(v$E,file = "matrix_norm_ENSRNOT.txt", sep = "\t")

##Plot unnormalized vs normalized data
par(mfrow=c(1,2),mar = c(7.2, 4.1, 4.1, 2.1))
boxplot(counts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(counts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")



################################################################################

#DECONVOLUTION PLOTS
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS/deconvolution/CIBERSORTx_Job27_output")

data<-read.table(file = "CIBERSORTxGEP_Job27_Fractions.TXT", sep = "\t", header = T)
ncol(data)
data<-data[,-c(24,25,26)]
data <- data[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 15, 16, 17, 20, 21, 14, 23, 22)]

colnames(data)
df<-melt(data, id.vars = unique(c("Mixture")))

unique(df$Mixture)
head(df)

ggplot(data=df, aes(x=factor(Mixture, level = c("RnD1.104","RnD1.132","RnD1.134","RnD1.162","RnD1.171","RnD1.189","RnD1.195","RnD1.207","RnD1.214","RnD1.36","RnD1.81","RnD1.98","RnD3.168","RnD3.203","RnD3.212","RnD3.232","RnD3.233","RnD3.71","RnD3.92","RnD6.132","RnD6.152","RnD6.179","RnD6.204","RnD12.145","RnD12.178","RnD12.184","RnD12.185.1","RnD12.185.2","RnD12.229","RnD12.231","RnD12.49","RnD12.51","RnPr109","RnPr140","RnPr146","RnPr214","44514","44541")),y=value, fill=variable))+
  geom_bar(stat="identity", position = "stack")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(element_blank())
  
control<-data[c(1,2),]
D1<-data[c(3,4,5,6,7,8,9,10,12,21,22,23),]
D3<-data[c(24,25,26,27,28,29,30),]
D6<-data[c(31,32,33,34),]
D12<-data[c(11,13,14,15,16,17,18,19,20),]
Pr<-data[c(35,36,37,38),]

lcols <- c("#B3B3B3","#B3B3B3","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#3b528b","#3b528b","#3b528b","#3b528b","#3b528b","#3b528b","#3b528b","#21908c","#21908c","#21908c","#21908c","#fde725","#fde725","#fde725","#fde725")
ctrllcols <- c("#B3B3B3","#B3B3B3")
D1lcols <- c("#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154","#440154")
D3lcols <- c("#3b528b","#3b528b","#3b528b","#3b528b","#3b528b","#3b528b","#3b528b")
D6lcols <- c("#21908c","#21908c","#21908c","#21908c")
D12lcols <- c("#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863","#5dc863")
Prlcols <- c("#fde725","#fde725","#fde725","#fde725")

ggradar(data,
        group.line.width = 0.3,
        group.point.size = 0.5,
        values.radar = c("0%","25%","50%"),
        font.radar = "calibri",
        grid.mid = 0.25,
        grid.max = 0.5,
        gridline.mid.colour = "gray",
        background.circle.transparency = 0.1,
        group.colours = lcols,
        legend.position = "none")

ctrl_radar<-ggradar(control,
        group.line.width = 0.3,
        group.point.size = 0.5,
        values.radar = c("0%","25%","50%"),
        font.radar = "calibri",
        grid.mid = 0.25,
        grid.max = 0.5,
        gridline.mid.colour = "gray",
        background.circle.transparency = 0.1,
        group.colours = ctrllcols,
        legend.position = "none",
        fill = T,
        fill.alpha = 0.1,
        axis.label.size = 3.5,
        grid.label.size = 5)

D1_radar<-ggradar(D1,
        group.line.width = 0.3,
        group.point.size = 0.5,
        values.radar = c("0%","25%","50%"),
        font.radar = "calibri",
        grid.mid = 0.25,
        grid.max = 0.5,
        gridline.mid.colour = "gray",
        background.circle.transparency = 0.1,
        group.colours = D1lcols,
        legend.position = "none",
        fill = T,
        fill.alpha = 0.1,
        axis.label.size = 3.5,
        grid.label.size = 5)

D3_radar<-ggradar(D3,
        group.line.width = 0.3,
        group.point.size = 0.5,
        values.radar = c("0%","25%","50%"),
        font.radar = "calibri",
        grid.mid = 0.25,
        grid.max = 0.5,
        gridline.mid.colour = "gray",
        background.circle.transparency = 0.1,
        group.colours = D3lcols,
        legend.position = "none",
        fill = T,
        fill.alpha = 0.1,
        axis.label.size = 3.5,
        grid.label.size = 5)

D6_radar<-ggradar(D6,
        group.line.width = 0.3,
        group.point.size = 0.5,
        values.radar = c("0%","25%","50%"),
        font.radar = "calibri",
        grid.mid = 0.25,
        grid.max = 0.5,
        gridline.mid.colour = "gray",
        background.circle.transparency = 0.1,
        group.colours = D6lcols,
        legend.position = "none",
        fill = T,
        fill.alpha = 0.1,
        axis.label.size = 3.5,
        grid.label.size = 5)

D12_radar<-ggradar(D12,
        group.line.width = 0.3,
        group.point.size = 0.5,
        values.radar = c("0%","25%","50%"),
        font.radar = "calibri",
        grid.mid = 0.25,
        grid.max = 0.5,
        gridline.mid.colour = "gray",
        background.circle.transparency = 0.1,
        group.colours = D12lcols,
        legend.position = "none",
        fill = T,
        fill.alpha = 0.1,
        axis.label.size = 3.5,
        grid.label.size = 5)

Pr_radar<-ggradar(Pr,
        group.line.width = 0.3,
        group.point.size = 0.5,
        values.radar = c("0%","25%","50%"),
        font.radar = "calibri",
        grid.mid = 0.25,
        grid.max = 0.5,
        gridline.mid.colour = "gray",
        background.circle.transparency = 0.1,
        group.colours = Prlcols,
        legend.position = "none",
        fill = T,
        fill.alpha = 0.1,
        axis.label.size = 3.5,
        grid.label.size = 5)


grid.arrange(ctrl_radar, D1_radar, D3_radar, D6_radar, D12_radar, Pr_radar, ncol=3, nrow=2)






-data<-read.table(file = "cell_fractions_TCELL.TXT", sep = "\t", header = T)

df<-melt(data, id.vars = unique(c("Mixture")))

ggplot(data=df, aes(x=factor(Mixture, level = c("RnD1_104","RnD1_132","RnD1_134","RnD1_162","RnD1_171","RnD1_189","RnD1_195","RnD1_207","RnD1_214","RnD1_36","RnD1_81","RnD1_98","RnD3_168","RnD3_203","RnD3_212","RnD3_232","RnD3_233","RnD3_71","RnD3_92","RnD6_132","RnD6_152","RnD6_179","RnD6_204","RnD12_145","RnD12_178","RnD12_184","RnD12_185.1","RnD12_185.2","RnD12_229","RnD12_231","RnD12_49","RnD12_51","RnPr109","RnPr140","RnPr146","RnPr214","44514","44541")),y=value, fill=variable))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#636363", "red", "orange"))+
  ylab("Cell proportion")+
  xlab(element_blank())+
  ggtitle("T Cell fractions radon-rats RNAseq deconvolution")




data<-read.table(file = "cell_fractions_BCELL.TXT", sep = "\t", header = T)

df<-melt(data, id.vars = unique(c("Mixture")))


ggplot(data=df, aes(x=factor(Mixture, level = c("RnD1_104","RnD1_132","RnD1_134","RnD1_162","RnD1_171","RnD1_189","RnD1_195","RnD1_207","RnD1_214","RnD1_36","RnD1_81","RnD1_98","RnD3_168","RnD3_203","RnD3_212","RnD3_232","RnD3_233","RnD3_71","RnD3_92","RnD6_132","RnD6_152","RnD6_179","RnD6_204","RnD12_145","RnD12_178","RnD12_184","RnD12_185.1","RnD12_185.2","RnD12_229","RnD12_231","RnD12_49","RnD12_51","RnPr109","RnPr140","RnPr146","RnPr214","44514","44541")),y=value, fill=variable))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#636363", "#69cadd"))+
  ylab("Cell proportion")+
  xlab(element_blank())+
  ggtitle("B Cell fractions radon-rats RNAseq deconvolution")




data<-read.table(file = "cell_fractions_grouped.TXT", sep = "\t", header = T)

df<-melt(data, id.vars = unique(c("Mixture")))

dec<-ggplot(data=df, aes(x=factor(Mixture, level = c("RnD1_104","RnD1_132","RnD1_134","RnD1_162","RnD1_171","RnD1_189","RnD1_195","RnD1_207","RnD1_214","RnD1_36","RnD1_81","RnD1_98","RnD3_168","RnD3_203","RnD3_212","RnD3_232","RnD3_233","RnD3_71","RnD3_92","RnD6_132","RnD6_152","RnD6_179","RnD6_204","RnD12_145","RnD12_178","RnD12_184","RnD12_185.1","RnD12_185.2","RnD12_229","RnD12_231","RnD12_49","RnD12_51","RnPr109","RnPr140","RnPr146","RnPr214","44514","44541")),y=value, fill=variable))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#f05f5f","#aaa5a5","#636363","#69cadd","#ffce5c","orange", "red"))+
  ylab("Cell proportion")+
  xlab(element_blank())+
  ggtitle("Cell fractions radon-rats")
dec

histology<-c("Adenocarcinoma","Adenocarcinoma","Squamous","Adenocarcinoma","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Squamous","Squamous","Squamous","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Other","Adenocarcinoma","Other","Adenocarcinoma","Adenocarcinoma")
Y<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","ç","ñ","aa","ab","ac","ad","ae","af","ag","at")
hist<-data.frame(histology,Y)
head(hist)

histo<-ggplot(data = hist,aes(x=Y,y=1, fill=histology))+
  geom_bar(stat="identity", position = "fill")+
  theme_void()+
  scale_fill_manual(values = c("#ef40d0","#40efc8","#ef8240"))+
  theme(aspect.ratio=1/37) 

histo  

grid.arrange(histo,dec,nrow=2)


#
data<-read.table(file = "cell_fractions_gouped_dose_groups.TXT", sep = "\t", header = T)

df<-melt(data, id.vars = unique(c("Mixture")))
ggplot(data=df, aes(x=factor(Mixture, level = c("D1","D3","D6","D12","Pr","Control")),y=value, fill=variable))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#f05f5f","#aaa5a5","#636363","#69cadd","#ffce5c","orange", "red"))+
  ylab("Cell proportion")+
  xlab(element_blank())+
  ggtitle("Cell fractions groups")





data<-read.table(file = "cell_fractions_grouped_dose.TXT", sep = "\t", header = T)

df<-melt(data, id.vars = unique(c("Mixture")))
ggplot(data=df, aes(x=factor(Mixture, level = c("Control","high","low")),y=value, fill=variable))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("#f05f5f","#aaa5a5","#636363","#69cadd","#ffce5c","orange", "red"))+
  ylab("Cell proportion")+
  xlab(element_blank())+
  ggtitle("Cell fractions total dose")










data<-read.table(file = "cell_fractions_grouped_not_that_much.TXT", sep = "\t", header = T)

df<-melt(data, id.vars = unique(c("Mixture")))

ggplot(data=df, aes(x=factor(Mixture, level = c("RnD1_104","RnD1_132","RnD1_134","RnD1_162","RnD1_171","RnD1_189","RnD1_195","RnD1_207","RnD1_214","RnD1_36","RnD1_81","RnD1_98","RnD3_168","RnD3_203","RnD3_212","RnD3_232","RnD3_233","RnD3_71","RnD3_92","RnD6_132","RnD6_152","RnD6_179","RnD6_204","RnD12_145","RnD12_178","RnD12_184","RnD12_185.1","RnD12_185.2","RnD12_229","RnD12_231","RnD12_49","RnD12_51","RnPr109","RnPr140","RnPr146","RnPr214","44514","44541")),y=value, fill=variable))+
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank())+
  scale_fill_manual(values = c("#582E5E", "#CF2F2F", "#dd7e6b", "#fce5cd", "#aaa5a5", "#666666", "#c8a276", "#ba7d36", "#8e7cc3", "#d5a6bd", "#6aa84f", "#b6d7a8", "#69cadd", "#ffce5c", "#ffa500"))+
  ylab("Cell proportion")+
  xlab(element_blank())+
  ggtitle("Cell fractions radon-rats deconvolution")













###################################################################################################################
setwd("C:/Users/Bedmar/Desktop/IDIBAPS/RADON-RATS")
# Load the TPM matrix from a file
tpm_matrix <- read.table("TPM.txt", header=TRUE, row.names=1)

# Convert TPM values to log2(TPM+1)
tpm_log2 <- log2(tpm_matrix + 1)

# Write the log2-transformed TPM matrix to a file
write.table(tpm_log2, file="tpm_matrix_log2.txt", sep="\t", quote=FALSE, col.names=NA)



# Load the log2(TPM+1) matrix from a file
tpm_log2 <- read.table("tpm_matrix_log2.txt", header=TRUE, row.names=1)

# Check that row names are present
rownames(tpm_log2)

# Round values to integers
tpm_int <- round(tpm_log2)

# Create the colData argument for DESeqDataSetFromMatrix
colData <- data.frame(row.names=colnames(tpm_int))

# Create DESeqDataSet object and normalize TPM values
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=tpm_int, colData=colData, design=~1)
dds <- DESeq(dds)
tpm_norm <- assay(rlog(dds))

# Write the normalized TPM matrix to a file
write.table(tpm_norm, file="tpm_matrix_norm.txt", sep="\t", quote=FALSE, col.names=NA)


########################################################################################################## NORMALIZE WITH LIMMA

counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)
counts<- counts[,-1:-2]
counts <- data.matrix(counts)
dge <- DGEList(counts=counts)

#################design COMULATIVE
comulative_dose <- factor(colnames(counts))
# Create a vector of the same length as comulative_dose, initialized to "high_comulative"
dose_type <- rep("high_comulative", length(comulative_dose))
# Set the dose type to "low_comultative" for levels starting with "RnD6"
dose_type[grep("^RnD6", comulative_dose)] <- "low_comultative"
dose_type[grep("^X445", comulative_dose)] <- "CONTROL"
comulative_dose<-factor(dose_type)
design_comulative <- model.matrix(~ 0 + comulative_dose)



keep <- filterByExpr(dge, design_comulative)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

fit <- lmFit(dge$counts, design_comulative)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res<-topTable(fit2, coef=1, adjust="BH", number = nrow(fit2))
head(res)


library(TxDb.Rnorvegicus.UCSC.rn6.ncbiRefSeq)


topTable()

results <- decideTests(fit2)




v <- voom(dge, design_comulative, plot=TRUE)

min(v$E)
v$E<-v$E-min(v$E)






write.table(v$E,file = "count_matrix_norm_WLM.txt", sep = "\t")

#################design PAEC
PAEC_dose <- factor(colnames(counts))
dose_type <- rep("high_PAEC", length(PAEC_dose))
dose_type[grep("^RnD3", PAEC_dose)] <- "low_PAEC"
dose_type[grep("^RnD6", PAEC_dose)] <- "low_PAEC"
dose_type[grep("^RnD12", PAEC_dose)] <- "low_PAEC"
dose_type[grep("^X445", PAEC_dose)] <- "CONTROL"

PAEC_dose<-factor(dose_type)
design_PAEC <- model.matrix(~ 0 + PAEC_dose)

keep <- filterByExpr(dge, design_PAEC)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design_PAEC, plot=TRUE)

min(v$E)
v$E<-v$E-min(v$E)
write.table(v$E,file = "count_matrix_norm_PAEC.txt", sep = "\t")

#################design both
WLM_and_PAEC <- factor(colnames(counts))
dose_type <- rep("high", length(WLM_and_PAEC))
dose_type[grep("^RnD3", WLM_and_PAEC)] <- "intermediate_high"
dose_type[grep("^RnD6", WLM_and_PAEC)] <- "low_cumulative_dose"
dose_type[grep("^RnD12", WLM_and_PAEC)] <- "intermediate_low"
dose_type[grep("^X445", WLM_and_PAEC)] <- "CONTROL"

WLM_and_PAEC<-factor(dose_type)
design_WLM_and_PAEC <- model.matrix(~ 0 + WLM_and_PAEC)


keep <- filterByExpr(dge, design_WLM_and_PAEC)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design_WLM_and_PAEC, plot=TRUE)

min(v$E)
v$E<-v$E-min(v$E)
write.table(v$E,file = "count_matrix_norm_WLM_and_PAEC.txt", sep = "\t")










################################################# VENN DIAGRAM
library(limma)
library(edgeR)
counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)
counts<- counts[,-1:-2]

# An appropriate design matrix can be created and a linear model fitted using
design_comulative <- model.matrix(~ 0 + comulative_dose)
head(design_comulative)
colnames(design_comulative) <- c("high", "low")

dge <- DGEList(counts=counts)
fit <- lmFit(dge$counts, design_comulative)

volcanoplot(fit)
contrast.matrix <- makeContrasts(low-high,
                                 levels=design_comulative)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)

vennDiagram(results, include=c("up","down"))

topTable(fit2, number=30)
title("COMMULATIVE")





################################################# VENN DIAGRAM
counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)

# An appropriate design matrix can be created and a linear model fitted using
design_PAEC <- model.matrix(~ 0 + PAEC_dose)
head(design_PAEC)
colnames(design_PAEC) <- c("high", "low")

dge <- DGEList(counts=counts)
fit <- lmFit(dge$counts, design_PAEC)

contrast.matrix <- makeContrasts(low-high,
                                 levels=design_PAEC)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)

vennDiagram(results, include=c("up","down"))
title("PAEC")

topTable(fit2, number=30)




################################################# VENN DIAGRAM
counts <- read.table("counts_matrix.txt", header=TRUE, row.names=1)

# An appropriate design matrix can be created and a linear model fitted using

WLM_and_PAEC <- factor(colnames(counts))
dose_type <- rep("high", length(WLM_and_PAEC))
dose_type[grep("^RnD3", WLM_and_PAEC)] <- "intermediate_high"
dose_type[grep("^RnD6", WLM_and_PAEC)] <- "low_cumulative_dose"
dose_type[grep("^RnD12", WLM_and_PAEC)] <- "intermediate_low"

WLM_and_PAEC<-factor(dose_type)
design_WLM_and_PAEC <- model.matrix(~ 0 + WLM_and_PAEC)

head(design_WLM_and_PAEC)
colnames(design_WLM_and_PAEC) <- c("high", "intermediate_high", "intermediate_low", "low_comulative_dose")

dge <- DGEList(counts=counts)
fit <- lmFit(dge$counts, design_WLM_and_PAEC)

contrast.matrix <- makeContrasts(high-intermediate_high, high-intermediate_low, high-low_comulative_dose, levels=design_WLM_and_PAEC)
contrast.matrix <- makeContrasts(intermediate_high-intermediate_low, intermediate_high-low_comulative_dose, intermediate_low-low_comulative_dose,
                                 levels=design_WLM_and_PAEC)


fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)

vennDiagram(results)

topTable(fit2, number=30)

