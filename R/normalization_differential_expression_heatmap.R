setwd('/Users/alexb/OneDrive/Escritorio/10-05-2022_Cibersort')

library("limma")
library("edgeR")

counts <- read.csv(file = "matrix_controls_no_norm.csv", sep = ";", header = T)

rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- data.matrix(counts)

dge <- DGEList(counts=counts)
design2 <- model.matrix(~ 0 + factor(colnames(counts)))

keep <- filterByExpr(dge, design2)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design2, plot=TRUE)
#write.table(v$E,file = "matrix_controls_norm.txt", sep = "\t")

##Plot unnormalized vs normalized data
par(mfrow=c(1,2))
boxplot(counts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(counts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")

################################################################################
#Testing for differential expression
targets=factor(c("D1","D1","D1","D1","D1","D1","D1","D1","D1","D1","D3","D3","D3","D3","D3","D3","D3","D3","D6","D6","D6","D6","D6","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","Pr","Pr","Pr","Pr","Ctrl","Ctrl"))

f <- factor(targets, levels=c("D1","D3","D6","D12","Pr","Ctrl"))
design <- model.matrix(~0+f)
colnames(design) <- c("D1","D3","D6","D12","Pr","Ctrl")

fit <- lmFit(v$E, design)
  contrast.matrix <- makeContrasts(D1-Ctrl,D3-Ctrl,D6-Ctrl,D12-Ctrl,Pr-Ctrl,D1-D3,D1-D6,D1-D12,D1-Pr,D3-D6,D3-D12,D3-Pr,D6-D12,D6-Pr,D12-Pr, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
toptable.fit<-topTable(fit2, n=nrow(v$E))
results <- decideTests(fit2)
vennDiagram(results, 
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))

#Venn diagram all vs ctrl
#Venn diagram interesting comparisons

################################################################################
# Heatmap with differentialy expressed genes
library(devtools)
library(ComplexHeatmap)
library("RColorBrewer")
library(limma)
library(circlize)
library(ggfortify)
library(dplyr)

de=as.data.frame(results)
de=de %>% filter_all(any_vars(. %in% c(1,-1)))
de_genes=rownames(de)

normalized_mat=v$E+abs(min(v$E))
diff_exp_genes=subset(normalized_mat, rownames(normalized_mat) %in% de_genes)
nrow(diff_exp_genes)


mat=data.matrix(diff_exp_genes)
rownames(mat)=rownames(diff_exp_genes)
Sample_Group=c("D1","D1","D1","D1","D1","D1","D1","D1","D1","D1","D3","D3","D3","D3","D3","D3","D3","D3","D6","D6","D6","D6","D6","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","D12","Pr","Pr","Pr","Pr","Ctrl","Ctrl")
#Histology=c("Squamous","Squamous","Adenocarcinoma","Adenocarcinoma","Squamous","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Malign_tumor","Malign_tumor","Malign_tumor","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Osteosarcoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Osteosarcoma","Mixte","Malign_tumor","Adenocarcinoma","Malign_tumor","Adenocarcinoma","Adenocarcinoma","Squamous","Adenocarcinoma","Malign_tumor","Malign_tumor","Adenocarcinoma","Adenocarcinoma","Malign_tumor","Adenocarcinoma","Adenocarcinoma","Adenocarcinoma","Ctrl","Ctrl")

col=brewer.pal(length(unique(Sample_Group)),"RdYlBu")
#col2=brewer.pal(length(unique(Histology)),"PiYG")

names(col)=unique(Sample_Group)
col_fun = colorRamp2(c(0,3,6,9,12,15,18), brewer.pal(n = 7, name = 'RdBu'))
#names(col)=unique(Histology)
#col_fun = colorRamp2(c(0,20,40,60,80,100), brewer.pal(n = 6, name = 'Greens'))

Heatmap(mat, name = "expression", row_km = 4, column_km = 6, col = col_fun,show_row_names=F,
        top_annotation = HeatmapAnnotation(Sample_Group = Sample_Group,col=list(Sample_Group=col)))

################################################################################
#Correlation
round(cor(counts), digits = 3)

corrplot2 <- function(data,
                      method = "pearson",
                      sig.level = 0.05,
                      order = "original",
                      diag = FALSE,
                      type = "upper",
                      tl.srt = 90,
                      number.font = 1,
                      number.cex = 1,
                      mar = c(0, 0, 0, 0)) {
  library(corrplot)
  data_incomplete <- data
  data <- data[complete.cases(data), ]
  mat <- cor(data, method = method)
  cor.mtest <- function(mat, method) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = method)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(data, method = method)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(mat,
           method = "color", col = col(200), number.font = number.font,
           mar = mar, number.cex = number.cex,
           type = type, order = order,
           addCoef.col = "black", # add correlation coefficient
           tl.col = "black", tl.srt = tl.srt, # rotation of text labels
           # combine with significance level
           p.mat = p.mat, sig.level = sig.level, insig = "blank",
           # hide correlation coefficiens on the diagonal
           diag = diag
  )
}
corrplot2(
  data = counts,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75
)

corrplot2(
  data = v$E,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75
)

matrix_controls_norm<-v$E
matrix_controls_norm=matrix_controls_norm+abs(min(v$E))
matrix_controls_norm
write.table(matrix_controls_norm,file = "matrix_controls_norm_summed.txt", sep = "\t")
