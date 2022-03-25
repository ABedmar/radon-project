setwd('/Users/alexb/OneDrive/Escritorio/IDIBAPS')

library("limma")
library("edgeR")

counts <- read.csv(file = "Matrix.csv", sep = ";", header = T)

rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- data.matrix(counts)

dge <- DGEList(counts=counts)

design <- model.matrix(~ 0+factor(colnames(counts)))

keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

v <- voom(dge, design, plot=TRUE)

write.table(v$E,file = "GSEA_matrix_normExpr_Rsubread_correct.txt", sep = "\t")

boxplot(counts)
boxplot(v$E)
