# radon-project
Retrospective analyis of radon induced lung cancer rats


## Pipeline
### Quality control
[MultiQC](https://github.com/ewels/MultiQC) is a tool used to create a single report with interactive plots for multiple bioinformatics analyses across many samples. We can analyze our whole set of samples, included the trimmed ones with this tool.
You can use MultiQC by navigating to the desired analysis directory and running the tool:
```
multiqc .
```
The report is created in `multiqc_report.html`

### Trimming
Using trim galore...

### Normalization
```
library("limma")
library("edgeR")
counts <- read.csv(file = "matrix.csv", sep = ";", header = T)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- data.matrix(counts)

dge <- DGEList(counts=counts)
design <- model.matrix(~ 0+factor(colnames(counts)))
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)
topTable(fit, coef=ncol(design))
write.table(v$E,file = "counts_normalized.txt", sep = "\t") 
```

## Analysis

### Deconvolution


### Gene expression analyis


### TMB

Generate maf files
Mergethe maf files and use maftools to create noseque...

```
#module load gcc/8.1.0 pcre2/10.35 R/4.0.3

library(maftools)
files=dir(pattern=".maf$")
maf=merge_mafs(files)

laml.mutload = tcgaCompare(maf = maf, cohortName = 'Radon-LAML', logscale = TRUE, capture_size = 50)

rb=(laml.mutload$mutation_burden_perSample)
rb= rb[rb$cohort == "Radon-LAML",]
```

rn contains a table with the maf samples and all of their TMB.
