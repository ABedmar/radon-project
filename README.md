# radon-project
Retrospective analyis of radon induced lung cancer in rats. This analysis has the objective of describeing the tumor molecular profile related to the immune profile of radon-induced lung cancer vs. spontaneous lung cancer in rats and to characterize pathologically and transcriptomically all the TME in lung primary tumors with a particular interest in CD8+ T cells (inflamed-TME) in rats exposed or not to radon.

# Table of contents
- [Pipeline](https://github.com/ABedmar/radon-project/blob/main/README.md#pipeline)
    - [Quality control](https://github.com/ABedmar/radon-project/blob/main/README.md#Quality_control)
    - [Trimming](https://github.com/ABedmar/radon-project/blob/main/README.md#Trimming)
    - [Normalization](https://github.com/ABedmar/radon-project/blob/main/README.md#Normalization)
- [Analysis](https://github.com/ABedmar/radon-project/blob/main/README.md#Analysis)
    - [Deconvolution](https://github.com/ABedmar/radon-project/blob/main/README.md#Deconvolution)
    - [Gene expression analysis](https://github.com/ABedmar/radon-project/blob/main/README.md#Gene_expression_analysis)
    - [TMB](https://github.com/ABedmar/radon-project/blob/main/README.md#TMB)

# Pipeline
## Quality control
[MultiQC](https://github.com/ewels/MultiQC) is a tool used to create a single report with interactive plots for multiple bioinformatics analyses across many samples. We can analyze our whole set of samples, included the trimmed ones with this tool.
You can use MultiQC by navigating to the desired analysis directory and running the tool:
```
multiqc .
```
The report is created in `multiqc_report.html`

## Trimming
[Trim Galore](https://github.com/FelixKrueger/TrimGalore) is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt) and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to consistently apply adapter and quality trimming to FastQ files.
```
module load intel python java trimgalore

for bam in fastq/*.fastqc.gz; do

nameHeader=$(basename ${fastqc/_R1.fastq.gz/})

echo " 

trim_galore --illumina --paired --trim1 $nameHeader ${nameHeader/_R1/_R2} -o trimmed/$(basename ${nameHeader/_R1.fastq.gz/})

" >> /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;
sbatch /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;
sed -i '10,$d' /slgpfs/projects/idib57/GA_002_21_Radon/sbatch.12h.cpt1 ;
done
```

## Normalization

[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) is a package for the analysis of gene expression data arising from microarray or RNA-seq
technologies, and [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) implements a range of statistical methodology based on the negative binomial distributions, including empirical Bayes estimation, exact tests, generalized linear models and quasi-likelihood tests. With this packages we can normalize our data
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

![no normalized](https://github.com/ABedmar/radon-project/blob/main/images/no_norm.png)
<p align="center">
    Not normalized samples
</p>

<br/><br/>

![normalized](https://github.com/ABedmar/radon-project/blob/main/images/norm.png)
<p align="center">
    Normalized samples
</p>

<br/><br/>

# Analysis

## Deconvolution


## Gene expression analyis
[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) is a computational method that determines whether an a priori defined set of genes shows statistically
significant, concordant differences between two biological states. To compare enriched gene sets acording to radon exposure we create the file containing phenotypes `phenotype_labels_ALL_c.cls`.

## TMB
Tumor mutational burden is the total number of mutations (changes) found in the DNA of cancer cells.Tumor mutational burden is being used as a type of biomarker and in this analysis we use the following tool to obtain it.
[maftools](https://github.com/PoisonAlien/maftools) provides a comprehensive set of functions for processing [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) files and to perform most commonly used analyses in cancer genomics. 

The function `tcgaCompare` uses mutation load from TCGA [MC3](https://gdc.cancer.gov/about-data/publications/mc3-2017) for comparing muttaion burden against 33 TCGA cohorts. 
In this particular case we are only interested in the TMB of our samples:

```
#module load gcc/8.1.0 pcre2/10.35 R/4.0.3

library(maftools)
files=dir(pattern=".maf$")
maf=merge_mafs(files)

laml.mutload = tcgaCompare(maf = maf, cohortName = 'Radon-LAML', logscale = TRUE, capture_size = 50)

rb=(laml.mutload$mutation_burden_perSample)
rb= rb[rb$cohort == "Radon-LAML",]
```

`rb` is a table containing all samples and their TMB.
