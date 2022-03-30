# radon-project
Retrospective analyis of radon induced lung cancer in rats. This analysis has the objective of describeing the tumor molecular profile related to the immune profile of radon-induced lung cancer vs. spontaneous lung cancer in rats and to characterize pathologically and transcriptomically all the TME in lung primary tumors with a particular interest in CD8+ T cells (inflamed-TME) in rats exposed or not to radon.

# Table of contents
- [Pipeline](https://github.com/ABedmar/radon-project/blob/main/README.md#pipeline)
    - [Quality control](https://github.com/ABedmar/radon-project/blob/main/README.md#Quality_control)
    - [Trimming](https://github.com/ABedmar/radon-project/blob/main/README.md#Trimming)
    - [BAM files](https://github.com/ABedmar/radon-project/blob/main/README.md#bam-files)
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
## BAM files


## Normalization

[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) is a package for the analysis of gene expression data arising from microarray or RNA-seq
technologies, and [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) implements a range of statistical methodology based on the negative binomial distributions, including empirical Bayes estimation, exact tests, generalized linear models and quasi-likelihood tests. With this packages we can normalize our data.
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


<p align="center">
    <img src="https://github.com/ABedmar/radon-project/blob/main/images/no_norm.png" width="900" height="350">
</p>
<p align="center">
    Not normalized samples
</p>

<br/><br/>

<p align="center">
    <img src="https://github.com/ABedmar/radon-project/blob/main/images/norm.png" width="900" height="350">
</p>
<p align="center">
    Normalized samples
</p>

<br/><br/>

# Analysis

## Deconvolution
Deconvolution refers to separating a heterogeneous mixture signal into its components. We are interested in deconvolution because we want to define the pathological and transcriptomic immune profile of the TME of radon-induced vs. spontaneous lung cancer in rats, to identify the specific immune patterns associated with radon. [EcoTyper](https://www.sciencedirect.com/science/article/abs/pii/S0092867421010618?via%3Dihub) is a framework for systematically identifying cell states and cellular communities (ecotypes) from gene expression data.

## Gene expression analyis
[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) is a computational method that determines whether an a priori defined set of genes shows statistically
significant, concordant differences between two biological states. To compare enriched gene sets acording to radon exposure we create the file containing phenotypes `phenotype_labels_ALL_c.cls`.
GSEA Desktop is a free genomic analysis program written in the Java(tm) language implementing the GSEA method while providing preprocessing tools along with further analysis methods and visualizations.
Specifically, for this analysis we used the following parameters:
- Gene set database: `c2.all.v7.5.1.symbols.gmt`
- Number of permutations: 1000
- Collapse to gene symbols
- Chip platform: `Rat_Ensembl_Transcript_ID_Human_Orthologs_MSigDB.v7.5.1.chip`

## TMB
Tumor mutational burden is the total number of mutations (changes) found in the DNA of cancer cells.Tumor mutational burden is being used as a type of biomarker and in this analysis we use the following tool to obtain it.
[maftools](https://github.com/PoisonAlien/maftools) provides a comprehensive set of functions for processing [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) files and to perform most commonly used analyses in cancer genomics. 

To get the MAF files we convert [VCF](http://samtools.github.io/hts-specs/) files with [vcf2maf](https://github.com/mskcc/vcf2maf).

```
module load htslib/1.10.2 intel/2018.3 bcftools/1.10.2  perl/5.26.2 gcc/8.1.0 vep/97.3 samtools/1.9 vcf2maf

genome=/slgpfs/projects/idib57/data/genome/Rattus/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa

for vcf in *.vcf; do
echo "

vcf2maf --input-vcf $vcf --output-maf ${vcf/.vcf/.maf} --tumor-id ${vcf/.vcf/} --ref-fasta $genome --vep-overwrite

" >> /home/idib57/idib57798/scripts/sbatch/sbatch.1h.cpt1 ;
sbatch < /home/idib57/idib57798/scripts/sbatch/sbatch.1h.cpt1 ;
sed -i '10,$d' /home/idib57/idib57798/scripts/sbatch/sbatch.1h.cpt1

done
```

<br/><br/>

Now we merge all the MAF files using the function `merge_mafs`.

```
setwd("/slgpfs/projects/idib57/giancarlo/radon/WES/vep_merged_filtered_nopolimorph")
library(maftools)

files=dir(pattern=".maf$")maf=merge_mafs(files)
maf=merge_mafs(files)
```

<br/><br/>

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
