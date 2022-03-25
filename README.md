# radon-project
Retrospective analyis of radon induced lung cancer rats



## Quality control
With MultiQC...

## Trimming
Using trim galore...


## Deconvolution


## Gene expression analyis


## TMB

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
