# Single cell co-expression QTL analysis

This repository contains the code to generate the results and figures from

**Identification of genetic variants that impact gene co-expression relationships using large-scale single-cell data**

Shuang Li *, Katharina T. Schmid *, Dylan de Vries *\*, Maryna Korshevniuk *\*, Corinna Losert, Roy Oelen, Irene van Blokland, BIOS Consortium, sc-eQTLgen Consortium, Hilde E. Groot, Morris A. Swertz, Pim van der Harst, Harm-Jan Westra, Monique van der Wijst, Matthias Heinig †, Lude Franke †

\* These authors contributed equally<br>
** These authors contributed equally<br>
† These authors contributed equally

Preprint: https://www.biorxiv.org/content/10.1101/2022.04.20.488925v1

## Overview

The code for the analysis is separated in different steps, each in its own subdirectory including a README file:
* Exploring different different association metrics and other GRN construction tools, including pseudotemporal based ones and combing cells to meta cells [01_association_metrics/](01_association_metrics/)
* Comparing correlation between the different single cell data sets, with bulk data and CRISPR knock-out data; testing potential occurence of Simpson's paradox [02_correlation_evaluation/](02_correlation_evaluation/)
* Comparing correlation between different cell types and between different individuals within one cell type [03_celltype_individual_comparison/](03_celltype_individual_comparison/)
* Running eQTL and coeQTL mapping pipelines followed by replication in bulk and technical evaluation of the co-eQTLs, such as correlation distribution, sub cell type effects and effects of subsampling cells or donors [04_coeqtl_mapping/](04_coeqtl_mapping/)
* Interpretation of co-eQTL results based on GWAS annotation and different enrichment analyses (GO enrichment, TFBS enrichment using Remap database and GWAS enrichment using MAGMA) [05_coeqtl_interpretation/](05_coeqtl_interpretation/)

## Software requirements

Most code is implemented in R or python, the required packages are documented in the respective yaml files (`conda_env_R.yml` and `conda_env_python.yml`) and can be used to setup the respective [conda environments](https://docs.conda.io/en/latest/):

```
conda env create -f conda_env_R.yml
conda activate r_env
```

The following R packages are not part of conda and need to be added afterwards (if the respective code part should be run):

```
if (!require("BiocManager")) install.packages('BiocManager') 
BiocManager::install("tanaylab/metacell")

devtools::install_github("heiniglab/scPower")
```

Further external tools were used:

* To calculate the eQTLs See the documentation [here](https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook-for-RNA-seq-data).
* To calculate the co-eQTLs: [mbQTL](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL).
* To perform GWAS enrichment analysis on the co-eQTL genes: [MAGMA](https://ctg.cncr.nl/software/magma).





