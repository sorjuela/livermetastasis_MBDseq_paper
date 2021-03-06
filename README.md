### Sequencing of DNA captured with a methyl-CpG binding domain reveals a similar hypermethylator phenotype in colon cancers and liver metastases.

## Code and data analysis

This repository contains the code for the following paper:

* S Orjuela, M Menigatti, P Schraml, P Kambakamba, MD Robinson and G Marra: [The DNA hypermethylation phenotype of colorectal cancer liver metastases resembles that of the primary colorectal cancers](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-06777-6).

The structure of this repository is as follows:

* `Data/`
Metadata from samples analyzed and outputs from `csaw` analysis (`.RData` objects with table of identified DMRs).
* `R/`
R scripts with code for DMR detection analysis, and code for plotting most of figures in paper.
* `bash/`
shell scripts for pre-processing of raw data (FASTQ files).