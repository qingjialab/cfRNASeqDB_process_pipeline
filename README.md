# cfRNASeqDB: a cell-free RNA-Seq database of gene expression and microbiome analysis for liquid biopsy studies
---
## Introduction
cfRNASeqDB refers to gene expression and microbiome analysis of cell-free RNA-Seq database. It is a comprehensive database contained 2240 paired-end sequencing and 145 single-end sequencing datasets respectively, including 15 types of disease conditions and 3 types of health conditions from 10 open cfRNA-Seq database. cfRNASeqDB hosts all the data and provides a user-friendly website for users to freely explore, query and download.

## Run snakemake locally
RNA-seq raw data were downloaded from the NCBI Gene Expression Omnibus (GEO, [http://www.ncbi.nlm.nih.gov/geo/](http://www.ncbi.nlm.nih.gov/geo/) and Short Read Archive (SRA, http://www.ncbi.nlm.nih.gov/sra/), using the Aspera high-speed file transfer protocol ([http://asperasoft.com/)](http://asperasoft.com/)). The paired-end sequencing and single-end sequencing datasets were processed using paired_end_snakemake_pipeline and single_end_snakemake_pipeline separately. Submit the following command at the terminal to run the pipeline.

```
$ snakemake --use-conda -c all -s snakefile --max-jobs-per-second 5 --keep-going
```

## Usage
This bioinformatics pipeline processes RNA-seq raw data to extract gene expression, bacteria and virus abundance information simultaneously.  Users can adopt the pipeline to process other cfRNA-seq or RNA-Seq datasets to extract useful information.
