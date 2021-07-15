# RNAseq pipelines - current status

Created: 2020/09/24 15:19:42
Last modified: 2021/07/15 13:20:13

- **Aim:** Evaluate the current pipelines available for processing **RNA-seq** data. This will help us decide if there is a pipeline currently available for our use, one we could adapt, or if we will need to create a RNA-seq pipeline from scratch

## Table of contents

- [RNAseq pipelines - current status](#rnaseq-pipelines---current-status)
  - [Table of contents](#table-of-contents)
  - [Overview](#overview)
  - [Pipelines](#pipelines)
    - [VIPER: Visualization Pipeline for RNA-seq](#viper-visualization-pipeline-for-rna-seq)
    - [TRAPLINE](#trapline)
    - [HppRNA](#hpprna)
    - [DRAGEN RNA Pipeline](#dragen-rna-pipeline)
    - [rna-seq-star-deseq2](#rna-seq-star-deseq2)
    - [rna-seq-kallisto-sleuth](#rna-seq-kallisto-sleuth)
    - [single-cell-rna-seq](#single-cell-rna-seq)
    - [single-cell-drop-seq](#single-cell-drop-seq)
    - [National Cancer Institute mRNA quantification analysis pipeline](#national-cancer-institute-mrna-quantification-analysis-pipeline)
    - [CloseCall](#closecall)
    - [snDrop_prep](#sndrop_prep)
    - [RNAsik-pipe](#rnasik-pipe)
    - [exceRpt](#excerpt)
    - [GEO2RNAseq](#geo2rnaseq)
    - [GeneTEFlow](#geneteflow)
    - [Parabricks](#parabricks)
    - [DiMSum](#dimsum)
    - [scTyper](#sctyper)
    - [RiboMiner](#ribominer)
    - [SLFinder](#slfinder)
    - [Cell Ranger pipeline (10×Genomics)](#cell-ranger-pipeline-10genomics)
    - [Seurat](#seurat)
    - [irap](#irap)
    - [ymp](#ymp)
    - [rnaseq-pipeline](#rnaseq-pipeline)
    - [umitools](#umitools)
    - [RNA-seq](#rna-seq)
    - [snakemake_RNA-seq](#snakemake_rna-seq)
    - [snakemake_fastqc](#snakemake_fastqc)
    - [snakemake_deeptools](#snakemake_deeptools)
    - [omics_pipe](#omics_pipe)

## Overview

A summary of some of the RNA-seq data analysis tools by [Ruairi J MacKenzie, 2018](https://www.technologynetworks.com/genomics/articles/rna-seq-basics-applications-and-protocol-299461):

"Tools like Sailfish, RSEM and BitSeq12 will help you quantify your expression levels, whilst tools like MISO, which quantifies alternatively spliced genes, are available for more specialized analysis. There is a library of these tools out there, and reading reviews and roundups are your best way to find the right tool for your research."

- There also exists an R framework for comparing pipelines (see [here](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02136-7))
- [This paper](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02136-70) outlines a number of papers that compare single-cell RNA pipelines
- A [recent paper](https://www.embopress.org/doi/full/10.15252/msb.20188746) (2019) that describes the best practices in single-cell RNA-seq analysis
- A [recent paper](AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAqUwggKhBgkqhkiG9w0BBwagggKSMIICjgIBADCCAocGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMdbpTwGltlCc7WTjbAgEQgIICWEvhDFT) (2018) that describes a website that acts as a reference set/database of human long non-coding RNAs
- A [recent paper](https://www.nature.com/articles/s12276-018-0071-8.pdf) (2018) that describes Single-cell RNA sequencing technologies and bioinformatics pipelines

##  Pipelines

### VIPER: Visualization Pipeline for RNA-seq

- Paper [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2139-9) and github repo [here](https://github.com/hanfeisun/viper-rnaseq)
- Open source
- Workflow language - snakemake
- Package management with Conda
- Hasn't been updated in four years - may not be actively maintained/supported
- [STAR](https://github.com/alexdobin/STAR) for alignment
- [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks) for transcript assembly, normalisation of read counts etc.
- bams in BigWig format to allow visualisation in genome browser
- [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) for fusion gene discovery for paired-end data
- Read quality metrics
- Differential expression and pathway analysis - can use [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) or build in a preferred differential expression method
- Immunology analysis module available
- Whole-genome SNV calling (human and mouse) module available
- Viral analysis (human samples only) module available
- Batch effect correction module available
- This pipeline focuses on data visualisation, using Snakemake, parallelization/speed, ease of use (such as less options, easier to install)
- They say they differ from two other pipelines available by "the number of features included, package management software, and reporting functionalities"
- They say they focus on the best practice software rather than including many options (like HppRNA)
- See their comparison to other pipelines [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2139-9/tables/1)

### TRAPLINE

- Package management with Galaxy

### HppRNA

- github repo [here](https://github.com/NextGenBioinformatics/hppRNA)
- RNA-Seq analysis of numerous samples
- They describe as "parameter-free"
- Hasn't been updated in two years - may not be actively maintained/supported

### DRAGEN RNA Pipeline

- Website [here](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)
- Not open source?
- For research use only

### rna-seq-star-deseq2

- github repo [here](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
- Open source
- Workflow language - snakemake
- Package management with Conda
- Differential expression analysis with [STAR](https://github.com/alexdobin/STAR) and [Deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

### rna-seq-kallisto-sleuth

- github repo [here](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth)
- Open source
- Workflow language - snakemake
- Package management with Conda
- Differential expression analysis with [Kallisto](https://pachterlab.github.io/kallisto/) and [Sleuth](https://pachterlab.github.io/sleuth/)

### single-cell-rna-seq

- github repo [here](https://github.com/snakemake-workflows/single-cell-rna-seq)
- Open source
- Workflow language - snakemake
- Package management with Conda
- Single cell RNA-seq workflow

### single-cell-drop-seq

- github repo [here](https://github.com/snakemake-workflows/single-cell-drop-seq)
- Open source
- Workflow language - snakemake
- Package management with Conda
- Single cell RNA-seq workflow
- [STAR](https://github.com/alexdobin/STAR) for alignment

### [National Cancer Institute](https://www.cancer.gov/) mRNA quantification analysis pipeline

- Website [here](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-analysis-pipeline)
- [STAR](https://github.com/alexdobin/STAR) for alignment
- GDC gene fusion pipeline using [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
- Arriba gene fusion pipeline

### CloseCall

- github repo [here](https://github.com/StevenWingett/CloseCall) and paper [here](https://www.nature.com/articles/s41597-020-0372-3)
- Pipeline for processing RNA-RNA proximity data
- Open source
- Mapping and QC
- Monte Carlo Simulation to identify statistically significant RNA-RNA interactions
- Actively maintained (last commit 13 months ago)
- Written in perl, java and R (not in a workflow language?)

### snDrop_prep

- github repo [here](https://github.com/chensong611/snDrop_prep) and paper [here](https://www.nature.com/articles/s41467-019-10861-2)
- Open source
- Single-nucleus RNA-sequencing pipeline
- Python and bash (not in a workflow language)

### RNAsik-pipe

- github [here](https://github.com/MonashBioinformaticsPlatform/RNAsik-pipe) and website [here](https://monashbioinformaticsplatform.github.io/RNAsik-pipe/)
- Open source
- Written in BigDataScript (bds) (with underlying Java) (not in a workflow language)

### exceRpt

- github [here](https://github.com/gersteinlab/exceRpt), website [here](https://rna-seqblog.com/excerpt-a-comprehensive-analytic-platform-for-extracellular-rna-profiling/)
- Open source
- Extracellular RNA profiling
- Hasn't been updated in two years - may not be actively maintained/supported
- Written in java, R, bash, pl (not in a workflow language)

### GEO2RNAseq

- Hosted on Anaconda [here](https://anaconda.org/xentrics/r-geo2rnaseq), paper [here](https://www.biorxiv.org/content/10.1101/771063v1.full)
- Open source

### GeneTEFlow

- github [here](https://github.com/zhongw2/GeneTEFlow), paper [here](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0232994&type=printable)
- Differential expression analysis of genes and locus-specific transposable elements from RNA sequencing
- Open source
- Nextflow
- [Trimmomatic](https://github.com/timflutre/trimmomatic) for fastq trimming
- Actively maintained (last commit 3 months ago)

### Parabricks

- Some RNA tools available in the GPU-accelerated toolkit
- Not open source
- [Star](https://www.nvidia.com/en-us/docs/parabricks/star/) for alignment
- [Star-Fusion](https://www.nvidia.com/en-us/docs/parabricks/star-fusion/) for fusion gene discovery for paired-end data

### DiMSum

- github [here](https://github.com/lehner-lab/DiMSum), paper [here](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02091-3)
- Deep mutational scanning (DMS) enabling the multiplexed measurement of the effects of thousands of variants of proteins, RNAs, and regulatory elements
- Open source
- Written in R (not in a workflow language)

### scTyper

- github [here](https://github.com/omicsCore/scTyper), paper [here](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03700-5)
- Cell typing analysis of single-cell RNA-seq data
- Written in R (not in a workflow language)

### RiboMiner

- github [here](https://github.com/xryanglab/RiboMiner), paper [here](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03670-8)
- Open source
- Mining multi-dimensional features of the translatome with ribosome profiling data
- Four function parts: Quality Control, Metagene Analysis, Feature Analysis, Enrichment Analysis:
- Written in python (not a pipeline, not in a workflow language)

### SLFinder

- github [here](https://github.com/LBC-Iriarte/SLFinder), paper [here](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03610-6)
- Open source
- Novel identification of splice-leader sequences
- Actively maintained (last commit 3 months ago)
- Written in bash (not in a workflow language)

### Cell Ranger pipeline (10×Genomics)
- website [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
- Looks like it's open source
- A set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis
- t-SNE is implemented (according to [this paper](https://www.nature.com/articles/s12276-018-0071-8.pdf))

### Seurat
- website [here](https://satijalab.org/seurat/)
- QC, analysis, and exploration of single-cell RNA-seq data
- Written in R (not in a workflow language)
- t-SNE is implemented (according to [this paper](https://www.nature.com/articles/s12276-018-0071-8.pdf))

### irap

- github [here](https://github.com/nunofonseca/irap), paper [here](https://academic.oup.com/bioinformatics/article/31/5/665/2748143)
- Flexible RNA-seq analysis pipeline that allows the user to select and apply their preferred combination of existing tools for mapping reads, quantifying expression and testing for differential expression
- Open source
- Written in R, perl, bash (not in a workflow language)
- Hasn't been updated in two years - may not be actively maintained/supported

### ymp

- github [here](https://github.com/epruesse/ymp)
- Flexible omics pipeline (QC, trimming, contaminant removal), assemble metagenomes, annotate assemblies, or assemble and quantify RNA-Seq transcripts, offering a choice of tools for each of those processing stages
- Written in python (not in a workflow language)
- Open source
- Actively maintained (last commit a month ago)

### rnaseq-pipeline

- github [here](https://github.com/PavlidisLab/rnaseq-pipeline)
- Written in python, R (not in a workflow language)
- Open source
- Actively maintained (last commit 3 months ago)

### umitools

- github [here](https://github.com/weng-lab/umitools) and paper [here](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4933-1)
- A toolset for handling sequencing data with unique molecular identifiers (UMIs)
- Unique molecular identifiers (UMIs) are a type of molecular barcoding that provides error correction and increased accuracy during sequencing. These molecular barcodes are short sequences used to uniquely tag each molecule in a sample library
- Open source
- Written in python (not in a workflow language)
- Hasn't been updated in two years - may not be actively maintained/supported

### RNA-seq

- github [here](https://github.com/biowdl/RNA-seq)
- A Biowdl workflows usable for processing RNA-seq data. This pipeline will performs QC (including adapter clipping), mapping, variant-calling and expression quantification
- Open source
- Workflow Description Language
- Actively maintained (last commit 20 days ago)

### snakemake_RNA-seq

- github [here](https://github.com/WilliamJeong2/snakemake_RNA-seq)
- A snakemake pipeline for the analysis of RNA-seq data that makes use of hisat2 and Stringtie
- Snakemake
- Open source
- Actively maintained (last commit 25 days ago)

### snakemake_fastqc

- github [here](https://github.com/AngryMaciek/snakemake_fastqc)
- A small snakemake workflow for FastQC
- Snakemake
- Open source
- Actively maintained (last commit 25 days ago)

### snakemake_deeptools

- github [here](https://github.com/AngryMaciek/snakemake_deeptools)
- Snakemake pipeline for RNA-Seq data analysis with deepTools
- Snakemake
- Open source
- Actively maintained (last commit 25 days ago)

### omics_pipe

- pypi [here](https://pypi.org/project/omics_pipe/), paper [here](https://watermark.silverchair.com/btv061.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAsQwggLABgkqhkiG9w0BBwagggKxMIICrQIBADCCAqYGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMunJwtqi7vVQY1pnxAgEQgIICdwI8KS4jfan_7kqncenMG6a7TWMf1z3nUUF2P-z6JQjhHKWVHh33yTeb3hafkRzYW8VTwXPY3lNV7Rew8P_Dq5HLIjhyaq1lhGb6cRuOtcAZyCAda15jz5-vfRLXP29YIjZ4CYdoeqPcTgXVaa9pX54-wD6tLsSLOatfbCLPVT94obKI2soNh9FSN_UvAk3-gWUaCs0CGt9eJ6bxf2F9hh3o0FgbSAt4_gdaAiQoVaB3CF1MjUDf3gL-QnNd-x2_q4EqqsZ350idWcJMK0jTKsEuUbFnWuUZJw03gr-D8gL_MgtiQ71GXV45BGBJ9JEDFoKDPfTNIu5bB1MOdl76VBTU1DcCJNv1P_YSxJ0dJp_AtNGn4RhT2KnnQwDfbvxbjGcIPSD32pUwYrbZ_QgCwK1icMN_5JsBiMN2BzP9j2cHrMZucRCRmye77crMXRyFqdeTaocacnO0KZoFLg_2r_N5JyZ4v8hTPEFsT2IynksHoxbnL6CMpvfQFZP5wFdSjKYui4luoS5XgmrxzcVUdbah4SYEOPqC2ZInkidtBH-Gew5p71FAfJAfc12aM6iVpdrPLJATHG9h60CbmixkjBCdu8pO9lpaHX8wJpwABDymWNrZ_6r1zaLGGabYP3qQszMk9Jq_c_mQVZhLk76SKkL9iav5X-Yo3T6HFLt6fQWcuEOnbpDxHShDkJMw6L8mr-xjlP7F7dH03A7n7lYXnWfOarjtFyKCSh6PP6g-KvXfRbqfvwk109owgJBn16RsKUNEFBe-teN6jx-jLgbjahyfumE4klbF-Vbot8xLww-ZLjV5Dk3wFsQu7MqwIo5agToy15oq9XE)
- A computational framework that automates multi-omics data analysis pipelines on high performance compute clusters and in the cloud. It supports best practice published pipelines for RNA-seq, miRNA-seq, Exome-seq, Whole-Genome sequencing, ChIP-seq analyses and automatic processing of data from The Cancer Genome Atlas (TCGA)
- Older - paper published in 2015
