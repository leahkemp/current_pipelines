# small RNAseq pipelines - current status

Created: 2020/09/24 15:19:42
Last modified: 2021/07/15 13:20:14

- **Aim:** Evaluate the current pipelines available for processing **small RNA-seq** data. This will help us decide if there is a pipeline currently available for our use, one we could adapt, or if we will need to create a small RNA-seq pipeline from scratch

## Table of contents

## Sipelines

### smrnaseq

- github [here](https://github.com/nf-core/smrnaseq) and paper [here](https://www.biorxiv.org/content/10.1101/610741v1)
- Best-practice analysis pipeline used for small RNA sequencing data
- Open source
- Nextflow
- Deployable to SLURM and AWS
- Actively maintained (last commit 12 months ago)

### sports1.1

- github [here](https://github.com/junchaoshi/sports1.1) and paper [here](https://www.sciencedirect.com/science/article/pii/S1672022918300445)
- Small non-coding RNA annotation Pipeline Optimized for rRNA- and tRNA-Derived Small RNAs
- Open source
- Written in R and perl (not in a workflow language)
- Actively maintained (last commit 2 months ago)

### smallseq

- github [here](https://github.com/eyay/smallseq)
- Analyze small RNAs from single-cells
- Open source
- Written in python (not in a workflow language)
- Hasn't been updated in three years - may not be actively maintained/supported

### SnapT

- github [here](https://github.com/ursky/SnapT)
- A Small non-coding RNA annotation pipeline for Transcriptomic or metatranscriptomic data
- Written in python (not in a workflow language)
- Open source
- Actively maintained (last commit 6 months ago)

### short-ncrna-annotation

- github [here](https://github.com/SimonSchafferer/short-ncrna-annotation)
- Annotation of commonly used interval/range based data such as the browser extended display format (BED), or the general feature format (GFF). In addition, it provides sequence based annotation, by employing the NCBI blast+ software
- Written in R (not in a workflow language)
- Open source
- Hasn't been updated in six years - may not be actively maintained/supported

### ncRNA_Pipeline

- github [here](https://github.com/navygit/ncRNA_Pipeline)
- Open source
- Written in perl (not in a workflow language)
- Hasn't been updated in five years - may not be actively maintained/supported

### FlaiMapper

- github [here](https://github.com/yhoogstrate/flaimapper), paper [here](https://academic.oup.com/bioinformatics/article/31/5/665/2748143)
- Annotation of small ncRNA-derived fragments using RNA-seq high-throughput data
- Open source
- Written in bash (not in a workflow language)
- Hasn't been updated in four years - may not be actively maintained/supported

### exceRNApipeline

- github [here](https://github.com/zhuchcn/exceRNApipeline)
- Data processing pipeline for extracellular small RNA-seq from human specimen.The pipeline is designated for running on HPC with the job management system SLURM.
- Preprocess: remove adapters and trimming low quality nucleotides, using HTStream.
- UniVec: map to the NCBI's UniVec database to remove vector origins
- RiboRNA: map to the rRNA sequences
- Human Genome: map to human genome
- Repetitive Elements: map to RepeatMasker's repetitive elements sequences
- SILVA: map to SILVA's ribosomal rRNA gene of bacteria, archaea, and fungi
- Bacteria: map to all bacteria genomes in ensemble's
- Only supports single end sequencing data
- Snakemake
- Open source
- Actively maintained (last commit 6 months ago)

### sRNAflow

- github [here](https://github.com/zajakin/sRNAflow)
- Analysis of small RNA that fulfills the specific needs for samples derived from biofluids
- Nextflow
- Open source
- Actively maintained (last commit 28 days ago)

### gorap

- github [here](https://github.com/koriege/gorap), website [here](https://www.rna.uni-jena.de/research/software/)
- Screens genomic sequences for all non-coding RNAs present in the Rfam database
- Provides ncRNA based reconstruction of phylogenetic trees and is able to perform de novo predictions including TPM calculations from RNA-Seq experiments
- Open source
- Written in perl (not in a workflow language)
- Actively maintained (last commit 7 months ago)

## Other things of possible interest

- Prediction and annotation of tRNA-derived small RNAs [here](https://github.com/wangqinhu/tsRFinder)
- A JBrowse plugin to support small RNA visualization [here](https://github.com/bhofmei/jbplugin-smallrna)
- A set of tools related to the Rfam (database of non-coding RNA families https://rfam.org/) production pipeline [here](https://github.com/Rfam/rfam-production)
- Script to map annotated ncRNAs (and other elements) to a chromosome and visualise [here](https://github.com/fanagislab/draw_annotation/tree/master/bin)
- [HPC-REDItools - A tool for large-scale RNA-editing analysis](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03562-x) (RNA-editing is a molecular process through which some cells can make discrete changes to specific nucleotide sequences within an RNA molecule after it has been generated by RNA polymerase)

## Notes

-	[smrnaseq](https://github.com/nf-core/smrnaseq) is still the most robust pipeline I’ve come across in terms of:
    - being in a workflow language
    - being specific to both small and non-coding RNA
    - being deployable to cluster/cloud environments
-	There are a plethora of pipelines available for RNA-seq (more than I’ve covered in this document) 
- There seems to very few pipelines dedicated to both small and non-coding RNA-seq
-	There are a bunch of interesting scripts available that might end up being useful to incorporate or compliment whichever pipeline we end up using, or be used for future analyses:
    - [CloseCall](https://github.com/StevenWingett/CloseCall) can identify RNA-RNA interactions 
    - [DiMSum](https://github.com/lehner-lab/DiMSum) uses deep mutational scanning to enable the multiplexed measurement of the effects of thousands of variants of proteins, RNAs, and regulatory elements
    - [HPC-REDItools](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03562-x) is a tool for large-scale RNA-editing analysis
    - [draw_annotation](https://github.com/fanagislab/draw_annotation) to map annotated ncRNAs (and other elements) to a chromosome and visualise
-	There are some other interesting pipelines: 
    - [sRNAflow](https://github.com/zajakin/sRNAflow) is a Nextflow pipeline to analyse small RNA that is tailored to meet the specific needs of samples derived from biofluids
    - [exceRNApipeline](https://github.com/zhuchcn/exceRNApipeline) is a Snakemake pipeline for extracellular small RNA-seq
-	There are a few interesting annotation pipelines:
    - [FlaiMapper](https://github.com/yhoogstrate/flaimapper )
    - [short-ncrna-annotation](https://github.com/SimonSchafferer/short-ncrna-annotation) 
    - However, these are typically older repos, and the scripts aren’t wrapped in a workflow language
    - Additionally it looks like there is some overlap in the analyses in smrnaseq and FlaiMapper 
-	Some interesting and more general RNA-Seq pipelines: 
    - [viper-rnaseq](https://github.com/hanfeisun/viper-rnaseq)
    - [rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
    - [rna-seq-kallisto-sleuth](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth)
    - [RNA-seq](https://github.com/biowdl/RNA-seq)
    - [snakemake_deeptools](https://github.com/AngryMaciek/snakemake_deeptools)
    - [omics_pipe](https://pypi.org/project/omics_pipe/)
