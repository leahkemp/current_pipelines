# RNAseq pipelines - current status

Created: 2020/09/24 15:19:42
Last modified: 2021/08/03 16:28:54

- **Aim:** Evaluate the current pipelines available for processing **RNA-seq** data

## Table of contents

- [RNAseq pipelines - current status](#rnaseq-pipelines---current-status)
  - [Table of contents](#table-of-contents)
  - [Pipelines](#pipelines)
  - [Notes](#notes)

## Pipelines

See [RNAseq pipelines](./RNAseq_pipelines.csv) for a spreadsheet of pipelines

## Notes

My current favorite pipelines:

- [nf-core/rnaseq](https://github.com/nf-core/rnaseq) because
  - It has many stars/forks/contributors - lots of people using and developing it
  - Open source, workflow language, package management, cloud computer support, resource management, good documentation
  - Good analysis - several QC steps, trimming, removal of ribosomal RNA, choice of multiple aligners/quantification tools (STAR, salmon, RSEM, HiSAT2, UMI-based deduplication, mark duplicates, coverage files)

- [rna_fq2bam](https://docs.nvidia.com/clara/parabricks/v3.5/text/rna.html#rna-fq2bam)
  - GPU accelerated with parabricks
  - Uses STAR, gatk SortSam, gatk MarkDuplicates under the hood
  - Looks like it would be fairly equivalent to part of [nf-core/rnaseq](https://github.com/nf-core/rnaseq) (which also uses STAR, picard MarkDuplicates (which is the basis of gatk MarkDuplicates). However, gatk SortSam would be slightly different to the samtools sort command the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline uses, but I can't imagine different sorting tools produce very different outputs since they all just sort the data
  - We could develop a subworkflow in [nf-core/rnaseq](https://github.com/nf-core/rnaseq) to run (part of the pipeline) in a GPU accelerated mode for those that have a parabricks license

- [star_fusion](https://docs.nvidia.com/clara/parabricks/v3.5/text/rna.html#star-fusion)
  - GPU accelerated with parabricks
  - Uses STAR-Fusion under the hood
  - A standalone STAR that is used in [rna_fq2bam](https://docs.nvidia.com/clara/parabricks/v3.5/text/rna.html#rna-fq2bam)

- [rna_pipeline](https://docs.nvidia.com/clara/parabricks/v3.5/text/rna_pipeline.html)
  - GPU accelerated with parabricks
  - "GATK best practice"
  - This is a variant discovery pipeline, not an transcript count pipeline such as [nf-core/rnaseq](https://github.com/nf-core/rnaseq)
  - This pipeline could be used alongside [nf-core/rnaseq](https://github.com/nf-core/rnaseq) if we're interested in variant discovery

- [snakemake-workflows/rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
  - It has many stars/forks/contributors - lots of people using and developing it
  - Open source, workflow language, package management, cloud computer support, resource management, good documentation
  - Good analysis - several QC steps (custom python scripts and multiqc), automatic download of ensemble sequences/annotation files, SRA fastq files, trimming with [cutadapt](https://cutadapt.readthedocs.io/en/stable/)), trimming, alignment with STAR and differential expression
  - Pro compared to [nf-core/rnaseq](https://github.com/nf-core/rnaseq) - does a differential expression/pca analysis (with the [deseq2 package](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) and the [affycoretool plotPCA function](https://www.rdocumentation.org/packages/affycoretools/versions/1.44.2/topics/plotPCA)) - I don't see this as a major benefit since we have all the code to run these analyses, and it would be much easier to customise our plots as well as run a differential expression analysis with additional tools (such as [limma/voom](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)) with our own code
  - Con compared to [nf-core/rnaseq](https://github.com/nf-core/rnaseq) - no multiple option for aligners/quantification tools, doesn't look like it does the UMI-based deduplication

- [epruesse/ymp](https://github.com/epruesse/ymp)
  - Not as popular (at the moment) as some of the other major pipelines but still has several people developing it
  - This is a flexible omics pipeline that has a focus on allowing the flexibility for the user to choose (and explore) the tools they want to analyse their data with
  - Based on Snakemake
  - This could be a really good way to go if it works well
