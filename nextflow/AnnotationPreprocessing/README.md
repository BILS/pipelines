# Annotation preprocessing pipeline

## Quickstart

```
nextflow run -profile nbis AnnotationPreprocessing.nf --genome_assembly '/path/to/genome_assembly.fasta'
```

## Parameters

* `genome_assembly`: The path to the genome assembly in quotes.
* `outdir`: The name of the results folder
* `min_length`: The minimum_length for fasta sequences in the assembly to be. 

## Stages

* Filter: Remove fasta sequences less than `min_length` bases
* Get Assembly Metrics: Calculate and plot summary metrics on the assembly

