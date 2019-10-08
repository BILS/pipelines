# Functional annotation preparation pipeline

## Quickstart

```
nextflow run -profile nbis,docker FunctionalAnnotationPreparation.nf \
  --genome '/path/to/genome_assembly.fasta' \
  --gff_annotation 'path/to/annotation.gff3'
```

## Parameters

* `genome`: The path to the genome assembly in quotes.
* `gff_annotation`: The path to the gff annotation in quotes.
* `outdir`: The name of the results folder

* `records_per_file`: The number of records per file. A parallelisation option 
to improve blast and interproscan speed.

* `blast_db`: The path to protein database files in quotes.

## Stages


