# Augustus training pipeline

## Quickstart

```
nextflow run -profile nbis,docker AugustusTraining.nf \
  --genome '/path/to/genome_assembly.fasta' \
  --gff_annotation 'path/to/annotation.gff3'
```

## Parameters

* `genome`: The path to the genome assembly in quotes.
* `gff_annotation`: The path to the gff annotation in quotes.
* `outdir`: The name of the results folder

* `gff_gene_model_filter_options`: Options to be passed to the filter by gene model script (default:'-c -r -d 500 -a 0.3').

* `codon_table`: The number of the codon table to use for translation.

* `test_size`: The size of the test data set
* `flank_region_size`: The fize of the flank region to include. 

## Stages


