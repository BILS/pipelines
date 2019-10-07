# Transcript assembly pipeline

## Quickstart

```
nextflow run -profile nbis TranscriptAssemblyHisat2Stringtie.nf \
  --reads '/path/to/reads*_{R1,R2}.fastq.gz' \
  --genome 'path/to/genome.fasta'
```

## Parameters

* `reads`: The path to the reads in quotes. Read pairs must include the `{}` notation to dictate what a read pair is.
* `genome`: The path to the genome assembly in quotes.
* `paired`: True if the reads are paired, or false if single end.
* `outdir`: The name of the results folder.

* `trimmomatic_adapter_path`: The path to the trimmomatic adapter file. (Default:'$TRIMMOMATIC_SHARE/adapters/TruSeq3-PE-2.fa' 
where '$TRIMMOMATIC_SHARE' is set in the Dockerfile).
* `trimmomatic_clip_options`: Read clipping options (Default:'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36').

* `hisat2_options`: Additional options for hisat2.

* `stringtie_options`: Additional options for stringtie.

## Stages


