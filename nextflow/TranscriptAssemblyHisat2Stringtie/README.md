# Transcript assembly pipeline

## Quickstart

```
nextflow run -profile nbis,docker TranscriptAssemblyHisat2Stringtie.nf \
  --reads '/path/to/reads*_{R1,R2}.fastq.gz' \
  --genome 'path/to/genome.fasta'
```

## Parameters

* General
    * `reads`: The path to the reads in quotes. The read pairs path must use the `{}` notation to define what a read pair is.
    * `genome`: The path to the genome assembly in quotes.
    * `paired`: True if the reads are paired, or false if single end.
    * `outdir`: The name of the results folder.
* Trimmomatic parameters
    * `trimmomatic_adapter_path`: The path to the trimmomatic adapter file. (Default:'$TRIMMOMATIC_SHARE/adapters/TruSeq3-PE-2.fa' 
    where '$TRIMMOMATIC_SHARE' is set in the Dockerfile).
    * `trimmomatic_clip_options`: Read clipping options (Default:'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36').
* Hisat2 parameters.
    * `hisat2_options`: Additional options for hisat2. E.g. strandedness (`--hisat2_options '--fr'`). 
    See the [Hisat2 Manual]() for the full range of options.
* StringTie parameters
    * `stringtie_options`: Additional options for stringtie.

## Stages

1.
    1. **FastQC**: Reads properties are summarised and checked for common issues relating to sequence content and quality.
    2. **Trimmomatic**: Reads are trimmed for adapter read-through and low quality regions.
    3. **Hisat2 Build**: Builds an index database for the input genome.
2. **Hisat2**: Align trimmed reads to the genome.
3. **StringTie**: Assemble transcripts from the aligned reads.
4. **MultiQC**: Provide a consolidated report of the statistics from each previous tool. 
