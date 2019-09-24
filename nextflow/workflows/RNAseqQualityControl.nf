// Pipeline to align FASTQ files to a genome using Bowtie2/Tophat
// and quantification of samples using the Cufflinks package

// about title: "RNA-seq pipeline to perform quality control and trimming on RNA-seq read files"
//
// inputs "fq.gz" : "RNA-seq reads in gzipped fastq format, paired-end data expected in the format %_*.fq.gz."
//
// load 'pipeline.config'
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	trimmomatic
	fastqc
}

// run { "%_*.fq.gz" * [ sample_dir_prepare.using(sample_dir:true)
// 	+ trimmomatic.using(sample_dir:true)
// 	+ fastqc.using(sample_dir:true)
//
// ]
// }
