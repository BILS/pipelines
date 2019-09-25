// A pipeline to create to assemble transcripts for annotation

// about title: "A pipeline to assemble transcripts from RNAseq data based on Hisat2/Stringtie (input ...)"
//
// //inputs "f*q.[gz]" : "Requires FastQ file(s) as input"
//
// load 'pipeline.config'
nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.reads = "$baseDir/test_data/reads*.fastq.gz"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

log.info """\
 Transcript assembly using Hisat2/Stringtie workflow
 ===================================================
 reads          : ${params.reads}
 outdir         : ${params.outdir}
 multiqc        : ${params.multiqc}
 """

include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	trimmomatic
	hisat2
	samtools_sam_to_bam
	samtools_sort_bam
	stringtie
}

// run { ~"(.*)_.+.f*q*[.gz]?" *
// 	[ verify_generic.using(binary:"hisat2")  + verify_generic.using(binary:"samtools")  + verify_generic.using(binary:"stringtie") + sample_dir_prepare.using(sample_dir:true) +
// 		trimmomatic +
// 		hisat2 +
// 		samtools_sam_to_bam +
// 		samtools_sort_bam +
// 		stringtie
// 	]
// }
