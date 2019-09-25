nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.reads = "$baseDir/test_data/reads*.fastq.gz"
params.outdir = "results"

log.info """\
 Transcript assembly using Hisat2/Stringtie workflow
 ===================================================
 reads          : ${params.reads}
 outdir         : ${params.outdir}
 """

include './../modules/annotation_modules'

workflow {

	main:
	transcript_assembly_hisat2_stringtie(params.reads)

	publish:
	transcript_assembly_hisat2_stringtie.out.fastqc to: "${params.outdir}/fastqc"
	transcript_assembly_hisat2_stringtie.out.trimmomatic to: "${params.outdir}/trimmomatic"
	transcript_assembly_hisat2_stringtie.out.hisat2 to: "${params.outdir}/hisat2"
	transcript_assembly_hisat2_stringtie.out.stringtie to: "${params.outdir}/stringtie"
	transcript_assembly_hisat2_stringtie.out.multiqc to: "${params.outdir}/multiqc"

}
