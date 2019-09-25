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

	main:
		transcript_assembly_hisat2_stringtie(params.reads)

	publish:

}
