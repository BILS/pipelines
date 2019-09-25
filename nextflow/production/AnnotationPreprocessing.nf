nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome_assembly = "$baseDir/test_data/test_assembly.fa"
params.outdir = "results"

log.info """\
 Annotation preprocessing workflow
 ===================================
 genome_assembly : ${params.genome_assembly}
 outdir          : ${params.outdir}
 """

include './../modules/annotation_modules'

workflow {

	main:
		annotation_preprocessing(params.genome_assembly)

	publish:

}
