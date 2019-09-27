nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "$baseDir/test_data/test.gff"
params.outdir = "results"

log.info """\
 _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Functional annotation input preparation workflow
 ===================================
 gff_annotation : ${params.gff_annotation}
 outdir         : ${params.outdir}
 """

include './../workflows/annotation_workflows' params(params)

workflow {

	main:
	functional_annotation_input_preparation(Channel.fromPath(params.gff_file, checkIfExists: true))

	publish:
	functional_annotation_input_preparation.out to: "${params.outdir}"
}

workflow.onComplete {
	log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}
