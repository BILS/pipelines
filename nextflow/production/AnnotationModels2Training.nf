nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "$baseDir/test_data/test.gff"
params.outdir = "results"

log.info """\
 Augustus training dataset workflow
 ===================================
 gff_annotation : ${params.gff_annotation}
 outdir         : ${params.outdir}
 """

include './../workflows/annotation_workflows' params(params)

workflow {

	main:
	augustus_training_dataset(gff_annotation)

	publish:
	gbk2augustus.out.dataset to: "${params.outdir}/augustus_training_dataset"

}
