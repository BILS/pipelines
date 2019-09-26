nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "$baseDir/test_data/test.gff"
params.outdir = "results"

params.gff_gene_model_filter_options = '-c -r -d 500 -a 0.3'
params.codon_table = 1
params.test_size = 100

log.info """\
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

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

workflow.onComplete {
	log.info ( workflow.success ? "\nAugustus training dataset complete!\n" : "Oops .. something went wrong\n" )
}
