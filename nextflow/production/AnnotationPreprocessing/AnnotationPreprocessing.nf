nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome_assembly = "$baseDir/test_data/test_assembly.fa"
params.outdir = "results"

params.chunk_size = 200

log.info """\
 _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Annotation preprocessing workflow
 ===================================
 genome_assembly : ${params.genome_assembly}
 outdir          : ${params.outdir}
 """

include './../workflows/annotation_workflows' params(params)

workflow {

	main:
	annotation_preprocessing(Channel.fromPath(params.genome_assembly, checkIfExists: true))

	publish:
	annotation_preprocessing.out.single_fasta_dir to: "${params.outdir}/scaffolds"
	annotation_preprocessing.out.chunk_fasta_dir to: "${params.outdir}/chunks"
	annotation_preprocessing.out.filtered_assembly to: "${params.outdir}/assembly"
	annotation_preprocessing.out.bowtie2_index to: "${params.outdir}/bowtie2-index"
	annotation_preprocessing.out.assembly_generate_stats to: "${params.outdir}/assembly_stats"
}

workflow.onComplete {
	log.info ( workflow.success ? "\nAnnotation preprocessing complete!\n" : "Oops .. something went wrong\n" )
}
