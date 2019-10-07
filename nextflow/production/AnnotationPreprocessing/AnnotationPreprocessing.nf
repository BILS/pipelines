// nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome_assembly = "$baseDir/test_data/test_assembly.fa"
params.outdir = "results"

params.min_length = 1000

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

// include './../workflows/annotation_workflows' params(params)
//
// workflow {
//
// 	main:
// 	annotation_preprocessing(Channel.fromPath(params.genome_assembly, checkIfExists: true))
//
// 	publish:
// 	annotation_preprocessing.out.single_fasta_dir to: "${params.outdir}/scaffolds"
// 	annotation_preprocessing.out.chunk_fasta_dir to: "${params.outdir}/chunks"
// 	annotation_preprocessing.out.filtered_assembly to: "${params.outdir}/assembly"
// 	annotation_preprocessing.out.bowtie2_index to: "${params.outdir}/bowtie2-index"
// 	annotation_preprocessing.out.assembly_generate_stats to: "${params.outdir}/assembly_stats"
// }

// workflow annotation_preprocessing {
//
// 	get:
// 		genome_assembly
//
// 	main:
// 		fasta_filter_size(genome_assembly)
// 		fasta_explode(fasta_filter_size.out)
// 		assembly_generate_stats(fasta_filter_size.out)
// 		bowtie2_index(fasta_filter_size.out)
// 		fastasplit(fasta_filter_size.out)
//
// 	emit:
// 		filtered_assembly = fasta_filter_size.out
// 		bowtie2_index = bowtie2_index.out
// 		single_fasta_dir = fasta_explode.out
// 		chunk_fasta_dir = fastasplit.out
// 		assembly_stats = assembly_generate_stats.out
//
// }

Channel.fromPath(params.genome, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
    .set { genome_for_filter }

process fasta_filter_size {

    tag "Filtering Fasta ${fasta_file} by min length ${params.min_length}"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    file fasta_file from genome_for_filter

    output:
    file "${fasta_file.baseName}_min${params.min_length}.fasta" into genome_for_stats

    script:
    """
    seqtk seq -A $fasta_file -L ${params.min_length} > ${fasta_file.baseName}_min${params.min_length}.fasta
    """

}

process assembly_generate_stats {

    tag "Generating statistics for ${fasta_file.simpleName}"
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    file fasta_file from genome_for_stats

    output:
    file "${fasta_file.baseName}_assembly_report.txt"

    script:
    """
    fasta_statisticsAndPlot.pl --infile $fasta_file --output ${fasta_file.baseName}_assembly_report.txt
    """

}

workflow.onComplete {
    log.info ( workflow.success ? "\nAnnotation preprocessing complete!\n" : "Oops .. something went wrong\n" )
}
