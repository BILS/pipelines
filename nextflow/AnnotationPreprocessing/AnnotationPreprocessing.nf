// nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome_assembly = "$baseDir/test_data/test_assembly.fa"
params.outdir = "results"

params.min_length = 1000

log.info """
NBIS
 _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Annotation preprocessing workflow
 ===================================

 General parameters
     genome_assembly : ${params.genome_assembly}
     outdir          : ${params.outdir}

 Filtering parameters
     min_length      : ${params.min_length}

 """

// include './../workflows/annotation_workflows' params(params)
//
// workflow {
//
// 	main:
// 	annotation_preprocessing(Channel.fromPath(params.genome_assembly, checkIfExists: true))
//
// 	publish:
// 	annotation_preprocessing.out.filtered_assembly to: "${params.outdir}/assembly"
// 	annotation_preprocessing.out.assembly_generate_stats to: "${params.outdir}/assembly_stats"
// }

// workflow annotation_preprocessing {
//
// 	get:
// 		genome_assembly
//
// 	main:
// 		fasta_filter_size(genome_assembly)
// 		assembly_generate_stats(fasta_filter_size.out)
//
// 	emit:
// 		filtered_assembly = fasta_filter_size.out
// 		assembly_stats = assembly_generate_stats.out
//
// }

Channel.fromPath(params.genome_assembly, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
    .set { genome_for_filter }

process fasta_filter_size {

    tag "${fasta_file.baseName} ; min length = ${params.min_length}"
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

    tag "${fasta_file.simpleName}"
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    file fasta_file from genome_for_stats

    output:
    file "${fasta_file.baseName}_assembly_report"

    script:
    """
    fasta_statisticsAndPlot.pl --infile $fasta_file --output ${fasta_file.baseName}_assembly_report
    """
    // fasta_statisticsAndPlot.pl can be found in the NBIS GAAS repository
}

workflow.onComplete {
    log.info ( workflow.success ? "\nAnnotation preprocessing complete!\n" : "Oops .. something went wrong\n" )
}
