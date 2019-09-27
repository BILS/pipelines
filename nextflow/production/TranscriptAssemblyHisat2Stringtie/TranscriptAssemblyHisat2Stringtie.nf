nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.reads = "$baseDir/test_data/*.fastq.gz"
params.genome = "$baseDir/test_data/genome.fa"
params.paired = true
params.outdir = "results"

params.trimmomatic_adapter_path = '$TRIMMOMATIC_SHARE/adapters/TruSeq3-PE-2.fa'
params.trimmomatic_clip_options = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

params.hisat2_options = ''

params.stringtie_options = ''

log.info """\

  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Transcript assembly using Hisat2/Stringtie workflow
 ===================================================

 General Parameters
     genome                     : ${params.genome}
     reads                      : ${params.reads}
     paired                     : ${params.paired}
     outdir                     : ${params.outdir}

 Trimmomatic parameters
     trimmomatic_adapter_path   : ${params.trimmomatic_adapter_path}
     trimmomatic_clip_options   : ${params.trimmomatic_clip_options}

 Hisat2 parameters
     hisat2_options             : ${params.hisat2_options}

 StringTie parameters
     stringtie_options          : ${params.stringtie_options}

 """

include './../workflows/annotation_workflows' params(params)

workflow {

	main:
	transcript_assembly_hisat2_stringtie(
		Channel.fromFilePairs(params.reads, checkIfExists: true)
		.ifEmpty { exit 1, "Cannot find reads matching ${params.reads}!\n" },
		Channel.fromPath(params.genome, checkIfExists: true)
		.ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" })

	publish:
	transcript_assembly_hisat2_stringtie.out.fastqc to: "${params.outdir}/fastqc"
	transcript_assembly_hisat2_stringtie.out.trimmomatic to: "${params.outdir}/trimmomatic"
	transcript_assembly_hisat2_stringtie.out.hisat2 to: "${params.outdir}/hisat2"
	transcript_assembly_hisat2_stringtie.out.stringtie to: "${params.outdir}/stringtie"
	transcript_assembly_hisat2_stringtie.out.multiqc to: "${params.outdir}/multiqc"

}

workflow.onComplete {
	log.info ( workflow.success ? "\nTranscript assembly complete!\n" : "Oops .. something went wrong\n" )
}
