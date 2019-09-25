// Annotation platform pipeline to generate production data a genome

// about title: "A pipeline to generate all pre-annotation production output from a genome sequence"
//
// inputs "fa" : "Requires genome sequence in fasta format"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome_assembly = "$baseDir/test_data/test_assembly.fa"
params.outdir = "results"

log.info """\
 Annotation Preprocessing
 ===================================
 genome_assembly : ${params.genome_assembly}
 outdir          : ${params.outdir}
 """

include './../modules/annotation_modules'

workflow {
	fasta_filter_size(params.genome_assembly,1000)
	fasta_explode
	assembly_generate_stats
	bowtie2_index
	fastasplit
}


// run { "%.fa" * [ verify_annotation_preprocess + sample_dir_prepare.using(sample_dir:true)
// 	+ fasta_filter_size.using(size:1000,directory:"assembly") + [
// 		fasta_explode.using(directory:"scaffolds"),
// 		assembly_generate_stats,
// 		bowtie2_index.using(directory:"bowtie2-index") ,
// 		fastasplit
// 	]
// ] }
