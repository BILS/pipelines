// Annotation platform pipeline to generate production data a genome

// about title: "A pipeline to generate all pre-annotation production output from a genome sequence"
//
// inputs "fa" : "Requires genome sequence in fasta format"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	fasta_filter_size
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
