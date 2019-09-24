// Annotation platform pipeline to perform various statistical operations on a GFF3-formatted file

// about title: "Pipeline to perform post-processing/statistical evaluation of GFF3-formatted annotation files"
//
// inputs "gff" : "Genome annotation file in GFF3 format"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	gff_get_trna
	gffread_extract_sequences
	gff_annotation_stats
	genome_tools_gff_stats
	genome_tools_gff_to_gtf
	gff_longest_cds
}

// run { "%.gff" * [ verify_generic.using(binary:"gffread") + verify_generic.using(binary:"gt") + sample_dir_prepare.using(sample_dir:true)
// 	+ gff_get_trna
// 	+ gffread_extract_sequences
// 	+ [ gff_annotation_stats,genome_tools_gff_stats,genome_tools_gff_to_gtf,gff_longest_cds	]
// ]
// }
