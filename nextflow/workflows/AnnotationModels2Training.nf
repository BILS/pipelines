// Pipeline to generate augustus profile from GFF annotation

// about title: "Takes a genome anntation in GFF3 format and extracts data for training augustus profile"
//
// inputs "gff" : "A gene annotation file in GFF format"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	gff_filter_gene_models
	gff_longest_cds
	gff2protein
	blast_makeblastdb
	blast_recursive
	gff_filter_by_blast
	gff2gbk
	gbk2augustus
}


// run {  "%.gff" * [ verify_dependencies_annotation_models + sample_dir_prepare.using(sample_dir:true)
// 		+ gff_filter_gene_models
// 		+ gff_longest_cds
// 		+ gff2protein
// 		+ blast_makeblastdb
// 		+ blast_recursive.using(blast_outfmt:6)
// 		+ gff_filter_by_blast
// 		+ gff2gbk.using(flank:"500")
// 		+ gbk2augustus.using(test_size:100)
// 	]
// }
