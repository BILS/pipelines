// A pipeline to retrieve functional annotation from Interproscan and Blast

// about title: "A pipeline to execute interproscan searches and Blastp"
//
// inputs "fa" : "Requires a fasta file as input (.fa)"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	fastasplit
	blastp
	merge_blast_tab
	interpro
	merge_interpro_tsv
	merge_interpro_xml
}


// run {  "%.fa" * [ verify_generic.using(binary:"fastasplit") + sample_dir_prepare.using(sample_dir:true) +
// 	fastasplit +
//         [
// 		[ "%" * [ blastp.using(outfmt:6)] + merge_blast_tab ],
// 		[ "%" * [ interpro ] + [  "*.tsv" * [ merge_interpro_tsv] , "*.xml" * [ merge_interpro_xml ] ] ]
//
//         ]
// ]
// }
