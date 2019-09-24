// A pipeline for running interproscan from protein fasta file

// about title: "A pipeline to execute interproscan searches"
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
	interpro ?
}

// run {  "%.fa" * [ verify_generic.using(binary:"fastasplit") + sample_dir_prepare.using(sample_dir:true) +
// 	fastasplit +
//         [
// 		[ "%" * [ interpro ] + [  "*.tsv" * [ merge_interpro_tsv] , "*.xml" * [ merge_interpro_xml ] ] ]
//
//         ]
// ]
// }
