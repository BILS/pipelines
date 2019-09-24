// about title: "A pipeline to execute blastp and interproscan searches and pass them to Blast2Go"
//
// inputs "gff" : "Requires an annotation in GFF format as input (.gff)"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	gff2protein
	fastasplit
	blastp
	merge_blast_xml
	interpro
	merge_interpro_tsv
	merge_interpro_xml
	run_blast2go
	blast2go2gff
	interpro2gff
}

// run { "%.gff" * [
// 	 sample_dir_prepare.using(sample_dir:true) +
//         gff2protein.using(sample_dir:true) +
//         fastasplit.using(sample_dir:true) +
//         [
//                 ["%" * [blastp.using(sample_dir:true)] +
//                         "*.blastp" * [merge_blast_xml.using(sample_dir:true)]
//                 ],
//                 [ "%" * [interpro.using(sample_dir:true)] +
//                         [
//                                 ["*.tsv" * [merge_interpro_tsv.using(sample_dir:true)] ],
//                                 ["*.xml" * [merge_interpro_xml.using(sample_dir:true)] ]
//                         ]
//                 ]
//
//         ] +
//         run_blast2go.using(sample_dir:true) +
//         blast2go2gff.using(sample_dir:true) +
//         interpro2gff.using(sample_dir:true)
// 	]
// }
