// A pipeline to create the required input for the functional annotation script gff3_sp_manage_functional_annotation

// about title: "A pipeline to execute blastp and interproscan searches to create input for the subsequent functional annotation script gff3_sp_manage_functional_annotation.pl."
//
// inputs "gff" : "Requires an annotation in GFF format as input (.gff)"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "$baseDir/test_data/test.gff"
params.outdir = "results"

log.info """\
 Functional annotation input preparation workflow
 ===================================
 gff_annotation : ${params.gff_annotation}
 outdir         : ${params.outdir}
 """

include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.gff_file)
	gff2protein
	fastasplit
	blastp
	merge_blast_tab
	interpro
	merge_interpro_tsv
	merge_interpro_xml
}

// run {  "%.gff*" * [ verify_generic.using(binary:"fastasplit") + sample_dir_prepare.using(sample_dir:true) +
// 	gff2protein +
// 	fastasplit +
//         [
//                 [ "%" * [ blastp.using(outfmt:6)] + merge_blast_tab ],
// 		[ "%" * [ interpro ] + [  "*.tsv" * [ merge_interpro_tsv] , "*.xml" * [ merge_interpro_xml ] ] ]
//
//         ]
// ]
// }
