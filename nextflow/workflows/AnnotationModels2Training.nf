/* A workflow to generate an augustus profile from GFF annotation

about title: "Takes a genome anntation in GFF3 format and extracts data for training augustus profile"

 inputs "gff" : "A gene annotation file in GFF format"
*/
nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "$baseDir/test_data/test.gff"
params.outdir = "results"

log.info """\
 Augustus training dataset workflow
 ===================================
 gff_annotation : ${params.gff_annotation}
 outdir         : ${params.outdir}
 """

include './../modules/annotation_modules' params(params)

workflow augustus_training_dataset {

	get:
		gff_annotation

	main:
		gff_filter_gene_models(gff_annotation)
		gff_longest_cds(gff_filter_gene_models.out)
		gff2protein(gff_longest_cds.out)
		blast_makeblastdb(gff2protein.out,"prot")
		blast_recursive(gff2protein.out,blast_makeblastdb.out,6)
		gff_filter_by_blast(gff_annotation,blast_recursive.out)
		gff2gbk(gff_filter_by_blast.out,500)
		gbk2augustus(gff2gbk.out,100)

	emit:
		dataset = gbk2augustus.out

	publish:
		gbk2augustus.out to: 'results/augustus_training_dataset'
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
