nextflow.preview.dsl=2

include '../modules/annotation_modules' params(params)

workflow annotation_preprocessing {

	get:
		genome_assembly

	main:
		fasta_filter_size(genome_assembly,1000)
		fasta_explode
		assembly_generate_stats
		bowtie2_index
		fastasplit

	emit:

	// run { "%.fa" * [ verify_annotation_preprocess + sample_dir_prepare.using(sample_dir:true)
	// 	+ fasta_filter_size.using(size:1000,directory:"assembly") + [
	// 		fasta_explode.using(directory:"scaffolds"),
	// 		assembly_generate_stats,
	// 		bowtie2_index.using(directory:"bowtie2-index") ,
	// 		fastasplit
	// 	]
	// ] }

}

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

}

workflow functional_annotation_input_preparation {

	get:
		gff_file

	main:
		gff2protein
		fastasplit
		blastp
		merge_blast_tab
		interpro
		merge_interpro_tsv
		merge_interpro_xml

	emit:

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

}

workflow transcript_assembly_hisat2_stringtie {

	get:
		reads
		genome

	main:
		trimmomatic
		hisat2
		samtools_sam_to_bam
		samtools_sort_bam
		stringtie

	emit:

	// run { ~"(.*)_.+.f*q*[.gz]?" *
	// 	[ verify_generic.using(binary:"hisat2")  + verify_generic.using(binary:"samtools")  + verify_generic.using(binary:"stringtie") + sample_dir_prepare.using(sample_dir:true) +
	// 		trimmomatic +
	// 		hisat2 +
	// 		samtools_sam_to_bam +
	// 		samtools_sort_bam +
	// 		stringtie
	// 	]
	// }

}
