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
		
}
