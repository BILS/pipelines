nextflow.preview.dsl=2

include '../modules/annotation_modules' params(params)

workflow annotation_preprocessing {

	get:
		genome_assembly

	main:
		fasta_filter_size(genome_assembly)
		fasta_explode(fasta_filter_size.out)
		assembly_generate_stats(fasta_filter_size.out)
		bowtie2_index(fasta_filter_size.out)
		fastasplit(fasta_filter_size.out)

	emit:
		filtered_assembly = fasta_filter_size.out
		bowtie2_index = bowtie2_index.out
		single_fasta_dir = fasta_explode.out
		chunk_fasta_dir = fastasplit.out
		assembly_stats = assembly_generate_stats.out

}

workflow augustus_training_dataset {

	get:
		gff_annotation

	main:
		gff_filter_gene_models(gff_annotation)
		gff_longest_cds(gff_filter_gene_models.out)
		gff2protein(gff_longest_cds.out)
		blast_makeblastdb(gff2protein.out)
		blast_recursive(gff2protein.out,blast_makeblastdb.out)
		gff_filter_by_blast(gff_annotation,blast_recursive.out)
		gff2gbk(gff_filter_by_blast.out)
		gbk2augustus(gff2gbk.out)

	emit:
		dataset = gbk2augustus.out

}

workflow functional_annotation_input_preparation {

	get:
		gff_file

	main:
		gff2protein(gff_file)
		blastp(gff2protein.out.splitFasta(by: params.chunk_size))
		merge_blast_tab(blastp.out.collect())
		interpro(gff2protein.out.splitFasta(by: params.chunk_size))
		merge_interpro_tsv(interpro.out.collect())
		merge_interpro_xml(interpro.out.collect())

	emit:
		blast_results = merge_blast_tab.out
		interpro_tsv = merge_interpro_tsv.out
		interpro_xml = merge_interpro.xml.out

}

workflow transcript_assembly_hisat2_stringtie {

	get:
		reads
		genome

	main:
		fastqc(reads)
		trimmomatic(reads)
		hisat2_index(genome)
		hisat2(trimmomatic.out[0].mix(trimmomatic.out[2]),hisat2_index.out.collect())
		stringtie(hisat2.out)
		multiqc(fastqc.out.mix(trimmomatic.out).mix(hisat2.out).mix(stringtie.out))

	emit:
		fastqc = fastqc.out
		trimmomatic = trimmomatic.out
		hisat2 = hisat2.out
		stringtie = stringtie.out
		multiqc = multiqc.out

}
