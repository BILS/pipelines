process assembly_generate_stats {

	tag "Generating statistics for ${fasta_file.simpleName}"

	input:
	path fasta_file

	output:
	path "${fasta_file.baseName}_assembly_report.txt"

	script:
	"""
	fasta_statisticsAndPlot.pl --infile $fasta_file --output ${fasta_file.baseName}_assembly_report.txt
	"""

}

process blastp {

	tag "Running blastp on ${query_fasta.simpleName} against ${blastp_dbpath}"

	input:
	path query_fasta
	path blastp_dbpath

	output:
	path "${query_fasta.baseName}_blast.tsv"

	script:
	"""
	blastp -db $blastp_dbpath -query $query_fasta -outfmt 6 \\
		-num_threads ${task.cpus} -out ${query_fasta.baseName}_blast.tsv
	"""

}

process merge_blast_tab {

	tag "Merging Blast results"

	input:
	path fragments

	output:
	path 'merged_blast_results.tsv'

	script:
	"""
	cat $fragments > merged_blast_results.tsv
	"""

}

process blast_makeblastdb {

	tag "Making Blast database: from: ${fasta_file.baseName} type: $dbtype"

	input:
	path fasta_file

	output:
	path "*.*"

	script:
	dbtype = "${params.dbtype}" ? "${params.dbtype}" ? 'prot'
	"""
	makeblastdb -in $fasta_file -dbtype $dbtype
	"""

}

process blast_recursive {

	tag "Performing recursive blast"

	input:
	path fasta_file
	path dbpath

	output:
	path "${fasta_file.baseName}_blast.tsv"

	script:
	"""
	blastp -query $fasta_file -db $db_path -num_threads ${task.cpus} \\
		-outfmt 6 -out ${fasta_file.baseName}_blast.tsv
	"""

}

process bowtie2_index {

	tag "Creating Bowtie2 Index: ${genome.baseName}"

	input:
	path genome

	output:
	path '*.bt2'

	script:
	"""
	bowtie2-build $genome ${genome.baseName}
	"""

}

process fasta_explode {

	// should not be necessary

	tag "Exploding Fasta"

	input:
	path input

	output:
	path output

	script:
	"""
	fastaexplode -d $output.dir -f $input 2> $output
	"""

}

process fasta_filter_size {

	tag "Filtering Fasta ${fasta_file} by min length ${params.min_length}"

	input:
	path fasta_file

	output:
	path "${fasta_file.baseName}_min${params.min_length}.fasta"

	script:
	"""
	seqtk seq -A $fasta_file -L ${params.min_length} > ${fasta_file.baseName}_min${params.min_length}.fasta
	"""

}

process fastasplit {

	// This should not be necessary.
	tag "Split fasta into chunks"

	input:
	path input

	output:
	path output

	script:
	"""
	fastasplit -f $input -o ${output.dir} -c $CHUNKS
	"""

}

process fastqc {

	tag "FastQC on $sample_id"

	input:
	tuple sample_id, path(reads)

	output:
	path "fastqc_${sample_id}_logs"

	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -t ${task.cpus} -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""

}

process gbk2augustus {

	tag "Converting Genbank to Augustus"

	input:
	path input

	output:
	path output

	script:
	"""
	randomSplit.pl $input ${params.test_size}
	"""

}

process gff2gbk {

	tag "Converting GFF to Genbank format"

	input:
	path gff_file
	path genome_gbk

	output:
	path "${gff_file.baseName}_f.gbk"

	script:
	"""
	gff2gbSmallDNA.pl $gff_file $genome_gbk ${params.flank_region_size} ${gff_file.baseName}_f.gbk
	"""

}

process gff2protein {

	tag "Converting GFF to protein sequence"

	input:
	path gff_file
	path genome_fasta

	output:
	path "${gff_file}_proteins.fasta"

	script:
	"""
	TMP_FASTA = $(mktemp -u --suffix ".fa" )
	gff3_sp_extract_sequences.pl -o \$TMP_FASTA -f $genome_fasta \\
		-p -cfs -cis -ct ${params.codon_table} --gff $gff_file
	fix_fasta.rb \$TMP_FASTA > ${gff_file}_proteins.fasta
	"""

}

process gff_filter_by_blast {

	tag "Filtering GFF by Blast results (outfmt:6)"

	input:
	path gff_file
	path blast_file

	output:
	path "${gff_file.baseName}_blast-filtered.gff3"

	script:
	"""
	gff_filter_by_mrna_id.pl --gff $gff_file --blast $blast_file \\
		--outfile ${gff_file.baseName}_blast-filtered.gff3
	"""

}

process gff_filter_gene_models {

	tag "Filter gene models by GFF"

	input:
	path gff3_file
	path genome_fasta

	output:
	path "${gff_file.baseName}_model-filtered.gff3"

	script:
	"""
	filter_sort.pl -f $gff3_file -F $genome_fasta \\
		-o ${gff_file.baseName}_model-filtered.gff3 ${params.gff_gene_model_filter_options}
	"""
}

process gff_longest_cds {

	tag "Retaining longest CDS sequences"

	input:
	path gff3_file

	output:
	path "${gff3_file.baseName}_longest_cds.gff3"

	script:
	"""
	find_longest_CDS.pl -f $gff3_file -o ${gff3_file.baseName}_longest_cds.gff3
	"""
}

process hisat2_index {

		tag "Indexing ${genome}"

		input:
		path genome_fasta

		output:
		tuple "${genome_fasta.baseName}.hisat2_index", path('*.ht2')

		script:
		"""
		hisat2-build -p ${task.cpus} $genome_fasta ${genome_fasta.baseName}.hisat2_index
		"""
}

process hisat2 {

	tag "Aligning reads (${sample_id}) to genome"

	input:
	tuple sample_id, path(reads)
	tuple hisat2_basename, path(hisat2_index_files)

	output:
	path "${sample_id}_sorted_alignment.bam"
	path 'splicesite.txt'

	script:
	if (params.paired){
		"""
		hisat2 ${params.hisat2_options} --novel-splicesite-outfile splicesite.txt
			-p ${task.cpus} -x $hisat2_basename -1 ${reads[0]} -2 ${reads[1]} | \\
			samtools sort -@ ${task.cpus} -o sorted_alignment.bam -
		"""
	} else {
		"""
		hisat2 ${params.hisat2_options} --novel-splicesite-outfile splicesite.txt
			-p ${task.cpus} -x $hisat2_basename -U $reads | \\
			samtools sort -@ ${task.cpus} -o sorted_alignment.bam -
		"""
	}

}


process interpro {

	tag "Running InterProScan: Protein function classification"

	input:
	path protein_fasta
	path interpro_dbpath

	output:
	path '*.gff3'
	path '*.xml'
	path '*.tsv'

	script:
	"""
	interproscan $interpro_dbpath -i $protein_fasta -d $output -iprlookup -goterms -pa -dp
	"""

}

process merge_interpro_xml {

	tag "Merging InterProScan XML results"

	input:
	path xml_files

	output:
	path 'interpro_search.xml'

	// This code is not robust at all. Need to rewrite (e.g. the -v "xml" already excludes "protein-matches" lines because of the "xmlns" attributes)
	script:
	"""
	head -n 2 ${xml_files[0]} > interpro_search.xml
	for XML in $xml_files; do
		grep -v "xml" \$XML | grep -v "protein-matches" >> interpro_search.xml
	done
	tail -n 1 ${xml_files[0]} >> interpro_search.xml
	"""
}

process merge_interpro_tsv {

	tag "Merging InterProScan TSV results"

	input:
	path tsv_files

	output:
	path 'interpro_search.tsv'

	script:
	"""
	cat $tsv_files > interpro_search.tsv
	"""

}

process stringtie {

	tag "StringTie Transcript assembly from RNA-seq alignments: ${sorted_bam_file.name}"

	input:
	path sorted_bam_file

	output:
	path 'transcripts.gtf'

	script:
	"""
	stringtie ${sorted_bam_file} -l ${sorted_bam_file.name} -o transcripts.gtf -p ${task.cpus} ${params.stringtie_options}
	"""

}

process trimmomatic {

	tag "Adapter-trimming reads: ${sample_id}"

	input:
	tuple sample_id, path(reads)

	output:
	tuple sample_id, path('*_paired_*.fastq.gz')
	tuple sample_id, path('*_unpaired_*.fastq.gz')
	tuple sample_id, path('*_trimmed.fastq.gz')

	script:
	if (params.paired) {
		"""
		java -jar ${params.trimmomatic_jar} PE -threads ${task.cpus} $reads \\
		 	${sample_id}_paired_1.fastq.gz ${sample_id}_unpaired_1.fastq.gz \\
			${sample_id}_paired_2.fastq.gz ${sample_id}_unpaired_2.fastq.gz \\
			ILLUMINACLIP:$${params.trimmomatic_adapter_path}:2:30:10 \\
			${params.trimmomatic_clip_options}
		"""
	} else {
		"""
		java -jar ${params.trimmomatic_jar} SE -threads ${task.cpus} $reads \\
		 	${sample_id}_trimmed.fastq.gz \\
			ILLUMINACLIP:$${params.trimmomatic_adapter_path}:2:30:10 \\
			${params.trimmomatic_clip_options}
		"""
	}

}
