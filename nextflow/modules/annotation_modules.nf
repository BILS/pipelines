
process assembly_generate_stats {

	input:
	path fasta_file

	output:
	path assembly_report

	script:
	"""
	fasta_statisticsAndPlot.pl --infile $fasta_file --output $assembly_report
	"""

}

process blastp {

	input:
	path query_fasta
	path blastp_dbpath
	val outfmt // 5
	val evalue

	output:
	path blast_results

	script:
	"""
	blastp -db $blastp_dbpath -query $query_fasta -outfmt $outfmt \
		-num_threads ${task.cpus} $evalue -out $blast_results
	"""

}

process merge_blast_tab {

	input:
	path fragments

	output:
	path output //branch.sample + "_blast.out"

	script:
	"""
	cat $fragments > $output
	"""

}

process blast_makeblastdb {

	input:
	path input
	val dbtype // "prot"

	output:
	path output

	script:
	"""
	makeblastdb -in $input -dbtype $dbtype
	"""

}

process blast_recursive {

	input:
	path fasta_file
	path dbpath
	val outfmt // 5

	output:
	path output

	script:
	"""
	blastp -query $fasta_file -db $db_path -num_threads ${task.cpus} -outfmt $outfmt -out $output
	"""

}

process bowtie2_index {

	input:
	tuple species_id, path(genome)

	output:
	path '*.bt2'

	script:
	"""
	bowtie2-build $genome $species_id
	"""

}

process fasta_explode {

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

	input:
	path fasta_file
	val min_length

	output:
	path 'filtered.fasta'

	script:
	"""
	# fasta_filter_size.rb -i $fasta_file -s $min_length -o $output
	seqtk seq -A $fasta_file -L $min_length > filtered.fasta
	"""

}

process fastasplit {

	// REPLACE with channel operator

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

	tag "FASTQC on $sample_id"
    publishDir params.outdir

	input:
	tuple sample_id, path(reads)

	output:
	path "fastqc_${sample_id}_logs"

	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""

}

process gbk2augustus {

	input:
	path input
	val test_size // 100

	output:
	path output

	script:
	"""
	randomSplit.pl $input $test_size
	"""

}

process gff2gbk {

	input:
	path gff_file
	path genome
	val flank_region_size // 1000

	output:
	path output

	script:
	"""
	gff2gbSmallDNA.pl $gff_file $genome $flank_region_size $output
	"""

}

process gff2protein {

	input:
	path input

	output:
	path output

	script:
	"""
	gff3_sp_extract_sequences.pl -o $input.prefix.tmp \
		-f $GENOME_FA -p -cfs -cis -ct $CODON_TABLE --gff $input && \
		$BPIPE_BIN/fix_fasta.rb $input.prefix.tmp.fa > $output && \
		rm $input.prefix.tmp.fa
	"""

	// doc about: "A module to extract protein sequences from GFF annotations",
	// description: "Reports protein sequences from GFF annotations",
	// author: "marc.hoeppner@bils.se, jacques.dainat@bils.se"
	//
	// var directory : "protein"
	//
	// if (branch.sample_dir) {
	// 	output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
	// } else {
	// 	if (directory.length() > 0) {
	// 		output.dir = directory
	// 	}
	// }

	// requires here
	// requires GENOME_FA : "Must provide a genome sequence (GENOME_FA)"
	// requires CODON_TABLE :  "Must specify the translation table (CODON_TABLE)"

	// Running a command
	// transform(".gff") to (".proteins.fa") {
	// 	exec "$BPIPE_BIN/gff3_sp_extract_sequences.pl -o $input.prefix"+".tmp -f $GENOME_FA -p -cfs -cis -ct $CODON_TABLE --gff $input && $BPIPE_BIN/fix_fasta.rb $input.prefix"+".tmp.fa > $output && rm $input.prefix"+".tmp.fa"
	// }

}

process gff_annotation_stats {

	input:
	path input

	output:
	path output

	script:
	"""
	$KAHARI_ANNOTATION_STATS --gff $input --fasta $GENOME_FA > $output
	"""

	// doc about: "A module to generate annotation statistics, based on the code by A. Kähäri",
	// author: "marc.hoeppner@bils.se"
	//
	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// requires KAHARI_ANNOTATION_STATS : "Must provide path to annotation_stats.pl (KAHARI_ANNOTATION_STATS)"
	// requires GENOME_FA : "Must provide path to genome sequence in FASTA format (GENOME_FA)"
	//
	// if (sample_dir) { output.dir = branch.outdir }

	// transform(".gff") to (".stats") {
	// 	exec "$KAHARI_ANNOTATION_STATS --gff $input --fasta $GENOME_FA > $output"
	// }

}

process gff_filter_by_blast {

	input:
	path input

	output:
	path output

	script:
	"""
	perl $BPIPE_BIN/gff_filter_by_mrna_id.pl --gff $input --blast $blast_file --outfile $output
	"""

	// doc about: "Remove GFF entries if they are found in a recursive blast (redundancy removal)",
	// description: "Takes a recursive blast input (format 6) and removes redundant entries from GFF file",
	// constraints: "Requires the blast file to be in outfmt 6",
	// author: "marc.hoeppner@bils.se"
	//
	// var directory : "nonredundant"
	//
	// if (branch.sample_dir) {
	// 	output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
	// } else {
	// 	if (directory.length() > 0) {
	// 		output.dir = directory
	// 	}
	// }

	// requires here

	// Running a command

	// blast_file = input
	// input will change due to the expression below, so we capture it here

	// produce(branch.name + ".nr.gff") {
	// 	from("*.gff") {
	// 		exec "perl $BPIPE_BIN/gff_filter_by_mrna_id.pl --gff $input --blast $blast_file --outfile $output"
	// 	}
	// }
}

process gff_filter_gene_models {

	input:
	path input

	output:
	path output

	script:
	"""
	perl $KAHARI_GIT/scripts/GFF/filter_sort.pl -f $input -F $GENOME_FA -o $output $options
	"""

	// doc about: "Wrapper module around A. Kaharis' gff filter_sort.pl script",
	// description: "Filters transcript models by a set of criteria",
	// constraints: "Requires the Kahari code base to be present",
	// author: "marc.hoeppner@bils.se"
	//
	// var directory : "filter"
	//
	// var filter_c : "-c"
	// var filter_a : "0.3"
	// var filter_r : "-r"
	// var filter_d : "500"
	//
	// def options = filter_c + " " + filter_r + " -d " + filter_d + " -a " + filter_a
	//
	// if (branch.sample_dir) {
	// 	output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
	// } else {
	// 	if (directory.length() > 0) {
	// 		output.dir = directory
	// 	}
	// }
	//
	// // requires here
	// requires KAHARI_GIT : "Must set path to Kahari git directory (KAHARI_GIT)"
	// requires GENOME_FA : "Must provide a genome sequence in FASTA format"
	//
	// // Running a command
	//
	// filter("filter") {
	// 	exec "perl $KAHARI_GIT/scripts/GFF/filter_sort.pl -f $input -F $GENOME_FA -o $output $options"
	// }
}

process gff_get_trna {

	input:
	path input

	output:
	path output

	script:
	"""
	grep "trnascan" -v $input.gff > $output
	"""

	// doc about: "Simple operation to excise tRNAs from a gff annotation",
	// constraints: "Expects a the annotation to come from Maker and contain mRNAs and tRNAs, nothing else",
	// author: "marc.hoeppner@bils.se"
	//
	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) { output.dir = branch.outdir}

	// transform(".gff") to (".mRNA.gff") {
	// 	exec "grep \"trnascan\" -v $input.gff > $output"
	// }

}

process gff_longest_cds {

	input:
	path input

	output:
	path output

	script:
	"""
	$KAHARI_GIT/scripts/GFF/find_longest_CDS.pl -f $input -o $output
	"""

	// doc about: "A module to limit gene loci in a GFF3-formatted annotation to the longest product",
	// description: "Parses a GFF3-formatted annotation and filters out all transcripts except those with the longest CDS",
	// author: "marc.hoeppner@bils.se"

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = branch.sample_dir }
	//
	// requires KAHARI_GIT : "Must set path to Kahari git directory (KAHARI_GIT)"
	//

	// Defining output directory
	// if (sample_dir) {
	// 	output.dir = branch.outdir
	// }

	// transform(".gff") to (".longest_cds.gff") {
	// 	exec "$KAHARI_GIT/scripts/GFF/find_longest_CDS.pl -f $input -o $output"
	// }

}

process gff_stable_id {

	input:
	path input

	output:
	path output

	script:
	"""
	$BPIPE_BIN/gff_create_stable_id.pl --gff $input --id_trunk $ID_TRUNK --outfile $output
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc title: "Generates stable IDs for all features in the GFF file"
	//
	// requires ID_TRUNK : "Must set variable ID_TRUNK"

	// if (sample_dir) {
	// 	output.dir = branch.outdir
	// }

	// filter("stable_id") {
	// 	exec "$BPIPE_BIN/gff_create_stable_id.pl --gff $input --id_trunk $ID_TRUNK --outfile $output"
	// }

	// branch.gff_file_with_ids output

}

process gffread_extract_sequences {

	input:
	path input

	output:
	path output

	script:
	"""
	gffread -x $output1 -g $GENOME_FA $input
	gffread -y $output2 -g $GENOME_FA $input
	"""

	// doc about: "A module to extract CDS and protein sequences from an annotation using gffread",
	// author: "marc.hoeppner@bils.se"
	//
	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// requires GENOME_FA : "Must provide genome sequence in FASTA format (GENOME_FA)"
	//
	// if (sample_dir) { output.dir = branch.outdir }
	//
	// produce(input.prefix+".cds.fa",input.prefix+".proteins.fa") {
	// 	exec """
	// 	gffread -x $output1 -g $GENOME_FA $input ; gffread -y $output2 -g $GENOME_FA $input
	// 	"""
	// }
}


// Hisat2 module
process hisat2 {

	input:
	path input

	output:
	path output

	script:
	if (PAIRED.toBoolean()){
	"""
	hisat2 $options --novel-splicesite-outfile $output.dir/splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -1 $input1 -2 $input2 >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5
	"""
	} else {
	"""
	hisat2 $options --novel-splicesite-outfile $output.dir/splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -U $input >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5
	"""
	}

	// doc title: "Align RNA-seq reads against a reference with hisat2",
	// desc: """
	// Uses hisat2/Bowtie2 to align reads against a genome
	//
	// Stage options:
	// sample_dir : create sample-specific output folrder (boolean)
	//
	// Required variables are:
	// PAIRED : paired-end data (boolean)
	// HISAT2 : this should point to the location of the tophat executable
	// LIBRARY_METHOD : This specifies the type of sequencing library: fr/rf/ff
	// PHRED_SCORE : Specifies the quality encoding, e.g. usually this will be solexa-quals/phred33/phred64
	// HISAT2_INDEX : The location of a hisat2-formatted genome index
	// // SAMTOOLS : This specifies the path to samtools
	// //  BOWTIE2 : This specifies the path to bowtie2
	// """,
	// constraints: """
	// Works with fq.gz and fq input files, but assumes paired-end
	// reads (if not, set paired to false).
	// The first argument is expected to be the left mate, the second
	// argument needs to be the right mate
	// """,
	// author: "jacques.dainat@nbis.se"
	//
	// // Exposed options
	// var sample_dir : false
	//
	// options = ""
	// // Handle --rna-strandness option
	// LIBRARY_METHOD=LIBRARY_METHOD.toLowerCase()
	// if (LIBRARY_METHOD == "r" || LIBRARY_METHOD == "rf" || LIBRARY_METHOD == "fr-firstrand" || LIBRARY_METHOD == "fr-firststrand") {
	// 	if (PAIRED.toBoolean()) {
	// 		options += "--rna-strandness RF"
	// 	}
	// 	else{ // single end reads
	// 		options += "--rna-strandness R"
	// 	}
	// }
	// else if (LIBRARY_METHOD == "f" || LIBRARY_METHOD == "fr" || LIBRARY_METHOD == "fr-secondstrand" ) {
	// 	if (PAIRED.toBoolean()) {
	// 		options += "--rna-strandness FR"
	// 	}
	// 	else{ // single end reads
	// 		options += "--rna-strandness F"
	// 	}
	// }
	// else{ // default value is unstranded
	// 	println "LIBRARY_METHOD is unstranded"
	// }
	//
	// // Handle Quality score
	// if (PHRED_SCORE == "phred33") {
	// 	options += " --phred33"
	// }
	// else if (PHRED_SCORE == "phred64") {
	// 	options += " --phred64"
	// }
	// else if (PHRED_SCORE == "solexa-quals") {
	// 	options += " --solexa-quals"
	// }
	//
	// // Configuring the output directory
	// if (branch.sample_dir) { sample_dir = true }
	//
	// requires HISAT2 : "Must set the TOPHAT variable to point to tophat location"
	// requires LIBRARY_METHOD : "Must specify a sequencing library method (LIBRARY_METHOD)"
	// requires PHRED_SCORE : "Must set a phred score (PHRED_SCORE)"
	// requires HISAT2_INDEX : "Must specify a Bowtie2 index (BWT2_INDEX)"
	// requires PAIRED : "Must specify if the sample is stranded or not (true or false)"
	//
	// // We subsequently need to keep track of folders
	// // Here we set a name accessible to all subsequent modules.
	// if (sample_dir) {
	// 	output.dir = branch.outdir + "/hisat2"
	// } else {
	// 	output.dir = "hisat2/" + branch.sample
	// }
	//
	// // If a basename for this branch was set further upstream
	//
	// produce("accepted_hits.sam") {
	// 	uses(threads:16) {
	// 		if (PAIRED.toBoolean()) {
	// 			exec "$HISAT2 $options --novel-splicesite-outfile $output.dir/splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -1 $input1 -2 $input2 >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5","hisat2"
	// 		} else {
	// 			exec "$HISAT2 $options --novel-splicesite-outfile $output.dir/splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -U $input >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5","hisat2"
	// 		}
	// 	}
	// }
	//
	// check {
	// 	exec "[ -s $output ]"
	// } otherwise {
	// 	succeed "The HISAT2 output ($output) is empty. Stopping this branch ($branch.name)"
	// }
}

process htseq_count {

	input:
	path input

	output:
	path output

	script:
	"""
	HTSEQ_COUNT -f bam -t exon $input.bam $GENOME_GTF > $output
	"""

	// doc title: "HTseq-count counts reads overlapping an annotation",
	//
	// desc: """
	// HTseq-count is a python tool designed to accept a read alignment in BAM or SAM
	// format and counts the reads overlapping annotated exons provided in GTF format.
	// """,
	//
	// constraints: "Requires a read alignment in BAM/SAM format and GTF-formatted annotation",
	//
	// author: "marc.hoeppner@bils.se"
	//
	// var sample_dir : true
	//
	// requires GENOME_GTF : "Must provide an annotation in GTF format (GENOME_GTF)"
	// requires HTSEQ_COUNT : "Must provide path to htseq-count (HTSEQ_COUNT)"
	//
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) {
	// 	output.dir = branch.outdir + "/htseq_count"
	// } else {
	// 	output.dir = "htseq_count"
	// }
	//
	// transform(".bam") to (".htseqcount") {
	// 	exec "HTSEQ_COUNT -f bam -t exon $input.bam $GENOME_GTF > $output"
	// }

}

process igvcount {

	input:
	path input

	output:
	path output

	script:
	"""
	igvtools count $input.bam $output hg19
	"""

	// exec "igvtools count $input.bam $output hg19"
}


process interpro {

	input:
	path input

	output:
	path output

	script:
	"""
	interproscan $db_list -i $input -d ${output.dir} -iprlookup -goterms -pa -dp > /dev/null 2> /dev/null
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc "Runs a protein fasta file against InterPro"
	//
	// if (sample_dir) { output.dir = branch.outdir }
	//
	// requires INTERPROSCAN : "Must specify the location of interproscan (INTERPROSCAN)"
	// requires INTERPRO_DB_LIST : "Must specify the list of DB to use (INTERPRO_DB_LIST)"
	//
	// // set dblist if needed - if db_list is empty, all the DB will be used
	// var db_list : ""
	// if (INTERPRO_DB_LIST.toLowerCase() != "all"){
	// 	db_list = "-appl $INTERPRO_DB_LIST"
	// }

	// produce(input+".gff3",input+".tsv",input+".xml") {
	// 	exec "$INTERPROSCAN $db_list -i $input -d ${output.dir} -iprlookup -goterms -pa -dp > /dev/null 2> /dev/null ","interpro"
	// }
}


process merge_interpro_xml {

	input:
	path input

	output:
	path output

	script:
	"""
	head -n 2 $first_file > $temp_file
	for (i in inputs) {
		grep -v xml $i | grep -v protein-matches >> $temp_file
	}
	tail -n 1 $first_file >> $temp_file
	mv $temp_file $output
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc "Merges the XML output from multiple InterPro searches"
	//
	// if (sample_dir) { output.dir = branch.outdir }
	//
	// produce(branch.sample+"_interpro.xml") {

	// This is a really stupid hack to merge XML files
	// We take the first lines and the last line from the first file
	// And squeeze the other stuff from all files in between

	// def first_file = inputs[0]
	// def temp_file = branch.sample + "_interpro.tmp"
	//
	// exec "head -n 2 $first_file > $temp_file"
	//
	// for (i in inputs) {
	// 	exec "grep -v xml $i | grep -v protein-matches >> $temp_file"
	// }
	//
	// exec "tail -n 1 $first_file >> $temp_file"
	//
	// exec "mv $temp_file $output"
	// }
	//
	// branch.ipr_xml = output
}

process merge_interpro_tsv {

	input:
	path input

	output:
	path output

	script:
	"""
	cat $inputs >> $output
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc "Merges the TSV output from multiple InterPro searches"
	//
	// if (sample_dir) { output.dir = branch.outdir }
	//
	// produce(branch.sample+"_interpro.tsv") {

	// 	exec "cat $inputs >> $output"
	//
	// }

	// We save the output name since we need it for a later pipeline stage.

	// branch.iprtsv = output

	// Pass on the name of the output file
}


process interpro2gff {

	input:
	path input

	output:
	path output

	script:
	"""
	ipr_update_gff $input $iprtsv > $output
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc "Converts InterPro results into GFF format"
	//
	// if (sample_dir) { output.dir = branch.outdir }


	// This is a horrible way of doing this...
	// We need the name of a file set in an earlier stage,
	// but the variable doesn't carry over, so:

	// if (sample_dir) {
	// 	iprtsv = branch.outdir + "/" + branch.sample + "_interpro.tsv"
	// } else {
	// 	iprtsv = branch.sample + "_interpro.tsv"
	// }

	// filter("interpro") {
	// 	exec "ipr_update_gff $input $iprtsv > $output"
	// }
}

// Maker module

process cufflinks2maker {

	input:
	path input

	output:
	path output

	script:
	"""
	cufflinks2gff3 $input > $output
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc "Converts cufflinks-reconstructed transcripts to Maker-compatible GFF3 format"
	//
	// if (sample_dir) {
	// 	output.dir = branch.outdir + "/maker_tracks"
	// } else {
	// 	output.dir = "maker_tracks"
	// }
	//
	// transform("gtf") to (assembly_method + ".maker.gff3") {
	// 	exec "cufflinks2gff3 $input > $output"
	// }

}

process dedup {

	input:
	path input

	output:
	path output

	script:
	"""
	java -Xmx6g -Djava.io.tmpdir=$TMP -jar $PICARD_HOME/MarkDuplicates.jar
	INPUT=$input.bam
	REMOVE_DUPLICATES=true
	VALIDATION_STRINGENCY=LENIENT
	AS=true
	METRICS_FILE=${output}.mark_duplicates.log
	OUTPUT=$output.bam
	2> /dev/null && md5sum $output.bam > ${output}.md5sum
	"""

	// doc "Marks duplicate reads in BAM file using Picard"
	//
	// var sample_dir : false
	//
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/bam"
	// } else {
	// output.dir="bam/${branch.sample}"
	// }
	//
	// requires TMP : "Must specify a temporary directory (TMP)"
	// requires PICARD_HOME : "Must provide path to picard (PICARD_HOME)"
	//
	// transform("bam") to ("dedup.bam") {
	//
	// exec """
	// java -Xmx6g -Djava.io.tmpdir=$TMP -jar $PICARD_HOME/MarkDuplicates.jar
	// INPUT=$input.bam
	// REMOVE_DUPLICATES=true
	// VALIDATION_STRINGENCY=LENIENT
	// AS=true
	// METRICS_FILE=${output}.mark_duplicates.log
	// OUTPUT=$output.bam
	// 2> /dev/null && md5sum $output.bam > ${output}.md5sum
	// """
	//
	// }

}

@filter("merge")

process mergeBams {

	input:
	path input

	output:
	path output

	script:
	"""
	java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/MergeSamFiles.jar
	${inputs.bam.split().collect { "INPUT="+it }.join(' ')}
	USE_THREADING=true
	VALIDATION_STRINGENCY=LENIENT
	AS=true
	OUTPUT=$output.bam
	"""

	// doc "Merge BAM files from multiple lanes or samples together. BAM files should have unique sample names and / or read groups"
	//
	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (branch.outdir) {
	// output.dir = branch.outdir + "/bam"
	// } else {
	// output.dir="bam"
	// }
	//
	// requires TMP : "Must specify a temporary directory (TMP)"
	// requires PICARD_HOME : "Must provide path to picard (PICARD_HOME)"
	//
	// exec """
	//
	// java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/MergeSamFiles.jar
	// ${inputs.bam.split().collect { "INPUT="+it }.join(' ')}
	// USE_THREADING=true
	// VALIDATION_STRINGENCY=LENIENT
	// AS=true
	// OUTPUT=$output.bam
	// """
}

@transform("bam")

process samToSortedBam {

	input:
	path input

	output:
	path output

	script:
	"""
	java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/SortSam.jar
	VALIDATION_STRINGENCY=LENIENT
	INPUT=$input.sam
	OUTPUT=$output.bam
	SORT_ORDER=coordinate
	"""

	// doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"
	//
	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (branch.outdir) {
	// output.dir = branch.outdir + "/bam"
	// } else {
	// output.dir="bam"
	// }
	//
	// requires TMP : "Must specify a temporary directory (TMP)"
	// requires PICARD_HOME : "Must provide path to picard (PICARD_HOME)"
	//
	// exec """
	// java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/SortSam.jar
	// VALIDATION_STRINGENCY=LENIENT
	// INPUT=$input.sam
	// OUTPUT=$output.bam
	// SORT_ORDER=coordinate
	// """
}

process repeatmodeler_format_genome {

	input:
	path input

	output:
	path output

	script:
	"""
	BuildDatabase -name $input.prefix -engine ncbi $input.fa
	"""

	// doc about: "A module to generate the input for RepeatModeler",
	// description: "Formats a genome sequence to be used by RepeatModeler",
	// author: "marc.hoeppner@bils.se"
	//
	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = branch.sample_dir }
	//
	// // requires here
	//
	// // Defining output directory
	// if (sample_dir) {
	// output.dir = branch.outdir + "/repeats"
	// } else {
	// output.dir = "repeats"
	// }
	//
	// transform("nhr") {
	// exec "BuildDatabase -name $input.prefix -engine ncbi $input.fa"
	// }
	//
	// forward input
}


process repeatmodeler_run {

	input:
	path input

	output:
	path output

	script:
	"""
	RepeatModeler -database $input.prefix -engine ncbi -pa $threads && cp RM_*/consensi.fa.classified $output 2>/dev/null
	"""

	// doc about: "A module to model repeats in a genome assembly",
	// description: "Runs RepatModeler on a genome assembly to predict repeats",
	// author: "marc.hoeppner@bils.se"
	//
	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = branch.sample_dir }
	//
	// // requires here
	//
	// // Defining output directory
	// if (sample_dir) {
	// output.dir = branch.outdir + "/repeats"
	// } else {
	// output.dir = "repeats"
	// }
	// uses(threads:16) {
	// exec "RepeatModeler -database $input.prefix -engine ncbi -pa $threads && cp RM_*/consensi.fa.classified $output 2>/dev/null"
	// }

}


repeatmodeler = segment {
	repeatmodeler_format_genome + repeatmodeler_run
}

process rseqc_read_quality {

	input:
	path input

	output:
	path output

	script:
	"""
	read_quality.py -i $input -o $target &>/dev/null
	"""

	// doc title: "RseQC - Analyzing read quality",
	//
	// desc: """
	// RseQC read_quality.py - analyzing and plotting
	// the quality of RNA-seq reads
	// """,
	//
	// constraints: "Requires a BAM file as input",
	//
	// author: "marc.hoeppner@bils.se"
	//
	// var paired : true
	// var sample_dir : false
	//
	// requires RSEQC_READ_QUALITY : "Must provide path to RseQC read_quality.py script"
	//
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/rseqc_read_quality"
	// } else {
	// output.dir = "rseqc_read_quality"
	// }
	//
	// // The output command line option is not identical to any of the
	// // produced files, so we need to work around that:
	// target = output.dir + "/" + input.prefix
	// product = target + ".qual.r"
	//
	// produce(product) {
	// exec "$RSEQC_READ_QUALITY -i $input -o $target &>/dev/null"
	// }
}

process rseqc_bam_stat {

	input:
	path input

	output:
	path output

	script:
	"""
	bam_stat.py --input-file=$input > $output.dir/rseqc_bam_stat.out 2> $output.dir/rseqc_bam_stat.err && md5sum $outputs >$output.dir/rseqc_bam_stat.md5
	"""

	// var sample_dir : false
	//
	// requires RSEQC_BAM_STAT : "Must provide path to RseQC bam_stat.py script"
	//
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/rseqc_bam_stat"
	// } else {
	// output.dir = "rseqc_bam_stat"
	// }
	//
	// input_extension = ".bam"
	//
	// products = [ "$output.dir/rseqc_bam_stat.out",
	// "$output.dir/rseqc_bam_stat.err" ]
	//
	// produce(products) {
	// exec "$RSEQC_BAM_STAT --input-file=$input >$output.dir/rseqc_bam_stat.out 2>$output.dir/rseqc_bam_stat.err && md5sum $outputs >$output.dir/rseqc_bam_stat.md5","rseqc"
	// }
}

process rseqc_read_distribution {

	input:
	path input

	output:
	path output

	script:
	"""
	rseqc --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED >$output.dir/rseqc_read_distribution.out 2>$output.dir/rseqc_read_distribution.err && md5sum $outputs >$output.dir/rseqc_read_distribution.md5
	"""

	// var sample_dir : false
	//
	// requires RSEQC_READ_DISTRIBUTION : "Must provide path to RseQC read_distribution.py script"
	// requires RSEQC_REF_GENE_MODEL_BED : "Must provide reference gene model in BED format"
	//
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/rseqc_read_distribution"
	// } else {
	// output.dir = "rseqc_read_distribution"
	// }
	//
	// input_extension = ".bam"
	//
	// products = [ "$output.dir/rseqc_read_distribution.out",
	// "$output.dir/rseqc_read_distribution.err" ]
	//
	// produce(products) {
	// exec "$RSEQC_READ_DISTRIBUTION --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED >$output.dir/rseqc_read_distribution.out 2>$output.dir/rseqc_read_distribution.err && md5sum $outputs >$output.dir/rseqc_read_distribution.md5","rseqc"
	// }
}

process rseqc_junction_annotation {

	input:
	path input

	output:
	path output

	script:
	"""
	rseqc --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED --out-prefix=$output.dir/$product_prefix >$output.dir/rseqc_junction_annotation.out 2>$output.dir/rseqc_junction_annotation.err && md5sum $outputs >$output.dir/rseqc_junction_annotation.md5
	"""

	// var sample_dir : false
	//
	// requires RSEQC_JUNCTION_ANNOTATION : "Must provide path to RseQC junction_annotation.py script"
	// requires RSEQC_REF_GENE_MODEL_BED : "Must provide reference gene model in BED format"
	//
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/rseqc_junction_annotation"
	// } else {
	// output.dir = "rseqc_junction_annotation"
	// }
	//
	// input_extension = ".bam"
	// product_prefix = "$output.dir/$input".replaceAll(/.*\//,"") - input_extension
	//
	// products = [ "$product_prefix" + ".junction_plot.r",
	// "$product_prefix" + ".junction.xls",
	// "$product_prefix" + ".splice_events.pdf",
	// "$product_prefix" + ".splice_junction.pdf" ]
	//
	// produce(products) {
	// exec "$RSEQC_JUNCTION_ANNOTATION --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED --out-prefix=$output.dir/$product_prefix >$output.dir/rseqc_junction_annotation.out 2>$output.dir/rseqc_junction_annotation.err && md5sum $outputs >$output.dir/rseqc_junction_annotation.md5","rseqc"
	// }
}

process sample_dir_prepare {

	input:
	path input

	output:
	path output

	script:
	"""
	"""

	// doc title: "A generic module to prepare the output options for a pipeline",
	//
	// desc:"""
	// This pipeline module can be used to set global variables required across all subsequent modules
	//
	// Options:
	//
	// sample_dir : false / true
	//
	// If set to true, all subsequent pipeline stages will be written into a subfolder based on the
	// name of the input file for a particular branch. The default is 'false'.
	//
	// For this to work, the module needs to check for branch.sample_dir.
	//
	// """,
	//
	// author: "mphoeppner@gmail.com"
	//
	// var sample_dir : false
	//
	// // Permanently store the sample name (to preserve it across forks)
	// branch.sample = branch.name
	// // Set the sample name as output directory
	// branch.outdir = branch.sample
	//
	// branch.sample_dir = sample_dir
	//
	// if (sample_dir) {
	// output.dir = branch.outdir
	// }

	// produce("sample.txt") {
	//	exec "echo $input > $output"
	// }

	// forward inputs
}

process indexBam {

	input:
	path input

	output:
	path output

	script:
	"""
	samtools index $input.bam
	"""

	// doc "A function to index a BAM file"
	//
	// requires SAMTOOLS : "Must provide path to samtools"
	//
	// transform("bam") to ("bam.bai") {
	// 	exec "$SAMTOOLS index $input.bam"
	// }
	// forward input
}

process flagstat {

	input:
	path input

	output:
	path output

	script:
	"""
	samtools flagstat $input.bam > $output
	"""

	// requires SAMTOOLS : "Must provide path to samtools"
	//
	// exec "$SAMTOOLS flagstat $input.bam > $output"
}

process samtools_filter_quality {

	input:
	path input

	output:
	path output

	script:
	"""
	samtools view -bq$quality -o $output $input && md5sum $output > ${output}.md5sum
	"""

	// doc "Filters BAM files by quality score"
	//
	// var sample_dir : false
	// var quality : "15"
	//
	// if(branch.sample_dir) { sample_dir = true }
	//
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/bam"
	// } else {
	// output.dir = "bam/${branch.sample}"
	// }
	//
	// requires SAMTOOLS : "Must provide path to samtools (SAMTOOLS)"
	//
	// transform("bam") to ("filtered.bam") {
	// exec "$SAMTOOLS view -bq$quality -o $output $input && md5sum $output > ${output}.md5sum"
	// }

}

process samtools_sort_bam {

	input:
	path input

	output:
	path output

	script:
	"""
	samtools sort -o $output $input && md5sum $output > ${output}.md5sum
	"""

	// doc "Sort BAM files. requires samtools version 1.2 or higher"
	//
	// requires SAMTOOLS : "Must provide path to samtools (SAMTOOLS)"
	//
	// // create output_name. Will be the input where we remove extension (should be .sam) and add .bam
	// String output_name = "$input" - ~/.[^.]*$/
	// String output = "${output_name}.sorted.bam"
	//
	// produce("$output") {
	// 	exec "$SAMTOOLS sort -o $output $input && md5sum $output > ${output}.md5sum"
	// }
}


process samtools_sam_to_bam {

	input:
	path input

	output:
	path output

	script:
	"""
	samtools view -bS -o $output $input && md5sum $output > ${output}.md5sum
	"""

	// doc "Convert SAM to BAM files."
	//
	// requires SAMTOOLS : "Must provide path to samtools (SAMTOOLS)"
	//
	// // create output_name. Will be the input where we remove extension (should be .sam) and add .bam
	// String output_name = "$input" - ~/.[^.]*$/
	// String output = "${output_name}.bam"
	//
	// produce("$output") {
	// 	exec "$SAMTOOLS view -bS -o $output $input && md5sum $output > ${output}.md5sum"
	// }
}
// Sequence conversion module

process cdna2protein {

	input:
	path input

	output:
	path output

	script:
	"""
	transeq -sequence $input -outseq $input"+".tmp -clean -trim
	sed s/_1//g $input"+".tmp > $output
	rm $input"+".tmp
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc "Converts cDNA sequences to protein"
	//
	// if (sample_dir) { output.dir = branch.name }
	//
	// filter("protein") {
	// exec "transeq -sequence $input -outseq $input"+".tmp -clean -trim"
	// exec "sed s/_1//g $input"+".tmp > $output"
	// exec "rm $input"+".tmp"
	// }
}


process stringtie {

	input:
	path input

	output:
	path output

	script:
	"""
	stringtie -l ${branch.sample} -o $output.dir"+"/transcripts.gtf -p $threads $options $input 2> $output"+".log
	"""

	// doc title: "Build transcripts from aligned reads using stringtie",
	// desc: """
	// Reconstructs transcripts from aligned RNA-seq reads in BAM format.
	//
	// Input:
	// A read alignment in BAM format
	//
	// Output:
	// transcripts.gtf
	//
	// Requires:
	// STRINGTIE : Path to the cufflinks binary
	//
	// Stage options:
	// sample_dir : true or false, determines whether output will be written to
	// a sample-specific subfolder (true) or not (false)
	// f : set the isoform fraction cut-off (default 0.1)
	// j : minimum junction coverage (default: 1)
	// G : set the reference annotation to use for guiding the assembly process (GTF/GFF3)
	// """,
	// author: "jacques.dainat@bils.se"
	//
	// var sample_dir : false
	// var stringtie_f : ""
	// var stringtie_j : ""
	// var stringtie_G : ""
	//
	// requires STRINGTIE : "Must provide path to stringtie (STRINGTIE)"
	//
	// def options = ""
	// def nameOut = "default"
	// // Checking which variables are set and build options string
	// if (stringtie_f) {
	// options += "-f $stringtie_f "
	// nameOut += "f" + stringtie_f
	// }
	// if (stringtie_j) {
	// options += "-j $stringtie_j "
	// nameOut += "j" + stringtie_j
	// }
	// if (stringtie_G) {
	// options += "-G $stringtie_G "
	// }
	//
	// branch.assembly_method = nameOut
	//
	// if (branch.sample_dir) { sample_dir = true }
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/stringtie/" + assembly_method
	// } else {
	// output.dir = "stringtie/" + branch.sample + "_" + assembly_method
	// }
	//
	// // The file to pass on is generally 'transcripts.gtf' - we use it as output.
	//
	// produce("transcripts.gtf") {
	// uses(threads:16) {
	// exec "$STRINGTIE -l ${branch.sample} -o $output.dir"+"/transcripts.gtf -p $threads $options $input 2> $output"+".log","stringtie"
	// }
	// }


}
// Tophat module

process tophat {

	input:
	path input

	output:
	path output

	script:
	"""
	if (PAIRED.toBoolean()) {
	exec "$TOPHAT $PHRED_SCORE $options -o $output.dir -p $threads --library-type=$LIBRARY_METHOD $BWT2_INDEX $input1 $input2 >$output.dir/tophat.out 2>$output.dir/tophat.err && md5sum $output >$output.dir/tophat.md5","tophat"
	} else {
	exec "$TOPHAT $PHRED_SCORE $options -o $output.dir -p $threads --library-type=$LIBRARY_METHOD $BWT2_INDEX $input >$output.dir/tophat.out 2>$output.dir/tophat.err && md5sum $output >$output.dir/tophat.md5","tophat"
	}
	"""

	// doc title: "Align RNA-seq reads against a reference with tophat",
	// desc: """
	// Uses Tophat2/Bowtie2 to align reads against a genome
	//
	// Stage options:
	// sample_dir : create sample-specific output folrder (boolean)
	//
	// Required variables are:
	// PAIRED : paired-end data (boolean)
	// TOPHAT : this should point to the location of the tophat executable
	// LIBRARY_METHOD : This specifies the type of sequencing library, e.g. fr-reverse
	// PHRED_SCORE : Specifies the quality encoding, e.g. usually this will be solexa-quals/phred33
	// BWT2_INDEX : The location of a bowtie2-formatted genome index
	// SAMTOOLS : This specifies the path to samtools
	// BOWTIE2 : This specifies the path to bowtie2
	// """,
	// constraints: """
	// Works with fq.gz and fq input files, but assumes paired-end
	// reads (if not, set paired to false).
	// The first argument is expected to be the left mate, the second
	// argument needs to be the right mate
	// """,
	// author: "marc.hoeppner@bils.se"
	//
	// // Exposed options
	// var sample_dir : false
	// var tophat_r : 50       // mate inner distance
	// var tophat_i : 50       // minimum intron length
	// var tophat_I : 500000   // maximum intron length
	// var GENOME_GTF : ""
	// var TRANSCRIPTOME_INDEX : ""
	// var tophat_T : false
	//
	// use_transcriptome = false
	//
	// options = "-r $tophat_r -i $tophat_i -I $tophat_I"
	//
	// // Check if an annotation file OR transcriptome index is passed and
	// // modify options
	// if (GENOME_GTF.length() > 0) {
	// options += " -G $GENOME_GTF"
	// use_transcriptome = true
	// } else if (TRANSCRIPTOME_INDEX.length() > 0) {
	// options += " --transcriptome-index $TRANSCRIPTOME_INDEX"
	// use_transcriptome = true
	// }
	//
	// // We enable quantifcation only against known transcripts but only
	// // if transcripts were provided
	// if (tophat_T && use_transcriptome) {
	// options += " -T"
	// }
	//
	// // Configuring the output directory
	// if (branch.sample_dir) { sample_dir = true }
	//
	// requires TOPHAT : "Must set the TOPHAT variable to point to tophat location"
	// requires LIBRARY_METHOD : "Must specify a sequencing library method (LIBRARY_METHOD)"
	// requires PHRED_SCORE : "Must set a phred score (PHRED_SCORE)"
	// requires BWT2_INDEX : "Must specify a Bowtie2 index (BWT2_INDEX)"
	// requires BOWTIE2 : "Must specify path to Bowtie2 (BOWTIE2)"
	// requires SAMTOOLS : "Must specify path to samtools (SAMTOOLS)"
	// requires PAIRED : "Must specify if the sample is stranded or not (true or false)"
	//
	// // We subsequently need to keep track of folders
	// // Here we set a name accessible to all subsequent modules.
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/tophat"
	// } else {
	// output.dir = "tophat/" + branch.sample
	// }
	//
	// // If a basename for this branch was set further upstream
	//
	// produce("accepted_hits.bam") {
	// uses(threads:16) {
	// if (PAIRED.toBoolean()) {
	// exec "$TOPHAT $PHRED_SCORE $options -o $output.dir -p $threads --library-type=$LIBRARY_METHOD $BWT2_INDEX $input1 $input2 >$output.dir/tophat.out 2>$output.dir/tophat.err && md5sum $output >$output.dir/tophat.md5","tophat"
	// } else {
	// exec "$TOPHAT $PHRED_SCORE $options -o $output.dir -p $threads --library-type=$LIBRARY_METHOD $BWT2_INDEX $input >$output.dir/tophat.out 2>$output.dir/tophat.err && md5sum $output >$output.dir/tophat.md5","tophat"
	// }
	// }
	// }
	//
	// check {
	// 	exec "[ -s $output ]"
	// } otherwise {
	// 	succeed "The Tophat output is empty. Stopping this branch ($branch.name)"
	// }
}
// Trimmomatic module

process trimmomatic {

	input:
	path input

	output:
	path output

	script:
	"""
	java -jar $TM_JAR PE -threads $threads $input1 $input2 ${output1} ${output3} ${output2} ${output4} ILLUMINACLIP:$TM_PATH/adapters/$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >$output.dir/trimmomatic.out 2>$output.dir/trimmomatic.err && md5sum $outputs >$output.dir/trimmomatic.md5
	"""

	// var sample_dir : false
	// if (branch.sample_dir) { sample_dir = true }
	//
	// doc title: "Adapter trimming of read files using Trimmomatic",
	// desc: """
	// Performs adapter trimming on paired-end RNA-seq reads.
	//
	// Requires:
	// TM_PATH : Location of the Trimmomatic files
	// TM_JAR : The name of the Trimmomatic jar file
	// ADAPTER : The name of the Trimmomatic adapter file
	// PAIRED : Bolean to know if sample are paired or not
	// """,
	// contraints: """
	// Files can be compressed (.fq.gz) or uncompressed (.fq)
	// """,
	// author: "mphoeppner@gmail.com"
	//
	// requires TM_PATH : "Must set TM_PATH variable to point to Trimmomatic folder"
	// requires TM_JAR : "Must set TM_JAR variable to point to Trimmomatic java file"
	// requires ADAPTER : "Must set the type of adapters to use"
	// requires PAIRED : "Bolean to know if sample are paired or not (true or false)"
	//
	// // Determine whether to write this into a sub-folder or not
	//
	// if (sample_dir) {
	// output.dir = branch.outdir + "/trimmomatic"
	// } else {
	// output.dir = "trimmomatic"
	// }
	//
	// input_extension = ".fq.gz"
	//
	// def products
	// def command
	//
	// if (PAIRED.toBoolean()) {
	// println "sample is paired"
	// products = [
	// ("$input1".replaceAll(/.*\//,"") - input_extension + '_paired.fq.gz'),
	// ("$input2".replaceAll(/.*\//,"") - input_extension + '_paired.fq.gz'),
	// ("$input1".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz'),
	// ("$input2".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz')
	// ]
	// } else {
	// println "sample is not paired"
	// products = [
	// ("$input".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz')
	// ]
	// }
	//
	// if (PAIRED.toBoolean()) {
	// produce(products) {
	// 	uses(threads:16) {
	// 		exec "java -jar $TM_JAR PE -threads $threads $input1 $input2 ${output1} ${output3} ${output2} ${output4} ILLUMINACLIP:$TM_PATH/adapters/$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >$output.dir/trimmomatic.out 2>$output.dir/trimmomatic.err && md5sum $outputs >$output.dir/trimmomatic.md5","trimmomatic"
	// 	}
	// }
	// } else {
	// 	produce(products) {
	// 		uses(threads:16) {
	// 			exec "java -jar $TM_JAR SE -threads $threads $input $output ILLUMINACLIP:$TM_PATH/adapters/$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >$output.dir/trimmomatic.out 2>$output.dir/trimmomatic.err && md5sum $output >$output.dir/trimmomatic.md5","trimmomatic"
	// 		}
	// 	}
	// }
	//
	// check {
	// 	exec "[ -s $output1 ]"
	// 	} otherwise {
	// 		succeed "The Trimmomatic filtering left not reads, stopping this branch ($branch.name)"
	// 	}
}

process verify_generic {

	input:
	path input

	output:
	path output

	script:
	"""
	"""

	// var binary : ""
	//
	// if("which $binary".execute().waitFor()!=0) {
	// 	fail "The path provided to $binary could not be resolved"
	// }
	//
	// forward inputs
}


// Verifications do not require an execute, so they will run on the head node and
// can be grouped as a segment (else it would mean submitting 4 or more jobs...)

verify_dependencies_rnaseq = segment {
	verify_generic.using(binary:"tophat") + verify_generic.using(binary:"cufflinks") + verify_generic.using(binary:"samtools") + verify_generic.using(binary:"bowtie2")
}

verify_annotation_preprocess = segment {
	verify_generic.using(binary:"fastasplit") + verify_generic.using(binary:"Rscript") + verify_generic.using(binary:"cpgplot") + verify_generic.using(binary:"bowtie2-build")
}

verify_dependencies_annotation_models = segment {
	verify_generic.using(binary:"gffread") + verify_generic.using(binary:"makeblastdb") + verify_generic.using(binary:"gff2gbSmallDNA.pl")
}
