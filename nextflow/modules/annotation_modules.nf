
process assembly_generate_stats {

	tag "Generating statistics for ${fasta_file.simpleName}"

	input:
	path fasta_file

	output:
	path "assembly_report.txt"

	script:
	"""
	fasta_statisticsAndPlot.pl --infile $fasta_file --output assembly_report.txt
	"""

}

process blastp {

	tag "Running blastp on ${query_fasta} against ${blastp_dbpath}"

	input:
	path query_fasta
	path blastp_dbpath
	val outfmt // 5
	val evalue

	output:
	path "results.blast"

	script:
	"""
	blastp -db $blastp_dbpath -query $query_fasta -outfmt $outfmt \
		-num_threads ${task.cpus} $evalue -out results.blast
	"""

}

process merge_blast_tab {

	tag "Merging Blast results"

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

	tag "Making Blast database: from: ${fasta_file} type: ${dbtype}"

	input:
	path fasta_file
	val dbtype // "prot"

	output:
	path output

	script:
	"""
	makeblastdb -in $fasta_file -dbtype $dbtype
	"""

}

process blast_recursive {

	tag "Performing recursive blast"

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

	tag "Creating Bowtie2 Index: ${species_id}"

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

	tag "Filtering Fasta ${fasta_file} by min length ${min_length}"

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
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""

}

process gbk2augustus {

	tag "Converting Genbank to Augustus"

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

	tag "Converting GFF to Genbank format"

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

	tag "Converting GFF to protein sequence"

	input:
	tuple sample_id, path(gff_file)
	tuple genome_id, path(genome)
	val codon_table

	output:
	path output

	script:
	"""
	tmp_fasta = $(mktemp -u )
	gff3_sp_extract_sequences.pl -o $tmp_fasta \
		-f $GENOME_FA -p -cfs -cis -ct $codon_table --gff $gff_file && \
		$BPIPE_BIN/fix_fasta.rb $tmp_fasta.fa > $output && \
		rm $tmp_fasta.fa
	"""

}

process gff_filter_by_blast {

	tag "Filtering GFF by Blast results"

	input:
	path gff_file
	path blast_file

	output:
	path output

	script:
	"""
	gff_filter_by_mrna_id.pl --gff $gff_file --blast $blast_file --outfile $output
	"""

	// doc about: "Remove GFF entries if they are found in a recursive blast (redundancy removal)",
	// description: "Takes a recursive blast input (format 6) and removes redundant entries from GFF file",
	// constraints: "Requires the blast file to be in outfmt 6",
}

process gff_filter_gene_models {

	tag "Filter gene models by GFF"

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

// Hisat2 module
process hisat2_index {
}

process hisat2 {

	input:
	tuple sample_id, path(reads)
	path hisat_index

	output:
	path output

	script:
	if (PAIRED.toBoolean()){
	"""
	hisat2 $options --novel-splicesite-outfile splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -1 $input1 -2 $input2 >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5
	"""
	} else {
	"""
	hisat2 $options --novel-splicesite-outfile splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -U $input >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5
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


process interpro {

	input:
	path protein_fasta
	path interpro_dbpath

	output:
	path '*.{gff3,tsv,xml}'

	script:
	"""
	interproscan $interpro_dbpath -i $protein_fasta -d $output -iprlookup -goterms -pa -dp
	"""

}


process merge_interpro_xml {

	input:
	path xml_files

	output:
	path 'interpro_xml'

	script:
	"""
	TMPFILE=$(mktemp)
	head -n 2 ${xml_files[0]} > \$TMPFILE
	for XML in $xml_files; do
		grep -v "xml" \$XML | grep -v "protein-matches" >> \$TMPFILE
	done
	tail -n 1 ${xml_files[0]} >> \$TMPFILE
	mv \$TMPFILE interpro_xml
	"""
}

process merge_interpro_tsv {

	input:
	path tsv_files

	output:
	path interpro_search.tsv

	script:
	"""
	cat $tsv_files >> interpro_search.tsv
	"""

}

process samtools_sort_bam {

	// Merge with process above

	input:
	path bam_alignment

	output:
	path 'sorted_alignment.bam'
	path 'sorted_alignment.bam.md5sum'

	script:
	"""
	samtools sort -o sorted_alignment.bam $bam_alignment && md5sum sorted_alignment.bam > sorted_alignment.bam.md5sum
	"""
}


process samtools_sam_to_bam {

	// Merge with process above

	input:
	path input

	output:
	path output

	script:
	"""
	samtools view -bS -o $output $input && md5sum $output > ${output}.md5sum
	"""

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
