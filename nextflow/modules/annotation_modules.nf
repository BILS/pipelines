process annie_interpro = {

  doc about: "Module to run Annie's InterPro converter",
    description: "Annie is a tool to convert search results into tabular format for functional annotation",
    constraints: "Requires the annie package to be available",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here
    requires ANNIE_ROOT : "Must provide path to Annie installation folder (ANNIE_ROOT)"

    // Defining output directory
    if (sample_dir) {
      output.dir = branch.outdir
    }

    exec "python $ANNIE_ROOT/annie.py ipr $input $output"

}

process assembly_generate_stats = {

  doc about: "Generates statistics from a genome assembly",
    description: "Accepts a multi-fasta file with nucleotide sequences and computes statistics",
    author: "marc.hoeppner@bils.se"

    // requires here

    var directory : "stats"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    produce("assembly_report.txt") {
    	exec "$BPIPE_BIN/fasta_statisticsAndPlot.pl --infile $input --output $output"
    }
}

process assembly_remove_mitochondria = {

  doc about: "Parses a genome assembly and removes organellar contigs",
    description: "Matches contigs against a database of organeller proteins to find orgenellar contaminations",
    constraints: "Requires a blast database of reliable organellar proteins",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here


    // Defining output directory
//    if (sample_dir) {
    //  output.dir = branch.outdir
//    } else {

//    }

    // exec "..."

   forward input

}

process gbk2training = {

	doc title: "Converts a GenBank file to an Augustus test and training set",

	desc: """
		Splits a Genbank file into two random sets for training and testing of
		Augustus profile models.

		Required variables:
		TEST_SIZE : How many models are to be kept for testing.
	""",
	constraints : "None",

	author : "marc.hoeppner@bils.se"

	var sample_dir : false

	requires TEST_SIZE : "Must provide a value for the size of the test data set (TEST_SIZE)"


	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir
	}

	produce(input+".train") {
		exec "randomSplit.pl $input $TEST_SIZE"
	}

}

process bayesembler = {

  doc about: "Module to run the Bayesembler transcriptome assembler",
    description: "Bayesembler uses a bayesian approach to reconstruct transcripts from BAM-formatted read alignments",
    constraints: "Requires a Tophat-generated BAM file and index",
    author: "marc.hoeppner@bils.se"

    var baysembler_c : "0.5"
    var bayesembler_m : "-m"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    produce("assembly.gtf") {
	    exec "bayesembler -b $input -p $threads $baysembler_m -c $bayesembler_c"
    }

}

process blast = {

        doc title: "Run a blast+ search on a FASTA file",

        desc: "Takes a sequence fasta file as input and blasts it against a database",

        author: "marc.hoeppner@bils.se"

        var sample_dir : false
        var outfmt : 5
        var program : "blastp"
        var max_target_seqs : 10

        requires BLASTP_DB : "Must specify a blast+ formatted protein database (BLASTP_DB)"
	requires BLAST_CUTOFF_EVALUE : "Must specify a blast cutoff evalue (BLAST_CUTOFF_EVALUE)"

        if (branch.sample_dir) { sample_dir = true}

        if (sample_dir) {
                output.dir = branch.outdir + "/theVoid"
        }

	var evalue : ""
	if (BLAST_CUTOFF_EVALUE){
		evalue = "-evalue $BLAST_CUTOFF_EVALUE"
	}

	uses(threads:8) {
               	exec "$program -db $BLASTP_DB -query $input -outfmt $outfmt -max_target_seqs $max_target_seqs -num_threads $threads $evalue -out $output","blast"
        }
}

process blastp = {

	doc title: "Run a blastp+ search on a protein FASTA file",

	desc: "Takes a protein fasta file as input and blasts it against a protein database",

	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	var outfmt : 5

	requires BLASTP_DB : "Must specify a blast+ formatted protein database (BLASTP_DB)"
        requires BLASTP : "Must provide path to blastp binary"
	requires BLAST_CUTOFF_EVALUE : "Must specify a blast cutoff evalue (BLAST_CUTOFF_EVALUE)"

	if (branch.sample_dir) { sample_dir = true}

	if (sample_dir) {
		output.dir = branch.outdir + "/theVoid"
	}

	var evalue : ""
        if (BLAST_CUTOFF_EVALUE){
		evalue = "-evalue $BLAST_CUTOFF_EVALUE"
	}

	uses(threads:8) {
	        exec "$BLASTP -db $BLASTP_DB -query $input -outfmt $outfmt -num_threads $threads $evalue -out $output","blastp"
	}
}

process recursive_blastp = {

	doc "Runs a blastp search against itself (blast+)"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) { output.dir = branch.outdir + "/theVoid" }

	requires BLASTP : "Must provide path to blastp binary"
	requires BLAST_CUTOFF_EVALUE : "Must specify a blast cutoff evalue (BLAST_CUTOFF_EVALUE)"

	var evalue : ""
        if (BLAST_CUTOFF_EVALUE){
                evalue = "-evalue $BLAST_CUTOFF_EVALUE"
        }

	exec "$BLASTP -db $input -query $input.fa -outfmt 5 -num_threads $threads $evalue -out $output"
}

process protein2blast_db = {

	doc "Creates a blast database from a protein fasta file (blast+)"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir + "/theVoid"
	}

	produce(input+".phr") {
		exec "makeblastdb -in $input.fa -dbtype prot"
	}

	forward input
}

recursive_protein2blastp = segment { protein2blast_db + recursive_blastp }

process merge_blast_xml = {

	doc "Crude method to merge the XML output from multiple BLAST searches"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir
	}

	produce(branch.sample + "_blast.out") {
		exec "cat $inputs > $output"
	}

}

process merge_blast_tab = {

	doc "Crude method to merge the tabular output from multiple BLAST searches"

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (sample_dir) {
                output.dir = branch.outdir
        }

        produce(branch.sample + "_blast.out") {
                exec "cat $inputs > $output"
        }

}

process run_blast2go = {

	doc title: "Runs the Blast2Go pipeline using blast and interpro XML inputs"

	var sample_dir : false

	if (branch.sample_dir) { sample_dir = true }

	if (branch.sample_dir) { output.dir = branch.outdir }

	// Assumes that the Blast output is input1 and interpro is input2 - BAD IDEA...

	produce(branch.sample+"_blast2go.annot") {
		exec "java -Xmx20G -cp ${B2G4PIPEPATH}/*:${B2G4PIPEPATH}/ext/* es.blast2go.prog.B2GAnnotPipe -in $input1 -out $output.prefix -prop ${B2G4PIPEPATH}/b2gPipe.properties -ips $input2 -annot"
	}

}


process blast2go2gff = {

	var sample_dir : false

	if (branch.sample_dir) { sample_dir = true }

	doc title: "Updates a GFF file with meta data from Blast2Go"

	requires B2GOGFF_SCRIPT : "Specify the full path to the b2gogff script"

	if (branch.sample_dir) { output.dir = branch.outdir }

	produce(gff_file_with_ids.prefix+".description.gff") {
		exec "$B2GOGFF_SCRIPT --gff $input.gff --b2go $input > $output"
	}

}

process blast_makeblastdb = {

  doc about: "Transforms a sequence file into a blast+ database",
    description: "Runs makeblastdb to generate a blast database",
    constraints: "Assumes protein input unless specified otherwise (dbtype)",
    author: "marc.hoeppner@bils.se"

    var dbtype : "prot"

    // requires here

    // Running a command
    produce(input + ".phr") {
	    exec "makeblastdb -in $input -dbtype $dbtype"
    }

    forward input

}

process blast_recursive = {

  doc about: "A model to blast a sequence file against a db of itself",
    description: "Takes a sequence file and blasts it against a database of the same name (recursive)",
    constraints: "Expects a blast database to exist for the input file and in the same location",
    author: "marc.hoeppner@bils.se"

    var directory : "blast_recursive"
    var blast_prog : "blastp"
    var blast_outfmt : "5"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    exec "$blast_prog -query $input -db $input -num_threads $threads -outfmt $blast_outfmt -out $output"

    branch.blastfile = output
}

process bowtie2_index = {

  doc about: "A pipeline to generate a bowtie2 index from a genome sequence",
    description: "Takes a genome sequence in fasta format and generates a bowtie2 index",
    author: "marc.hoeppner@bils.se"

    var directory : "bowtie-index"

    if (branch.sample_dir) {
	output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
	if (directory.length() > 0) {
		output.dir = directory
	}
    }

    // requires here

    produce(input.prefix + ".1.bt2") {
    	exec "bowtie2-build $input $output.prefix.prefix >/dev/null"
    }

    forward input

}
@transform("sai")
alignBWA = {
        doc "Aligns using BWA. Note: assumes input file are gzipped"

        var sample_dir : false
	// Exposed variables with defaults
	var bwa_l : 32 // seed length
	var bwa_k : 2 // max differences in seed

	options = "-l $bwa_l -k $bwa_k"

        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

        requires BWA : "Must provide path to bwa (BWA)"
        requires REF : "Must provide PATH to reference (REF)"

        exec "$BWA aln $options -t $threads $BWA_INDEX $input.gz > $output.sai"
}

@transform("sam")

process alignToSamPE = {
    	doc "Create SAM files from BWA alignment. Note that these will be very large."

       	var sample_dir : false

	if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires BWA : "Must provide path to bwa (BWA)"
	requires PLATFORM : "Must specify a sequencing platform (PLATFORM)"

    	branch.lane = (input.sai =~ /.*L([0-9]*)_*R.*/)[0][1].toInteger()
    	branch.sample = branch.name

	exec """
        	$BWA sampe $REF -r "@RG\\tID:1\\tPL:$PLATFORM\\tPU:${branch.lane}\\tSM:${branch.sample}"  $input1.sai $input2.sai $input2.gz $input2.gz > $output.sam
    	"""
}

process cds2protein = {

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        doc "Converts cDNA sequences to protein"

        if (sample_dir) { output.dir = branch.name }

        filter("protein") {
		exec "transeq -sequence $input -outseq $output -clean -trim"
		//exec "$BPIPE_BIN/cds2protein.pl --infile $input --outfile $output 2>/dev/null"
        }
}

process cegma = {

  doc about: "A wrapper around the Cegma package",
    description: "Cegma checks the expected gene space coverage of an assembly",
    constraints: "Requires the Cegma package",
    author: "marc.hoeppner@bils.se"

    var directory : "cegma"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    produce("cegma.completeness_report") {
	uses(threads:16) {
	    	exec "cegma -T $threads -g $input -o $output.dir/cegma"
	}
    }

}

process cmsearch = {

  doc about: "Module to run Infernal's cmsearch algorithm",
    description: "Cmsearch identifies putative non-coding RNAs in nucleotide sequences",
    author: "marc.hoeppner@bils.se"

    var threads : 1
    var db : "/projects/references/databases/rfam/11.0/models_1_1/E_plus.cm"
    var directory : "rfam"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here


    exec "cmsearch --cpu $threads --rfam --cut_tc --tblout $output $db $input >/dev/null"
}

process cufflinks = {

	 doc title: "Build transcripts from aligned reads using cufflinks",
                desc: """
                        Reconstructs transcripts from aligned RNA-seq reads in BAM format.

                        Input:
                        A read alignment in BAM format

                        Output:
                        transcripts.gtf

                        Requires:
                        GENOME_FA : Genome sequence in FASTA format
			CUFFLINKS : Path to the cufflinks binary

                        Stage options:
                        sample_dir : true or false, determines whether output will be written to
                                a sample-specific subfolder (true) or not (false)
                        j : set the isoform fraction cut-off
                        F : set the intra-splice-junction read cut-off
                """,
	        author: "marc.hoeppner@bils.se"

	var sample_dir : false
        var cufflinks_j : "0.15"
        var cufflinks_F : "0.10"
	var cufflinks_g : false
	var cufflinks_min_intron_length : "50"
        var cufflinks_I : "300000"
	var GENOME_GTF : ""
	var LIBRARY_METHOD : "fr-unstranded"

	requires CUFFLINKS : "Must provide path to cufflinks (CUFFLINKS)"
        requires GENOME_FA : "Please set the GENOME_FA variable"

	def options

	// Checking which variables are set and build options string
	if (cufflinks_g && GENOME_GTF.length() > 0) {
		options = "-j $cufflinks_j -F $cufflinks_F -g"
	} else {
		options = "-j $cufflinks_j -F $cufflinks_F"
	}

	// Check if an annotation file is passed and modify options if so
	if (GENOME_GTF.length() > 0) {
  		options += " -G $GENOME_GTF"
	}

        branch.assembly_method = "j" + cufflinks_j + "_F" + cufflinks_F

	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir + "/cufflinks/" + assembly_method
        } else {
                output.dir = "cufflinks/" + branch.sample + "_" + assembly_method
        }

        // The file to pass on is generally 'transcripts.gtf' - we use it as output.

	produce("transcripts.gtf") {
		uses(threads:16) {
		        exec "$CUFFLINKS --library-type=$LIBRARY_METHOD -L ${branch.sample} -o $output.dir -p $threads -u -b $GENOME_FA $options $input 2> $output"+".log","cufflinks"
		}
	}


}

process emboss_cpgplot = {

  doc about: "A module to run Emboss' cpgplot tool",
    description: "cpgplot predicts cpg islands in genomic sequences",
    constraints: "Should only be used for vertebrate genomes",
    author: "marc.hoeppner@bils.se"

    var cpgisland_window : 100
    var cpgisland_minlen : 200
    var cpgisland_minoe : 0.6
    var cpgisland_minpc : 50.0

    var directory : "cpgplot"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    produce($input+".cpg.gff") {
    	exec "cpgplot -sequence $input -outfeat $output -outfile $output" + ".out -window $cpgisland_window -minlen $cpgisland_minlen -minoe $cpgisland_minoe -minpc $cpgisland_minpc -noplot -nocg"
    }
}

process fasta_explode = {

  doc about: "Explodes a fasta file into its individual sequences",
    description: "Creates one file per sequence in the input FASTA file",
    constraints: "Requires exonerate to be loaded",
    author: "marc.hoeppner@bils.se"

    var directory : "sequences"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    // Doesn't produce a pre-definable output, so we use a log file as dummy target
    produce("fastaexplode.log") {
	    exec "fastaexplode -d $output.dir -f $input 2> $output"
    }

    forward input
}

process fasta_filter_size = {

  doc about: "Parses a FASTA file and removes sequences smaller than the cutoff",
    description: "Filters fasta file by size",
    constraints: "Requires bioruby to be installed",
    author: "marc.hoeppner@bils.se"

    var size : 1000
    var directory : ""

    // requires here

    // Defining output directory
    if (branch.sample_dir) {
	output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
      if (directory.length() > 0) {
	      output.dir = directory
      }
    }

    filter("filtered") {
	    exec "$BPIPE_BIN/fasta_filter_size.rb -i $input -s $size -o $output"
    }

}
// File operation module

process fastasplit = {

	doc title: "Splits a FASTA file into chunks, using the CHUNKS variable",
		desc: """
			Splits a fasta-formatted file into chunks of similar
			size using the exonerate tool fastasplit.

			Input:
			A fasta-formatted sequence file

			Requires:
			CHUNKS : defines the number of chunks to split into.

		""",

		constraints: """
			The number of chunks should be chosen so that at least
			10 sequences are included in each chunk - else the splitting
			may fail.
		""",
		author: "marc.hoeppner@bils.se"

        requires CHUNKS : "Must set variable CHUNKS"

	if (branch.sample_dir) { output.dir = branch.outdir }

	def chunkfiles = []

	for ( value in (0..(CHUNKS.toInteger()-1)) ) {
		if (value < 10) {
			chunkfiles.push(input+"_chunk_000000${value}")
		} else if (value < 100) {
			chunkfiles.push(input+"_chunk_00000${value}")
		} else if (value < 1000) {
			chunkfiles.push(input+"_chunk_0000${value}")
		} else {
			chunkfiles.push(input+"_chunk_000${value}")
		}
	}

      	//  def chunkfiles = (0..(CHUNKS.toInteger()-1)).collect{ input+"_chunk_000000${it}" }

	produce(chunkfiles) {
                exec "fastasplit -f $input -o ${output.dir} -c $CHUNKS"
        }

}

process fastqc = {

	doc title: "Quality control of read data using FastQC",

	desc: """
		FastQC is a light-weight java tool that analyses RNAseq
		read data and reports statistics on the sequencing quality.
	""",

	constraints: "Requires one (single-end) or two (paired-end) gzipped fastq files",

	author: "marc.hoeppner@bils.se"

	var sample_dir : false // Write output to a sample-specific directory
	var paired : true // input data is paired

	requires FASTQC : "Must provide path to fastqc (FASTQC)"

	input_extension = ".gz"

	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir + "/fastqc"
	} else {
		output.dir = "fastqc"
	}

	def products

	if (paired) {
		products = [
			("$input1".replaceAll(/.*\//,"") - input_extension + '_fastqc.html'),
			("$input2".replaceAll(/.*\//,"") - input_extension + '_fastqc.html')
		]
	} else {
		products = [
			("$input".replaceAll(/.*\//,"") - input_extension + '_fastqc.html')
		]
	}

	if (paired) {
		produce(products) {
			multi "fastqc --outdir=${output.dir} $input1","fastqc --outdir=${output.dir} $input2"
		}
	} else {
		produce(products) {
			exec "fastqc --outdir=${output.dir} $input"
		}
	}
}

process realignIntervals = {

	// Hard-coded to take 2 known indels files right now

	doc "Identify realign intervals with GATK"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
    	requires REF : "Must provide a path to reference (REF)"
	requires GOLD_STANDARD_INDELS : "Must provide a path to gold standard indels (GOLD_STANDARD_INDELS)"
	requires INDELS_100G : "Must provide a path to indels 100G (INDELS_100G)"
	requires LOG : "Must provide a location for log file (LOG)"

	exec """
        	java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I $input.bam --known $GOLD_STANDARD_INDELS --known $INDELS_100G -log $LOG -o $output.intervals
   	 """
}

process realign = {

	doc "Realign indels with GATK"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"

    	exec """
        	java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $input.bam -targetIntervals $input.intervals -log $LOG -o $output.bam
    	"""
}


process baseQualRecalCount = {

    	doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"
	requires DBSNP : "Must provide a path to known SNP DB (DBSNP)"

	exec "java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $input.bam -R $REF --knownSites $DBSNP -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -log $LOG -o $output.counts"
}

process baseQualRecalTabulate = {

	doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"

	exec "java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -I $input.bam -BQSR $input.counts -R $REF -l INFO -log $LOG -o $output"
}

process callSNPs = {
    	doc "Call SNPs/SNVs using GATK Unified Genotyper"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"
	requires DBSNP : "Must provide path to refence SNP DB (DBSNP)"

 	exec """
            java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper
               -nt $threads
               -R $REF
               -I $input.bam
               --dbsnp $DBSNP
               -stand_call_conf 50.0 -stand_emit_conf 10.0
               -dcov 1600
               -l INFO
               -A AlleleBalance -A DepthOfCoverage -A FisherStrand
               -glm SNP -log $LOG
               -o $output.vcf
        """
}

process callIndels = {

	doc "Call variants using GATK Unified Genotyper"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"
	requires DBSNP : "Must provide path to refence SNP DB (DBSNP)"

    	exec """
        	java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper
             	-nt $threads
             	-R $REF
             	-I $input.bam
             	--dbsnp $DBSNP
             	-stand_call_conf 50.0 -stand_emit_conf 10.0
             	-dcov 1600
             	-l INFO
             	-A AlleleBalance -A DepthOfCoverage -A FisherStrand
             	-glm INDEL
             	-log $LOG -o $output.vcf
    	"""
}

@filter("filter")

process filterSNPs = {
    	// Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires LOG : "Must provide a location for log file (LOG)"

    	exec """
        	java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration
             	-R $REF
             	--filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
             	--filterName 'GATK_MINIMAL_FILTER'
             	--variant $input.vcf
             	-log $LOG
             	-o $output.vcf
    	"""
}

@filter("filter")

process filterIndels = {
    	doc """
            Filter data using very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
            If you have 10 or more samples GATK also recommends the filter InbreedingCoeff < -0.8
        	"""
	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"

    	exec """
        	java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration
                    -R $REF
                    --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0'
                    --filterName 'GATK_MINIMAL_FILTER' -log $LOG
                    --variant $input.vcf
                    -o $output.vcf
    	"""
}

process gbk2augustus = {

    doc about: "A generic module that needs a description",
    description: "Description here",
    constraints: "Information on constraints here",
    author: "marc.hoeppner@bils.se"

    var test_size : 100

    // requires here

    // Running a command
    produce(input+".train") {
	    exec "randomSplit.pl $input $test_size"
    }

}

process genome_sanity_check = {

  doc about: "A module to analyse and produce statistics for a genome assembly",
    description: "Verifies the integrity of a genome assembly for annotation",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // Defining output directory
    if (sample_dir) {
      output.dir = branch.outdir
    }

   // exec "..."

   forward input
}

process genome_tools_gff_sort = {

	doc about: "A module to sort a GFF3-formatted annotation by coordinates using GenomeTools",

	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires GENOME_TOOLS : "Must provide path to genometools (GENOMETOOLS)"

	if (sample_dir) { output.dir = branch.outdir }

	filter("sorted") {
		exec "$GENOME_TOOLS gff3 -sort $input > $output 2>/dev/null"
	}

}

process genome_tools_gff_stats = {

  doc title: "Generate statistics from a gff3 file using the genometools package",

  desc: "Genometools compiles feature counts from a gff3-formatted annotation file",

  author: "marc.hoeppner@bils.se"

  var sample_dir : false
  if (branch.sample_dir) { sample_dir = true }

  requires GENOME_TOOLS : "Must provide path to genometools (GENOME_TOOLS)"

  if (sample_dir) {
    output.dir = branch.outdir
  }

  transform(".gff") to (".gtcounts") {
    exec "gt stat $input > $output 2>/dev/null"
  }

}

process genome_tools_gff_to_gtf = {

	doc about: "A module to convert a GFF3 file to GTF using genometools",
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires GENOME_TOOLS : "Must provide path to genometools (GENOME_TOOLS)"

	if (sample_dir) { output.dir = branch.outdir }

	transform(".gff") to (".gtf") {
		exec "$GENOME_TOOLS gff3_to_gtf -o $output -force $input 2>/dev/null"
	}
}

process gff2cds = {

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        doc "Extracts cDNA sequences from the annotation/genome"

        requires GENOME_FA : "Must set variable GENOME_FA"

        if (sample_dir) {
                output.dir = branch.outdir
        }

        transform(".gff") to (".cds.fa") {
                exec "gffread -x $output -g $GENOME_FA $input"
        }
}

process gff2gbk = {

    doc "Converts a GFF3 formatted file into GenBank format"

    var directory : "gff2genbank"
    var flank : 1000

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here
    requires GENOME_FA : "Must provide a genome sequence in fasta format (GENOME_FA)"

    // Running a command

    transform("gbk") {
        exec "gff2gbSmallDNA.pl $input $GENOME_FA $flank $output"
    }
}

process gff2protein = {

    doc about: "A module to extract protein sequences from GFF annotations",
    description: "Reports protein sequences from GFF annotations",
    author: "marc.hoeppner@bils.se, jacques.dainat@bils.se"

    var directory : "protein"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here
    requires GENOME_FA : "Must provide a genome sequence (GENOME_FA)"
    requires CODON_TABLE :  "Must specify the translation table (CODON_TABLE)"

    // Running a command
    transform(".gff") to (".proteins.fa") {
            exec "$BPIPE_BIN/gff3_sp_extract_sequences.pl -o $input.prefix"+".tmp -f $GENOME_FA -p -cfs -cis -ct $CODON_TABLE --gff $input && $BPIPE_BIN/fix_fasta.rb $input.prefix"+".tmp.fa > $output && rm $input.prefix"+".tmp.fa"
    }

}

process gff_annotation_stats = {

	doc about: "A module to generate annotation statistics, based on the code by A. Kähäri",
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires KAHARI_ANNOTATION_STATS : "Must provide path to annotation_stats.pl (KAHARI_ANNOTATION_STATS)"
	requires GENOME_FA : "Must provide path to genome sequence in FASTA format (GENOME_FA)"

	if (sample_dir) { output.dir = branch.outdir }

	transform(".gff") to (".stats") {
		exec "$KAHARI_ANNOTATION_STATS --gff $input --fasta $GENOME_FA > $output"
	}


}

process gff_filter_by_blast = {

  doc about: "Remove GFF entries if they are found in a recursive blast (redundancy removal)",
    description: "Takes a recursive blast input (format 6) and removes redundant entries from GFF file",
    constraints: "Requires the blast file to be in outfmt 6",
    author: "marc.hoeppner@bils.se"

    var directory : "nonredundant"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    blast_file = input
    // input will change due to the expression below, so we capture it here

    produce(branch.name + ".nr.gff") {
	from("*.gff") {
                exec "perl $BPIPE_BIN/gff_filter_by_mrna_id.pl --gff $input --blast $blast_file --outfile $output"
	}
    }
}

process gff_filter_gene_models = {

  doc about: "Wrapper module around A. Kaharis' gff filter_sort.pl script",
    description: "Filters transcript models by a set of criteria",
    constraints: "Requires the Kahari code base to be present",
    author: "marc.hoeppner@bils.se"

    var directory : "filter"

    var filter_c : "-c"
    var filter_a : "0.3"
    var filter_r : "-r"
    var filter_d : "500"

    def options = filter_c + " " + filter_r + " -d " + filter_d + " -a " + filter_a

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here
    requires KAHARI_GIT : "Must set path to Kahari git directory (KAHARI_GIT)"
    requires GENOME_FA : "Must provide a genome sequence in FASTA format"

    // Running a command

    filter("filter") {
	    exec "perl $KAHARI_GIT/scripts/GFF/filter_sort.pl -f $input -F $GENOME_FA -o $output $options"
    }
}

process gff_get_trna = {

	doc about: "Simple operation to excise tRNAs from a gff annotation",
	constraints: "Expects a the annotation to come from Maker and contain mRNAs and tRNAs, nothing else",
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) { output.dir = branch.outdir}

	transform(".gff") to (".mRNA.gff") {
		exec "grep \"trnascan\" -v $input.gff > $output"
	}

}

process gff_longest_cds = {

  doc about: "A module to limit gene loci in a GFF3-formatted annotation to the longest product",
    description: "Parses a GFF3-formatted annotation and filters out all transcripts except those with the longest CDS",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    requires KAHARI_GIT : "Must set path to Kahari git directory (KAHARI_GIT)"


    // Defining output directory
    if (sample_dir) {
      output.dir = branch.outdir
    }

    transform(".gff") to (".longest_cds.gff") {
    	exec "$KAHARI_GIT/scripts/GFF/find_longest_CDS.pl -f $input -o $output"
    }

}

process gff_stable_id = {

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        doc title: "Generates stable IDs for all features in the GFF file"

        requires ID_TRUNK : "Must set variable ID_TRUNK"

        if (sample_dir) {
                output.dir = branch.outdir
        }

        filter("stable_id") {
                exec "$BPIPE_BIN/gff_create_stable_id.pl --gff $input --id_trunk $ID_TRUNK --outfile $output"
        }

        branch.gff_file_with_ids = output

}

process gffread_extract_sequences = {

	doc about: "A module to extract CDS and protein sequences from an annotation using gffread",
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires GENOME_FA : "Must provide genome sequence in FASTA format (GENOME_FA)"

	if (sample_dir) { output.dir = branch.outdir }

	produce(input.prefix+".cds.fa",input.prefix+".proteins.fa") {
		exec """
			gffread -x $output1 -g $GENOME_FA $input ; gffread -y $output2 -g $GENOME_FA $input
		"""
	}
}
// Hisat2 module

process hisat2 = {

    doc title: "Align RNA-seq reads against a reference with hisat2",
        desc: """
            Uses hisat2/Bowtie2 to align reads against a genome

            Stage options:
                sample_dir : create sample-specific output folrder (boolean)

            Required variables are:
                PAIRED : paired-end data (boolean)
		HISAT2 : this should point to the location of the tophat executable
                LIBRARY_METHOD : This specifies the type of sequencing library: fr/rf/ff
                PHRED_SCORE : Specifies the quality encoding, e.g. usually this will be solexa-quals/phred33/phred64
                HISAT2_INDEX : The location of a hisat2-formatted genome index
               // SAMTOOLS : This specifies the path to samtools
              //  BOWTIE2 : This specifies the path to bowtie2
        """,
        constraints: """
            Works with fq.gz and fq input files, but assumes paired-end
            reads (if not, set paired to false).
            The first argument is expected to be the left mate, the second
            argument needs to be the right mate
        """,
        author: "jacques.dainat@nbis.se"

    // Exposed options
    var sample_dir : false

    options = ""
    // Handle --rna-strandness option
    LIBRARY_METHOD=LIBRARY_METHOD.toLowerCase()
    if (LIBRARY_METHOD == "r" || LIBRARY_METHOD == "rf" || LIBRARY_METHOD == "fr-firstrand" || LIBRARY_METHOD == "fr-firststrand") {
        if (PAIRED.toBoolean()) {
		options += "--rna-strandness RF"
	}
	else{ // single end reads
		options += "--rna-strandness R"
	}
    }
    else if (LIBRARY_METHOD == "f" || LIBRARY_METHOD == "fr" || LIBRARY_METHOD == "fr-secondstrand" ) {
        if (PAIRED.toBoolean()) {
		options += "--rna-strandness FR"
	}
	else{ // single end reads
		 options += "--rna-strandness F"
	}
    }
    else{ // default value is unstranded
	    println "LIBRARY_METHOD is unstranded"
    }

    // Handle Quality score
    if (PHRED_SCORE == "phred33") {
        options += " --phred33"
    }
    else if (PHRED_SCORE == "phred64") {
        options += " --phred64"
    }
    else if (PHRED_SCORE == "solexa-quals") {
        options += " --solexa-quals"
    }

    // Configuring the output directory
    if (branch.sample_dir) { sample_dir = true }

    requires HISAT2 : "Must set the TOPHAT variable to point to tophat location"
    requires LIBRARY_METHOD : "Must specify a sequencing library method (LIBRARY_METHOD)"
    requires PHRED_SCORE : "Must set a phred score (PHRED_SCORE)"
    requires HISAT2_INDEX : "Must specify a Bowtie2 index (BWT2_INDEX)"
    requires PAIRED : "Must specify if the sample is stranded or not (true or false)"

    // We subsequently need to keep track of folders
    // Here we set a name accessible to all subsequent modules.
    if (sample_dir) {
        output.dir = branch.outdir + "/hisat2"
    } else {
        output.dir = "hisat2/" + branch.sample
    }

    // If a basename for this branch was set further upstream

    produce("accepted_hits.sam") {
        uses(threads:16) {
            if (PAIRED.toBoolean()) {
                exec "$HISAT2 $options --novel-splicesite-outfile $output.dir/splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -1 $input1 -2 $input2 >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5","hisat2"
            } else {
                exec "$HISAT2 $options --novel-splicesite-outfile $output.dir/splicesite.txt -S $output.dir/accepted_hits.sam -p $threads -x $HISAT2_INDEX -U $input >$output.dir/hisat2.out 2>$output.dir/hisat2.err && md5sum $output >$output.dir/hisat2.md5","hisat2"
            }
        }
    }

    check {
        exec "[ -s $output ]"
    } otherwise {
        succeed "The HISAT2 output ($output) is empty. Stopping this branch ($branch.name)"
    }
}

process htseq_count = {

	doc title: "HTseq-count counts reads overlapping an annotation",

	desc: """
		HTseq-count is a python tool designed to accept a read alignment in BAM or SAM
		format and counts the reads overlapping annotated exons provided in GTF format.
	""",

	constraints: "Requires a read alignment in BAM/SAM format and GTF-formatted annotation",

	author: "marc.hoeppner@bils.se"

	var sample_dir : true

	requires GENOME_GTF : "Must provide an annotation in GTF format (GENOME_GTF)"
	requires HTSEQ_COUNT : "Must provide path to htseq-count (HTSEQ_COUNT)"

	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir + "/htseq_count"
	} else {
		output.dir = "htseq_count"
	}

	transform(".bam") to (".htseqcount") {
		exec "HTSEQ_COUNT -f bam -t exon $input.bam $GENOME_GTF > $output"
	}

}

process igvcount = {
    exec "igvtools count $input.bam $output hg19"
}
// Interproscan Module


process interpro = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Runs a protein fasta file against InterPro"

	if (sample_dir) { output.dir = branch.outdir }

	requires INTERPROSCAN : "Must specify the location of interproscan (INTERPROSCAN)"
	requires INTERPRO_DB_LIST : "Must specify the list of DB to use (INTERPRO_DB_LIST)"

	// set dblist if needed - if db_list is empty, all the DB will be used
	var db_list : ""
	if (INTERPRO_DB_LIST.toLowerCase() != "all"){
		db_list = "-appl $INTERPRO_DB_LIST"
	}

	produce(input+".gff3",input+".tsv",input+".xml") {
		exec "$INTERPROSCAN $db_list -i $input -d ${output.dir} -iprlookup -goterms -pa -dp > /dev/null 2> /dev/null ","interpro"
        }
}


process merge_interpro_xml = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Merges the XML output from multiple InterPro searches"

	if (sample_dir) { output.dir = branch.outdir }

	produce(branch.sample+"_interpro.xml") {

		// This is a really stupid hack to merge XML files
		// We take the first lines and the last line from the first file
		// And squeeze the other stuff from all files in between

		def first_file = inputs[0]
		def temp_file = branch.sample + "_interpro.tmp"

		exec "head -n 2 $first_file > $temp_file"

		for (i in inputs) {
			exec "grep -v xml $i | grep -v protein-matches >> $temp_file"
		}

		exec "tail -n 1 $first_file >> $temp_file"

		exec "mv $temp_file $output"
	}

	branch.ipr_xml = output
}

process merge_interpro_tsv = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Merges the TSV output from multiple InterPro searches"

	if (sample_dir) { output.dir = branch.outdir }

	produce(branch.sample+"_interpro.tsv") {

		exec "cat $inputs >> $output"

	}

	// We save the output name since we need it for a later pipeline stage.

	branch.iprtsv = output

	// Pass on the name of the output file


}


process interpro2gff = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Converts InterPro results into GFF format"

	if (sample_dir) { output.dir = branch.outdir }


	// This is a horrible way of doing this...
	// We need the name of a file set in an earlier stage,
	// but the variable doesn't carry over, so:

	if (sample_dir) {
		iprtsv = branch.outdir + "/" + branch.sample + "_interpro.tsv"
	} else {
		iprtsv = branch.sample + "_interpro.tsv"
	}

	filter("interpro") {
		exec "ipr_update_gff $input $iprtsv > $output"
	}
}
// Maker module


process cufflinks2maker = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Converts cufflinks-reconstructed transcripts to Maker-compatible GFF3 format"

	if (sample_dir) {
		output.dir = branch.outdir + "/maker_tracks"
	} else {
		output.dir = "maker_tracks"
	}

	transform("gtf") to (assembly_method + ".maker.gff3") {
		exec "cufflinks2gff3 $input > $output"
	}

}

process dedup = {

        doc "Marks duplicate reads in BAM file using Picard"

        var sample_dir : false

        if (branch.sample_dir) { sample_dir = true }

        if (sample_dir) {
                output.dir = branch.outdir + "/bam"
        } else {
                output.dir="bam/${branch.sample}"
        }

        requires TMP : "Must specify a temporary directory (TMP)"
        requires PICARD_HOME : "Must provide path to picard (PICARD_HOME)"

        transform("bam") to ("dedup.bam") {

                exec """
                        java -Xmx6g -Djava.io.tmpdir=$TMP -jar $PICARD_HOME/MarkDuplicates.jar
                        INPUT=$input.bam
                        REMOVE_DUPLICATES=true
                        VALIDATION_STRINGENCY=LENIENT
                        AS=true
                        METRICS_FILE=${output}.mark_duplicates.log
                        OUTPUT=$output.bam
                        2> /dev/null && md5sum $output.bam > ${output}.md5sum
                 """

        }

}

@filter("merge")

process mergeBams = {

        doc "Merge BAM files from multiple lanes or samples together. BAM files should have unique sample names and / or read groups"

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/bam"
        } else {
                output.dir="bam"
        }

        requires TMP : "Must specify a temporary directory (TMP)"
        requires PICARD_HOME : "Must provide path to picard (PICARD_HOME)"

        exec """

                java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/MergeSamFiles.jar
                ${inputs.bam.split().collect { "INPUT="+it }.join(' ')}
                USE_THREADING=true
                VALIDATION_STRINGENCY=LENIENT
                AS=true
                OUTPUT=$output.bam
        """
}
@transform("bam")

process samToSortedBam = {

	doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/bam"
        } else {
                output.dir="bam"
        }

	requires TMP : "Must specify a temporary directory (TMP)"
	requires PICARD_HOME : "Must provide path to picard (PICARD_HOME)"

	exec """
	        java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/SortSam.jar
                    VALIDATION_STRINGENCY=LENIENT
                    INPUT=$input.sam
                    OUTPUT=$output.bam
                    SORT_ORDER=coordinate
    	"""
}

process repeatmodeler_format_genome = {

  doc about: "A module to generate the input for RepeatModeler",
    description: "Formats a genome sequence to be used by RepeatModeler",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here

    // Defining output directory
    if (sample_dir) {
      output.dir = branch.outdir + "/repeats"
    } else {
      output.dir = "repeats"
    }

    transform("nhr") {
    	exec "BuildDatabase -name $input.prefix -engine ncbi $input.fa"
    }

    forward input
}


process repeatmodeler_run = {

    doc about: "A module to model repeats in a genome assembly",
    description: "Runs RepatModeler on a genome assembly to predict repeats",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here

    // Defining output directory
    if (sample_dir) {
      output.dir = branch.outdir + "/repeats"
    } else {
      output.dir = "repeats"
    }
    uses(threads:16) {
    	exec "RepeatModeler -database $input.prefix -engine ncbi -pa $threads && cp RM_*/consensi.fa.classified $output 2>/dev/null"
    }

}


repeatmodeler = segment {

        repeatmodeler_format_genome + repeatmodeler_run

}

process rseqc_read_quality = {
    doc title: "RseQC - Analyzing read quality",

    desc: """
        RseQC read_quality.py - analyzing and plotting
        the quality of RNA-seq reads
    """,

    constraints: "Requires a BAM file as input",

    author: "marc.hoeppner@bils.se"

    var paired : true
    var sample_dir : false

    requires RSEQC_READ_QUALITY : "Must provide path to RseQC read_quality.py script"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_read_quality"
    } else {
        output.dir = "rseqc_read_quality"
    }

    // The output command line option is not identical to any of the
    // produced files, so we need to work around that:
    target = output.dir + "/" + input.prefix
    product = target + ".qual.r"

    produce(product) {
        exec "$RSEQC_READ_QUALITY -i $input -o $target &>/dev/null"
    }
}

process rseqc_bam_stat = {
    var sample_dir : false

    requires RSEQC_BAM_STAT : "Must provide path to RseQC bam_stat.py script"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_bam_stat"
    } else {
        output.dir = "rseqc_bam_stat"
    }

    input_extension = ".bam"

    products = [ "$output.dir/rseqc_bam_stat.out",
                 "$output.dir/rseqc_bam_stat.err" ]

    produce(products) {
        exec "$RSEQC_BAM_STAT --input-file=$input >$output.dir/rseqc_bam_stat.out 2>$output.dir/rseqc_bam_stat.err && md5sum $outputs >$output.dir/rseqc_bam_stat.md5","rseqc"
    }
}

process rseqc_read_distribution = {
    var sample_dir : false

    requires RSEQC_READ_DISTRIBUTION : "Must provide path to RseQC read_distribution.py script"
    requires RSEQC_REF_GENE_MODEL_BED : "Must provide reference gene model in BED format"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_read_distribution"
    } else {
        output.dir = "rseqc_read_distribution"
    }

    input_extension = ".bam"

    products = [ "$output.dir/rseqc_read_distribution.out",
                 "$output.dir/rseqc_read_distribution.err" ]

    produce(products) {
        exec "$RSEQC_READ_DISTRIBUTION --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED >$output.dir/rseqc_read_distribution.out 2>$output.dir/rseqc_read_distribution.err && md5sum $outputs >$output.dir/rseqc_read_distribution.md5","rseqc"
    }
}

process rseqc_junction_annotation = {
    var sample_dir : false

    requires RSEQC_JUNCTION_ANNOTATION : "Must provide path to RseQC junction_annotation.py script"
    requires RSEQC_REF_GENE_MODEL_BED : "Must provide reference gene model in BED format"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_junction_annotation"
    } else {
        output.dir = "rseqc_junction_annotation"
    }

    input_extension = ".bam"
    product_prefix = "$output.dir/$input".replaceAll(/.*\//,"") - input_extension

    products = [ "$product_prefix" + ".junction_plot.r",
                 "$product_prefix" + ".junction.xls",
                 "$product_prefix" + ".splice_events.pdf",
                 "$product_prefix" + ".splice_junction.pdf" ]

    produce(products) {
        exec "$RSEQC_JUNCTION_ANNOTATION --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED --out-prefix=$output.dir/$product_prefix >$output.dir/rseqc_junction_annotation.out 2>$output.dir/rseqc_junction_annotation.err && md5sum $outputs >$output.dir/rseqc_junction_annotation.md5","rseqc"
    }
}

process sample_dir_prepare = {

	doc title: "A generic module to prepare the output options for a pipeline",

	desc:"""
		This pipeline module can be used to set global variables required across all subsequent modules

		Options:

		sample_dir : false / true

		If set to true, all subsequent pipeline stages will be written into a subfolder based on the
		name of the input file for a particular branch. The default is 'false'.

		For this to work, the module needs to check for branch.sample_dir.

	""",

	author: "mphoeppner@gmail.com"

	var sample_dir : false

	// Permanently store the sample name (to preserve it across forks)
	branch.sample = branch.name
	// Set the sample name as output directory
	branch.outdir = branch.sample

	branch.sample_dir = sample_dir

	if (sample_dir) {
		output.dir = branch.outdir
	}

	// produce("sample.txt") {
	//	exec "echo $input > $output"
	// }

	forward inputs
}

process indexBam = {

	doc "A function to index a BAM file"

	requires SAMTOOLS : "Must provide path to samtools"

    	transform("bam") to ("bam.bai") {
        	exec "$SAMTOOLS index $input.bam"
    	}
	forward input
}

process flagstat = {

	requires SAMTOOLS : "Must provide path to samtools"

	exec "$SAMTOOLS flagstat $input.bam > $output"
}

process samtools_filter_quality = {

	doc "Filters BAM files by quality score"

	var sample_dir : false
	var quality : "15"

	if(branch.sample_dir) { sample_dir = true }


	if (sample_dir) {
		output.dir = branch.outdir + "/bam"
	} else {
		output.dir = "bam/${branch.sample}"
	}

	requires SAMTOOLS : "Must provide path to samtools (SAMTOOLS)"

	transform("bam") to ("filtered.bam") {

		exec "$SAMTOOLS view -bq$quality -o $output $input && md5sum $output > ${output}.md5sum"

	}

}

process samtools_sort_bam = {

        doc "Sort BAM files. requires samtools version 1.2 or higher"

	requires SAMTOOLS : "Must provide path to samtools (SAMTOOLS)"

        // create output_name. Will be the input where we remove extension (should be .sam) and add .bam
        String output_name = "$input" - ~/.[^.]*$/
        String output = "${output_name}.sorted.bam"

        produce("$output") {
		exec "$SAMTOOLS sort -o $output $input && md5sum $output > ${output}.md5sum"
	}
}


process samtools_sam_to_bam = {

	doc "Convert SAM to BAM files."

	requires SAMTOOLS : "Must provide path to samtools (SAMTOOLS)"

	// create output_name. Will be the input where we remove extension (should be .sam) and add .bam
	String output_name = "$input" - ~/.[^.]*$/
	String output = "${output_name}.bam"

	produce("$output") {
		exec "$SAMTOOLS view -bS -o $output $input && md5sum $output > ${output}.md5sum"
	}
}
// Sequence conversion module

process cdna2protein = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Converts cDNA sequences to protein"

	if (sample_dir) { output.dir = branch.name }

        filter("protein") {
                exec "transeq -sequence $input -outseq $input"+".tmp -clean -trim"
                exec "sed s/_1//g $input"+".tmp > $output"
                exec "rm $input"+".tmp"
        }
}


process stringtie = {

	 doc title: "Build transcripts from aligned reads using stringtie",
                desc: """
                        Reconstructs transcripts from aligned RNA-seq reads in BAM format.

                        Input:
                        A read alignment in BAM format

                        Output:
                        transcripts.gtf

                        Requires:
			STRINGTIE : Path to the cufflinks binary

                        Stage options:
                        sample_dir : true or false, determines whether output will be written to
                                a sample-specific subfolder (true) or not (false)
                        f : set the isoform fraction cut-off (default 0.1)
			j : minimum junction coverage (default: 1)
                        G : set the reference annotation to use for guiding the assembly process (GTF/GFF3)
                """,
	        author: "jacques.dainat@bils.se"

	var sample_dir : false
        var stringtie_f : ""
        var stringtie_j : ""
	var stringtie_G : ""

	requires STRINGTIE : "Must provide path to stringtie (STRINGTIE)"

	def options = ""
	def nameOut = "default"
	// Checking which variables are set and build options string
	if (stringtie_f) {
		options += "-f $stringtie_f "
		nameOut += "f" + stringtie_f
	}
	if (stringtie_j) {
		options += "-j $stringtie_j "
		nameOut += "j" + stringtie_j
	}
        if (stringtie_G) {
                options += "-G $stringtie_G "
        }

        branch.assembly_method = nameOut

	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir + "/stringtie/" + assembly_method
        } else {
                output.dir = "stringtie/" + branch.sample + "_" + assembly_method
        }

        // The file to pass on is generally 'transcripts.gtf' - we use it as output.

	produce("transcripts.gtf") {
		uses(threads:16) {
		        exec "$STRINGTIE -l ${branch.sample} -o $output.dir"+"/transcripts.gtf -p $threads $options $input 2> $output"+".log","stringtie"
		}
	}


}
// Tophat module

process tophat = {

    doc title: "Align RNA-seq reads against a reference with tophat",
        desc: """
            Uses Tophat2/Bowtie2 to align reads against a genome

            Stage options:
                sample_dir : create sample-specific output folrder (boolean)

            Required variables are:
                PAIRED : paired-end data (boolean)
		TOPHAT : this should point to the location of the tophat executable
                LIBRARY_METHOD : This specifies the type of sequencing library, e.g. fr-reverse
                PHRED_SCORE : Specifies the quality encoding, e.g. usually this will be solexa-quals/phred33
                BWT2_INDEX : The location of a bowtie2-formatted genome index
                SAMTOOLS : This specifies the path to samtools
                BOWTIE2 : This specifies the path to bowtie2
        """,
        constraints: """
            Works with fq.gz and fq input files, but assumes paired-end
            reads (if not, set paired to false).
            The first argument is expected to be the left mate, the second
            argument needs to be the right mate
        """,
        author: "marc.hoeppner@bils.se"

    // Exposed options
    var sample_dir : false
    var tophat_r : 50       // mate inner distance
    var tophat_i : 50       // minimum intron length
    var tophat_I : 500000   // maximum intron length
    var GENOME_GTF : ""
    var TRANSCRIPTOME_INDEX : ""
    var tophat_T : false

    use_transcriptome = false

    options = "-r $tophat_r -i $tophat_i -I $tophat_I"

    // Check if an annotation file OR transcriptome index is passed and
    // modify options
    if (GENOME_GTF.length() > 0) {
        options += " -G $GENOME_GTF"
        use_transcriptome = true
    } else if (TRANSCRIPTOME_INDEX.length() > 0) {
        options += " --transcriptome-index $TRANSCRIPTOME_INDEX"
        use_transcriptome = true
    }

    // We enable quantifcation only against known transcripts but only
    // if transcripts were provided
    if (tophat_T && use_transcriptome) {
        options += " -T"
    }

    // Configuring the output directory
    if (branch.sample_dir) { sample_dir = true }

    requires TOPHAT : "Must set the TOPHAT variable to point to tophat location"
    requires LIBRARY_METHOD : "Must specify a sequencing library method (LIBRARY_METHOD)"
    requires PHRED_SCORE : "Must set a phred score (PHRED_SCORE)"
    requires BWT2_INDEX : "Must specify a Bowtie2 index (BWT2_INDEX)"
    requires BOWTIE2 : "Must specify path to Bowtie2 (BOWTIE2)"
    requires SAMTOOLS : "Must specify path to samtools (SAMTOOLS)"
    requires PAIRED : "Must specify if the sample is stranded or not (true or false)"

    // We subsequently need to keep track of folders
    // Here we set a name accessible to all subsequent modules.

    if (sample_dir) {
        output.dir = branch.outdir + "/tophat"
    } else {
        output.dir = "tophat/" + branch.sample
    }

    // If a basename for this branch was set further upstream

    produce("accepted_hits.bam") {
        uses(threads:16) {
            if (PAIRED.toBoolean()) {
                exec "$TOPHAT $PHRED_SCORE $options -o $output.dir -p $threads --library-type=$LIBRARY_METHOD $BWT2_INDEX $input1 $input2 >$output.dir/tophat.out 2>$output.dir/tophat.err && md5sum $output >$output.dir/tophat.md5","tophat"
            } else {
                exec "$TOPHAT $PHRED_SCORE $options -o $output.dir -p $threads --library-type=$LIBRARY_METHOD $BWT2_INDEX $input >$output.dir/tophat.out 2>$output.dir/tophat.err && md5sum $output >$output.dir/tophat.md5","tophat"
            }
        }
    }

    check {
        exec "[ -s $output ]"
    } otherwise {
        succeed "The Tophat output is empty. Stopping this branch ($branch.name)"
    }
}
// Trimmomatic module

process trimmomatic = {

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = true }

    doc title: "Adapter trimming of read files using Trimmomatic",
        desc: """
            Performs adapter trimming on paired-end RNA-seq reads.

            Requires:
            TM_PATH : Location of the Trimmomatic files
            TM_JAR : The name of the Trimmomatic jar file
            ADAPTER : The name of the Trimmomatic adapter file
	    PAIRED : Bolean to know if sample are paired or not
	""",
        contraints: """
            Files can be compressed (.fq.gz) or uncompressed (.fq)
        """,
        author: "mphoeppner@gmail.com"

    requires TM_PATH : "Must set TM_PATH variable to point to Trimmomatic folder"
    requires TM_JAR : "Must set TM_JAR variable to point to Trimmomatic java file"
    requires ADAPTER : "Must set the type of adapters to use"
    requires PAIRED : "Bolean to know if sample are paired or not (true or false)"

    // Determine whether to write this into a sub-folder or not

    if (sample_dir) {
        output.dir = branch.outdir + "/trimmomatic"
    } else {
        output.dir = "trimmomatic"
    }

    input_extension = ".fq.gz"

    def products
    def command

    if (PAIRED.toBoolean()) {
	println "sample is paired"
        products = [
            ("$input1".replaceAll(/.*\//,"") - input_extension + '_paired.fq.gz'),
	    ("$input2".replaceAll(/.*\//,"") - input_extension + '_paired.fq.gz'),
            ("$input1".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz'),
            ("$input2".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz')
        ]
    } else {
	println "sample is not paired"
        products = [
            ("$input".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz')
       ]
    }

    if (PAIRED.toBoolean()) {
        produce(products) {
            uses(threads:16) {
                exec "java -jar $TM_JAR PE -threads $threads $input1 $input2 ${output1} ${output3} ${output2} ${output4} ILLUMINACLIP:$TM_PATH/adapters/$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >$output.dir/trimmomatic.out 2>$output.dir/trimmomatic.err && md5sum $outputs >$output.dir/trimmomatic.md5","trimmomatic"
            }
        }
    } else {
        produce(products) {
            uses(threads:16) {
                exec "java -jar $TM_JAR SE -threads $threads $input $output ILLUMINACLIP:$TM_PATH/adapters/$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >$output.dir/trimmomatic.out 2>$output.dir/trimmomatic.err && md5sum $output >$output.dir/trimmomatic.md5","trimmomatic"
            }
        }
    }

    check {
        exec "[ -s $output1 ]"
    } otherwise {
        succeed "The Trimmomatic filtering left not reads, stopping this branch ($branch.name)"
    }
}

process verify_generic = {

	var binary : ""

	if("which $binary".execute().waitFor()!=0) {
                fail "The path provided to $binary could not be resolved"
        }

        forward inputs
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
