// Hisat2 module

hisat2 = {

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
