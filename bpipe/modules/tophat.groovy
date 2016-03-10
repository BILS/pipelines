// Tophat module

tophat = {

    doc title: "Align RNA-seq reads against a reference with tophat",
        desc: """
            Uses Tophat2/Bowtie2 to align reads against a genome

            Stage options:
                sample_dir : create sample-specific output folrder (boolean)
                paired : paired-end data (boolean)

            Required variables are:
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

    var paired : true // paired-end data?

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
            if (paired) {
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
