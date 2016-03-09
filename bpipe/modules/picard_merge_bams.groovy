@filter("merge")
mergeBams = {

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
