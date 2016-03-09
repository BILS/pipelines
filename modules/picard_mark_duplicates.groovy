dedup = {

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
