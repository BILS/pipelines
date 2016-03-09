@transform("bam")
samToSortedBam = {
	
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

