@transform("sam")
alignToSamPE = {
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
