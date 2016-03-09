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
