// Maker module


cufflinks2maker = {
	
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
