run_blast2go = {
	
	doc title: "Runs the Blast2Go pipeline using blast and interpro XML inputs"

	var sample_dir : false

	if (branch.sample_dir) { sample_dir = true }

	if (branch.sample_dir) { output.dir = branch.outdir }
	
	// Assumes that the Blast output is input1 and interpro is input2 - BAD IDEA...

	produce(branch.sample+"_blast2go.annot") {
		exec "java -Xmx20G -cp ${B2G4PIPEPATH}/*:${B2G4PIPEPATH}/ext/* es.blast2go.prog.B2GAnnotPipe -in $input1 -out $output.prefix -prop ${B2G4PIPEPATH}/b2gPipe.properties -ips $input2 -annot"
	}

}


blast2go2gff = {

	var sample_dir : false

	if (branch.sample_dir) { sample_dir = true }

	doc title: "Updates a GFF file with meta data from Blast2Go"
	
	requires B2GOGFF_SCRIPT : "Specify the full path to the b2gogff script"
	
	if (branch.sample_dir) { output.dir = branch.outdir }

	produce(gff_file_with_ids.prefix+".description.gff") {
		exec "$B2GOGFF_SCRIPT --gff $gff_file_with_ids --b2go $input > $output"
	}

}
