gff_get_trna = {
	
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
