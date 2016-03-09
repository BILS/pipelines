gffread_extract_sequences = {

	doc about: "A module to extract CDS and protein sequences from an annotation using gffread",
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires GENOME_FA : "Must provide genome sequence in FASTA format (GENOME_FA)"

	if (sample_dir) { output.dir = branch.outdir }
	
	produce(input.prefix+".cds.fa",input.prefix+".proteins.fa") {
		exec """
			gffread -x $output1 -g $GENOME_FA $input ; gffread -y $output2 -g $GENOME_FA $input 
		"""
	}
}
