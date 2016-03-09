genome_tools_gff_to_gtf = {

	doc about: "A module to convert a GFF3 file to GTF using genometools",
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires GENOME_TOOLS : "Must provide path to genometools (GENOME_TOOLS)"

	if (sample_dir) { output.dir = branch.outdir }

	transform(".gff") to (".gtf") {
		exec "$GENOME_TOOLS gff3_to_gtf -o $output -force $input 2>/dev/null"
	}
}
