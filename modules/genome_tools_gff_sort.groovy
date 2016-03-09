genome_tools_gff_sort = {

	doc about: "A module to sort a GFF3-formatted annotation by coordinates using GenomeTools",
	
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires GENOME_TOOLS : "Must provide path to genometools (GENOMETOOLS)"

	if (sample_dir) { output.dir = branch.outdir }

	filter("sorted") {
		exec "$GENOME_TOOLS gff3 -sort $input > $output 2>/dev/null"	
	} 

}
