gff_annotation_stats = {

	doc about: "A module to generate annotation statistics, based on the code by A. Kähäri",
	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	requires KAHARI_ANNOTATION_STATS : "Must provide path to annotation_stats.pl (KAHARI_ANNOTATION_STATS)"
	requires GENOME_FA : "Must provide path to genome sequence in FASTA format (GENOME_FA)"

	if (sample_dir) { output.dir = branch.outdir }

	transform(".gff") to (".stats") {
		exec "$KAHARI_ANNOTATION_STATS --gff $input --fasta $GENOME_FA > $output"
	}


}
