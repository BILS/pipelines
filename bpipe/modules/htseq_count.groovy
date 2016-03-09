htseq_count = {

	doc title: "HTseq-count counts reads overlapping an annotation",

	desc: """
		HTseq-count is a python tool designed to accept a read alignment in BAM or SAM
		format and counts the reads overlapping annotated exons provided in GTF format.
	""",

	constraints: "Requires a read alignment in BAM/SAM format and GTF-formatted annotation",

	author: "marc.hoeppner@bils.se"

	var sample_dir : true
	
	requires GENOME_GTF : "Must provide an annotation in GTF format (GENOME_GTF)"
	requires HTSEQ_COUNT : "Must provide path to htseq-count (HTSEQ_COUNT)"

	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir + "/htseq_count"
	} else {
		output.dir = "htseq_count"
	}

	transform(".bam") to (".htseqcount") {
		exec "HTSEQ_COUNT -f bam -t exon $input.bam $GENOME_GTF > $output"
	}

}
