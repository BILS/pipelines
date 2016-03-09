gff_longest_cds = {

  doc about: "A module to limit gene loci in a GFF3-formatted annotation to the longest product",
    description: "Parses a GFF3-formatted annotation and filters out all transcripts except those with the longest CDS",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    requires KAHARI_CDS_FILTER : "Must provide path to Aaharis CDS filter script"


    // Defining output directory
    if (sample_dir) { 
      output.dir = branch.outdir 
    }

    transform(".gff") to (".longest_cds.gff") {
    	exec "$KAHARI_CDS_FILTER -f $input -o $output"
    }

}
