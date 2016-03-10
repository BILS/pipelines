gff_longest_cds = {

  doc about: "A module to limit gene loci in a GFF3-formatted annotation to the longest product",
    description: "Parses a GFF3-formatted annotation and filters out all transcripts except those with the longest CDS",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    requires KAHARI_GIT : "Must set path to Kahari git directory (KAHARI_GIT)"


    // Defining output directory
    if (sample_dir) { 
      output.dir = branch.outdir 
    }

    transform(".gff") to (".longest_cds.gff") {
    	exec "$KAHARI_GIT/scripts/GFF/find_longest_CDS.pl -f $input -o $output"
    }

}
