genome_sanity_check = {

  doc about: "A module to analyse and produce statistics for a genome assembly",
    description: "Verifies the integrity of a genome assembly for annotation",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // Defining output directory
    if (sample_dir) { 
      output.dir = branch.outdir
    } 

   // exec "..."

   forward input
}
