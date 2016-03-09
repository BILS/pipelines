assembly_remove_mitochondria = {

  doc about: "Parses a genome assembly and removes organellar contigs",
    description: "Matches contigs against a database of organeller proteins to find orgenellar contaminations",
    constraints: "Requires a blast database of reliable organellar proteins",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here


    // Defining output directory
//    if (sample_dir) { 
    //  output.dir = branch.outdir 
//    } else {
      
//    }

    // exec "..."
	
   forward input

}
