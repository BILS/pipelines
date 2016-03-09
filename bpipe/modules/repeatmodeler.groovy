repeatmodeler_format_genome = {

  doc about: "A module to generate the input for RepeatModeler",
    description: "Formats a genome sequence to be used by RepeatModeler",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here

    // Defining output directory
    if (sample_dir) { 
      output.dir = branch.outdir + "/repeats"
    } else {
      output.dir = "repeats"
    }

    transform("nhr") {
    	exec "BuildDatabase -name $input.prefix -engine ncbi $input.fa"
    }

    forward input
}


repeatmodeler_run = {

    doc about: "A module to model repeats in a genome assembly",
    description: "Runs RepatModeler on a genome assembly to predict repeats",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here

    // Defining output directory
    if (sample_dir) {
      output.dir = branch.outdir + "/repeats"
    } else {
      output.dir = "repeats"
    }
    uses(threads:16) {
    	exec "RepeatModeler -database $input.prefix -engine ncbi -pa $threads && cp RM_*/consensi.fa.classified $output 2>/dev/null"
    }

}


repeatmodeler = segment {

        repeatmodeler_format_genome + repeatmodeler_run

}
