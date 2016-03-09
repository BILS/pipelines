fasta_filter_size = {

  doc about: "Parses a FASTA file and removes sequences smaller than the cutoff",
    description: "Filters fasta file by size",
    constraints: "Requires bioruby to be installed",
    author: "marc.hoeppner@bils.se"

    var size : 1000
    var directory : ""

    // requires here

    // Defining output directory
    if (branch.sample_dir) { 
	output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
      if (directory.length() > 0) {
	      output.dir = directory
      }
    }

    filter("filtered") {
	    exec "$BPIPE_BIN/fasta_filter_size.rb -i $input -s $size -o $output"
    }

}
