fasta_explode = {

  doc about: "Explodes a fasta file into its individual sequences",
    description: "Creates one file per sequence in the input FASTA file",
    constraints: "Requires exonerate to be loaded",
    author: "marc.hoeppner@bils.se"

    var directory : "sequences"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    // Doesn't produce a pre-definable output, so we use a log file as dummy target
    produce("fastaexplode.log") {
	    exec "fastaexplode -d $output.dir -f $input 2> $output"
    }

    forward input
}
