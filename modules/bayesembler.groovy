bayesembler = {

  doc about: "Module to run the Bayesembler transcriptome assembler",
    description: "Bayesembler uses a bayesian approach to reconstruct transcripts from BAM-formatted read alignments",
    constraints: "Requires a Tophat-generated BAM file and index",
    author: "marc.hoeppner@bils.se"

    var baysembler_c : "0.5"
    var bayesembler_m : "-m"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    produce("assembly.gtf") {
	    exec "bayesembler -b $input -p $threads $baysembler_m -c $bayesembler_c"
    }

}
