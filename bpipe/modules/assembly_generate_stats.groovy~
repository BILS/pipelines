assembly_generate_stats = {

  doc about: "Generates statistics from a genome assembly",
    description: "Accepts a multi-fasta file with nucleotide sequences and computes statistics",
    author: "marc.hoeppner@bils.se"

    // requires here

    var directory : "stats"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }
    
    produce("assembly_report.txt") {
    	exec "$BPIPE_BIN/fasta_statisticsAndPlot.pl --infile $input --outfile $output"
    }
}
