gff_filter_by_blast = {

  doc about: "Remove GFF entries if they are found in a recursive blast (redundancy removal)",
    description: "Takes a recursive blast input (format 6) and removes redundant entries from GFF file",
    constraints: "Requires the blast file to be in outfmt 6",
    author: "marc.hoeppner@bils.se"

    var directory : "nonredundant"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    blast_file = input
    // input will change due to the expression below, so we capture it here

    produce(branch.name + ".nr.gff") {
	from("*.gff") {
                exec "perl $BPIPE_BIN/gff_filter_by_mrna_id.pl --gff $input --blast $blast_file --outfile $output"
	}
    }
}
