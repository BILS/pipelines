blast_recursive = {

  doc about: "A model to blast a sequence file against a db of itself",
    description: "Takes a sequence file and blasts it against a database of the same name (recursive)",
    constraints: "Expects a blast database to exist for the input file and in the same location",
    author: "marc.hoeppner@bils.se"

    var directory : "blast_recursive"
    var blast_prog : "blastp"
    var blast_outfmt : "5"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command
    
    exec "$blast_prog -query $input -db $input -num_threads $threads -outfmt $blast_outfmt -out $output"

}
