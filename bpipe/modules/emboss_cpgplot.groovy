emboss_cpgplot = {

  doc about: "A module to run Emboss' cpgplot tool",
    description: "cpgplot predicts cpg islands in genomic sequences",
    constraints: "Should only be used for vertebrate genomes",
    author: "marc.hoeppner@bils.se"

    var cpgisland_window : 100
    var cpgisland_minlen : 200
    var cpgisland_minoe : 0.6
    var cpgisland_minpc : 50.0
    
    var directory : "cpgplot"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    produce($input+".cpg.gff") {
    	exec "cpgplot -sequence $input -outfeat $output -outfile $output" + ".out -window $cpgisland_window -minlen $cpgisland_minlen -minoe $cpgisland_minoe -minpc $cpgisland_minpc -noplot -nocg"
    }
}
