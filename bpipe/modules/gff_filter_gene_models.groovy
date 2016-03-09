gff_filter_gene_models = {

  doc about: "Wrapper module around A. Kaharis' gff filter_sort.pl script",
    description: "Filters transcript models by a set of criteria",
    constraints: "Requires the Kahari code base to be present",
    author: "marc.hoeppner@bils.se"

    var directory : "filter" 

    var filter_c : "-c"
    var filter_a : "-a 0.3"
    var filter_r : "-r"
    var filter_d : "-d 500"

    def options = filter_c + " " + filter_r + " " + filter_d + " " + filter_a

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here 
    requires KAHARI_GIT : "Must set path to Kahari git directory (KAHARI_GIT)"
    requires GENOME_FA : "Must provide a genome sequence in FASTA format"

    // Running a command
    
    filter("filter") {
	    exec "perl $KAHARI_GIT/scripts/GFF/filter_sort.pl -f $input -F $GENOME_FA -o $output $options"
    }
}
