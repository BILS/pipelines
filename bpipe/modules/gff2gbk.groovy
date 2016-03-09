gff2gbk = {

    doc "Converts a GFF3 formatted file into GenBank format"

    var directory : "gff2genbank"
    var flank : 1000

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here
    requires GENOME_FA : "Must provide a genome sequence in fasta format (GENOME_FA)"
    
    // Running a command

    transform("gbk") {
        exec "gff2gbSmallDNA.pl $input $GENOME_FA $flank $output"
    }
}
