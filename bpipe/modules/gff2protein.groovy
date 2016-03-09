gff2protein = {

    doc about: "A module to extract protein sequences from GFF annotations",
    description: "Reports protein sequences from GFF annotations",
    constraints: "Only works with standard eukaryotic genetic code!",
    author: "marc.hoeppner@bils.se"

    var directory : "protein"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here
    requires GENOME_FA : "Must provide a genome sequence (GENOME_FA)"

    // Running a command

    transform(".gff") to (".proteins.fa") {
            exec "gffread -y $input.prefix"+".tmp -g $GENOME_FA $input && $BPIPE_BIN/fix_fasta.rb $input.prefix"+".tmp > $output && rm $input.prefix"+".tmp"
    }

}
