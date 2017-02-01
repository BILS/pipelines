gff2protein = {

    doc about: "A module to extract protein sequences from GFF annotations",
    description: "Reports protein sequences from GFF annotations",
    author: "marc.hoeppner@bils.se, jacques.dainat@bils.se"

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
    requires CODON_TABLE :  "Must specify the translation table (CODON_TABLE)"

    // Running a command
    transform(".gff") to (".proteins.fa") {
            exec "$BPIPE_BIN/gff3_sp_extract_sequences.pl -o $input.prefix"+".tmp -f $GENOME_FA -p -cfs -cis -ct $CODON_TABLE --gff $input && $BPIPE_BIN/fix_fasta.rb $input.prefix"+".tmp.fa > $output && rm $input.prefix"+".tmp.fa"
    }

}
