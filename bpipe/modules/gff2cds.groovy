gff2cds = {

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        doc "Extracts cDNA sequences from the annotation/genome"

        requires GENOME_FA : "Must set variable GENOME_FA"

        if (sample_dir) {
                output.dir = branch.outdir
        }

        transform(".gff") to (".cds.fa") {
                exec "gffread -x $output -g $GENOME_FA $input"
        }
}

