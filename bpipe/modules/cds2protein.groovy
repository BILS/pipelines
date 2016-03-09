cds2protein = {

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        doc "Converts cDNA sequences to protein"

        if (sample_dir) { output.dir = branch.name }

        filter("protein") {
		exec "transeq -sequence $input -outseq $output -clean -trim"
		//exec "$BPIPE_BIN/cds2protein.pl --infile $input --outfile $output 2>/dev/null"
        }
}
