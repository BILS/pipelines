blast = {

        doc title: "Run a blast+ search on a FASTA file",

        desc: "Takes a sequence fasta file as input and blasts it against a database",

        author: "marc.hoeppner@bils.se"

        var sample_dir : false
        var outfmt : 5
        var program : "blastp"
        var max_target_seqs : 10

        requires BLASTP_DB : "Must specify a blast+ formatted protein database (BLASTP_DB)"

        if (branch.sample_dir) { sample_dir = true}

        if (sample_dir) {
                output.dir = branch.outdir + "/theVoid"
        }
        uses(threads:8) {
                exec "$program -db $BLASTP_DB -query $input -outfmt $outfmt -max_target_seqs $max_target_seqs -num_threads $threads -out $output","blast"
        }
}

blastp = {

	doc title: "Run a blastp+ search on a protein FASTA file",

	desc: "Takes a protein fasta file as input and blasts it against a protein database",

	author: "marc.hoeppner@bils.se"

	var sample_dir : false
	var outfmt : 5

	requires BLASTP_DB : "Must specify a blast+ formatted protein database (BLASTP_DB)"
        requires BLASTP : "Must provide path to blastp binary"
	
	if (branch.sample_dir) { sample_dir = true}

	if (sample_dir) {
		output.dir = branch.outdir + "/theVoid"
	}
	uses(threads:8) {
	        exec "$BLASTP -db $BLASTP_DB -query $input -outfmt $outfmt -num_threads $threads -out $output","blastp"
	}
}

recursive_blastp = {

	doc "Runs a blastp search against itself (blast+)"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) { output.dir = branch.outdir + "/theVoid" }

	requires BLASTP : "Must provide path to blastp binary"

	exec "$BLASTP -db $input -query $input.fa -outfmt 5 -num_threads $threads -out $output"

}

protein2blast_db = {

	doc "Creates a blast database from a protein fasta file (blast+)"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir + "/theVoid"
	}

	produce(input+".phr") {
		exec "makeblastdb -in $input.fa -dbtype prot"
	}

	forward input
}

recursive_protein2blastp = segment { protein2blast_db + recursive_blastp }

merge_blast_xml = {

	doc "Crude method to merge the XML output from multiple BLAST searches"

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir
	}
	
	produce(branch.sample + "_blast.out") {
		exec "cat $inputs > $output"	
	}

}

merge_blast_tab = {

	doc "Crude method to merge the tabular output from multiple BLAST searches"

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (sample_dir) {
                output.dir = branch.outdir
        }

        produce(branch.sample + "_blast.out") {
                exec "cat $inputs > $output"
        }

}
