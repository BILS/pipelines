// Interproscan Module


interpro = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Runs a protein fasta file against InterPro"

	if (sample_dir) { output.dir = branch.outdir }

	requires INTERPROSCAN : "Must specify the location of interproscan (INTERPROSCAN)"

	produce(input+".gff3",input+".tsv",input+".xml") {
	        exec "$INTERPROSCAN -appl PfamA-27.0,ProDom-2006.1,SuperFamily-1.75,PIRSF-2.84 -i $input -d ${output.dir} -iprlookup -goterms -pa -dp > /dev/null 2> /dev/null ","interpro"
        }
}


merge_interpro_xml = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Merges the XML output from multiple InterPro searches"

	if (sample_dir) { output.dir = branch.outdir }

	produce(branch.sample+"_interpro.xml") {

		// This is a really stupid hack to merge XML files
		// We take the first lines and the last line from the first file
		// And squeeze the other stuff from all files in between

		def first_file = inputs[0]
		def temp_file = branch.sample + "_interpro.tmp"

		exec "head -n 2 $first_file > $temp_file"

		for (i in inputs) {
			exec "grep -v xml $i | grep -v protein-matches >> $temp_file"
		}
		
		exec "tail -n 1 $first_file >> $temp_file"

		exec "mv $temp_file $output"
	}

	branch.ipr_xml = output
}

merge_interpro_tsv = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Merges the TSV output from multiple InterPro searches"

	if (sample_dir) { output.dir = branch.outdir }

	produce(branch.sample+"_interpro.tsv") {

		exec "cat $inputs >> $output"

	}

	// We save the output name since we need it for a later pipeline stage.

	branch.iprtsv = output	
	
	// Pass on the name of the output file


}


interpro2gff = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Converts InterPro results into GFF format"

	if (sample_dir) { output.dir = branch.outdir }


	// This is a horrible way of doing this...
	// We need the name of a file set in an earlier stage,
	// but the variable doesn't carry over, so:

	if (sample_dir) {
		iprtsv = branch.outdir + "/" + branch.sample + "_interpro.tsv"
	} else {
		iprtsv = branch.sample + "_interpro.tsv"
	}

	filter("interpro") {
		exec "ipr_update_gff $input $iprtsv > $output"
	}
}
