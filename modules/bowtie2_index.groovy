bowtie2_index = {

  doc about: "A pipeline to generate a bowtie2 index from a genome sequence",
    description: "Takes a genome sequence in fasta format and generates a bowtie2 index",
    author: "marc.hoeppner@bils.se"

    var directory : "bowtie-index"

    if (branch.sample_dir) { 
	output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
	if (directory.length() > 0) {
		output.dir = directory
	}
    }
	
    // requires here

    produce(input.prefix + ".1.bt2") {
    	exec "bowtie2-build $input $output.prefix.prefix >/dev/null"
    }

    forward input

}
