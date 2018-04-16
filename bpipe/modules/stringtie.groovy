stringtie = {

	 doc title: "Build transcripts from aligned reads using stringtie",
                desc: """
                        Reconstructs transcripts from aligned RNA-seq reads in BAM format.

                        Input:
                        A read alignment in BAM format

                        Output:
                        transcripts.gtf

                        Requires:
			STRINGTIE : Path to the cufflinks binary

                        Stage options:
                        sample_dir : true or false, determines whether output will be written to
                                a sample-specific subfolder (true) or not (false)
                        f : set the isoform fraction cut-off (default 0.1)
			j : minimum junction coverage (default: 1)
                        G : set the reference annotation to use for guiding the assembly process (GTF/GFF3)
                """,
	        author: "jacques.dainat@bils.se"

	var sample_dir : false
        var stringtie_f : ""
        var stringtie_j : ""
	var stringtie_G : ""

	requires STRINGTIE : "Must provide path to stringtie (STRINGTIE)"

	def options = ""
	def nameOut = "default"
	// Checking which variables are set and build options string
	if (stringtie_f) {	
		options += "-f $stringtie_f "
		nameOut += "f" + stringtie_f
	} 
	if (stringtie_j) {
		options += "-j $stringtie_j "
		nameOut += "j" + stringtie_j
	}
        if (stringtie_G) {
                options += "-G $stringtie_G "
        }

        branch.assembly_method = nameOut
	
	if (branch.sample_dir) { sample_dir = true }
 	
	if (sample_dir) {
		output.dir = branch.outdir + "/stringtie/" + assembly_method
        } else {
                output.dir = "stringtie/" + branch.sample + "_" + assembly_method
        }

        // The file to pass on is generally 'transcripts.gtf' - we use it as output.

	produce("transcripts.gtf") {
		uses(threads:16) {
		        exec "$STRINGTIE -l ${branch.sample} -o $output.dir"+"/transcripts.gtf -p $threads $options $input 2> $output"+".log","stringtie"
		}
	}


}
