// Sequence conversion module

cdna2protein = {

	var sample_dir : false
	if (branch.sample_dir) { sample_dir = true }

	doc "Converts cDNA sequences to protein"

	if (sample_dir) { output.dir = branch.name }

        filter("protein") {
                exec "transeq -sequence $input -outseq $input"+".tmp -clean -trim"
                exec "sed s/_1//g $input"+".tmp > $output"
                exec "rm $input"+".tmp"
        }
}


