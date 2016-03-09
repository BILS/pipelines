sample_dir_prepare = {
	
	doc title: "A generic module to prepare the output options for a pipeline",
	
	desc:"""
		This pipeline module can be used to set global variables required across all subsequent modules

		Options:

		sample_dir : false / true 

		If set to true, all subsequent pipeline stages will be written into a subfolder based on the
		name of the input file for a particular branch. The default is 'false'.

		For this to work, the module needs to check for branch.sample_dir. 
		
	""",

	author: "mphoeppner@gmail.com"

	var sample_dir : false

	// Permanently store the sample name (to preserve it across forks)
	branch.sample = branch.name
	// Set the sample name as output directory
	branch.outdir = branch.sample

	branch.sample_dir = sample_dir

	if (sample_dir) {
		output.dir = branch.outdir
	} 

	// produce("sample.txt") {
	//	exec "echo $input > $output"
	// }
	
	forward input
}




	
