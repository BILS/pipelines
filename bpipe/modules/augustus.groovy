gbk2training = {
	
	doc title: "Converts a GenBank file to an Augustus test and training set",

	desc: """
		Splits a Genbank file into two random sets for training and testing of
		Augustus profile models. 

		Required variables:
		TEST_SIZE : How many models are to be kept for testing.
	""",
	constraints : "None",

	author : "marc.hoeppner@bils.se"
		
	var sample_dir : false

	requires TEST_SIZE : "Must provide a value for the size of the test data set (TEST_SIZE)"


	if (branch.sample_dir) { sample_dir = true }

	if (sample_dir) {
		output.dir = branch.outdir
	}

	produce(input+".train") {
		exec "randomSplit.pl $input $TEST_SIZE"
	}

}
