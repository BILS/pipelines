gbk2augustus = {

    doc about: "A generic module that needs a description",
    description: "Description here",
    constraints: "Information on constraints here",
    author: "marc.hoeppner@bils.se"

    var test_size : 100

    // requires here

    // Running a command
    produce(input+".train") {
	    exec "randomSplit.pl $input $test_size"
    }

}
