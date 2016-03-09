cegma = {

  doc about: "A wrapper around the Cegma package",
    description: "Cegma checks the expected gene space coverage of an assembly",
    constraints: "Requires the Cegma package",
    author: "marc.hoeppner@bils.se"

    var directory : "cegma"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here

    // Running a command

    produce("cegma.completeness_report") {
	uses(threads:16) {
	    	exec "cegma -T $threads -g $input -o $output.dir/cegma"
	}
    }

}
