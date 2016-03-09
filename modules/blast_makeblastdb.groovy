blast_makeblastdb = {

  doc about: "Transforms a sequence file into a blast+ database",
    description: "Runs makeblastdb to generate a blast database",
    constraints: "Assumes protein input unless specified otherwise (dbtype)",
    author: "marc.hoeppner@bils.se"

    var dbtype : "prot"
    
    // requires here

    // Running a command
    produce(input + ".phr") {
	    exec "makeblastdb -in $input -dbtype $dbtype"
    }

    forward input

}
