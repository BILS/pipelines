annie_interpro = {

  doc about: "Module to run Annie's InterPro converter",
    description: "Annie is a tool to convert search results into tabular format for functional annotation",
    constraints: "Requires the annie package to be available",
    author: "marc.hoeppner@bils.se"

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = branch.sample_dir }

    // requires here
    requires ANNIE_ROOT : "Must provide path to Annie installation folder (ANNIE_ROOT)"

    // Defining output directory
    if (sample_dir) { 
      output.dir = branch.outdir 
    } 

    exec "python $ANNIE_ROOT/annie.py ipr $input $output"

}
