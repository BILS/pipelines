cmsearch = {

  doc about: "Module to run Infernal's cmsearch algorithm",
    description: "Cmsearch identifies putative non-coding RNAs in nucleotide sequences",
    author: "marc.hoeppner@bils.se"

    var threads : 1
    var db : "/projects/references/databases/rfam/11.0/models_1_1/E_plus.cm"
    var directory : "rfam"

    if (branch.sample_dir) {
        output.dir = (directory.length() > 0) ? branch.outdir + "/" + directory : branch.outdir
    } else {
        if (directory.length() > 0) {
                output.dir = directory
        }
    }

    // requires here


    exec "cmsearch --cpu $threads --rfam --cut_tc --tblout $output $db $input >/dev/null"
}
