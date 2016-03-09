gff_stable_id = {

        var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        doc title: "Generates stable IDs for all features in the GFF file"

        requires ID_TRUNK : "Must set variable ID_TRUNK"

        if (sample_dir) {
                output.dir = branch.outdir
        }

        filter("stable_id") {
                exec "$BPIPE_BIN/gff_create_stable_id.pl --gff $input --id_trunk $ID_TRUNK --outfile $output"
        }

        branch.gff_file_with_ids = output

}

