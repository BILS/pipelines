rseqc_read_quality = {
    doc title: "RseQC - Analyzing read quality",

    desc: """
        RseQC read_quality.py - analyzing and plotting
        the quality of RNA-seq reads
    """,

    constraints: "Requires a BAM file as input",

    author: "marc.hoeppner@bils.se"

    var paired : true
    var sample_dir : false

    requires RSEQC_READ_QUALITY : "Must provide path to RseQC read_quality.py script"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_read_quality"
    } else {
        output.dir = "rseqc_read_quality"
    }

    // The output command line option is not identical to any of the
    // produced files, so we need to work around that:
    target = output.dir + "/" + input.prefix
    product = target + ".qual.r"

    produce(product) {
        exec "$RSEQC_READ_QUALITY -i $input -o $target &>/dev/null"
    }
}

rseqc_bam_stat = {
    var sample_dir : false

    requires RSEQC_BAM_STAT : "Must provide path to RseQC bam_stat.py script"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_bam_stat"
    } else {
        output.dir = "rseqc_bam_stat"
    }

    input_extension = ".bam"

    products = [ "$output.dir/rseqc_bam_stat.out",
                 "$output.dir/rseqc_bam_stat.err" ]

    produce(products) {
        exec "$RSEQC_BAM_STAT --input-file=$input >$output.dir/rseqc_bam_stat.out 2>$output.dir/rseqc_bam_stat.err && md5sum $outputs >$output.dir/rseqc_bam_stat.md5","rseqc"
    }
}

rseqc_read_distribution = {
    var sample_dir : false

    requires RSEQC_READ_DISTRIBUTION : "Must provide path to RseQC read_distribution.py script"
    requires RSEQC_REF_GENE_MODEL_BED : "Must provide reference gene model in BED format"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_read_distribution"
    } else {
        output.dir = "rseqc_read_distribution"
    }

    input_extension = ".bam"

    products = [ "$output.dir/rseqc_read_distribution.out",
                 "$output.dir/rseqc_read_distribution.err" ]

    produce(products) {
        exec "$RSEQC_READ_DISTRIBUTION --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED >$output.dir/rseqc_read_distribution.out 2>$output.dir/rseqc_read_distribution.err && md5sum $outputs >$output.dir/rseqc_read_distribution.md5","rseqc"
    }
}

rseqc_junction_annotation = {
    var sample_dir : false

    requires RSEQC_JUNCTION_ANNOTATION : "Must provide path to RseQC junction_annotation.py script"
    requires RSEQC_REF_GENE_MODEL_BED : "Must provide reference gene model in BED format"

    if (branch.sample_dir) { sample_dir = true }

    if (sample_dir) {
        output.dir = branch.outdir + "/rseqc_junction_annotation"
    } else {
        output.dir = "rseqc_junction_annotation"
    }

    input_extension = ".bam"
    product_prefix = "$output.dir/$input".replaceAll(/.*\//,"") - input_extension

    products = [ "$product_prefix" + ".junction_plot.r",
                 "$product_prefix" + ".junction.xls",
                 "$product_prefix" + ".splice_events.pdf",
                 "$product_prefix" + ".splice_junction.pdf" ]

    produce(products) {
        exec "$RSEQC_JUNCTION_ANNOTATION --input-file=$input --refgene=$RSEQC_REF_GENE_MODEL_BED --out-prefix=$output.dir/$product_prefix >$output.dir/rseqc_junction_annotation.out 2>$output.dir/rseqc_junction_annotation.err && md5sum $outputs >$output.dir/rseqc_junction_annotation.md5","rseqc"
    }
}
