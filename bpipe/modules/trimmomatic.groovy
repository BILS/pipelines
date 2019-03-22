// Trimmomatic module

trimmomatic = {

    var sample_dir : false
    if (branch.sample_dir) { sample_dir = true }

    doc title: "Adapter trimming of read files using Trimmomatic",
        desc: """
            Performs adapter trimming on paired-end RNA-seq reads.

            Requires:
            TM_PATH : Location of the Trimmomatic files
            TM_JAR : The name of the Trimmomatic jar file
            ADAPTER : The name of the Trimmomatic adapter file
	    PAIRED : Bolean to know if sample are paired or not
	""",
        contraints: """
            Files can be compressed (.fq.gz) or uncompressed (.fq)
        """,
        author: "mphoeppner@gmail.com"

    requires TM_PATH : "Must set TM_PATH variable to point to Trimmomatic folder"
    requires TM_JAR : "Must set TM_JAR variable to point to Trimmomatic java file"
    requires ADAPTER : "Must set the type of adapters to use"
    requires PAIRED : "Bolean to know if sample are paired or not (true or false)"

    // Determine whether to write this into a sub-folder or not

    if (sample_dir) {
        output.dir = branch.outdir + "/trimmomatic"
    } else {
        output.dir = "trimmomatic"
    }

    input_extension = ".fq.gz"

    def products
    def command

    if (PAIRED.toBoolean()) {
	println "sample is paired"
        products = [
            ("$input1".replaceAll(/.*\//,"") - input_extension + '_paired.fq.gz'),
	    ("$input2".replaceAll(/.*\//,"") - input_extension + '_paired.fq.gz'),
            ("$input1".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz'),
            ("$input2".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz')
        ]
    } else {
	println "sample is not paired"
        products = [
            ("$input".replaceAll(/.*\//,"") - input_extension + '_unpaired.fq.gz')
       ]
    }

    if (PAIRED.toBoolean()) {
        produce(products) {
            uses(threads:16) {
                exec "java -jar $TM_JAR PE -threads $threads $input1 $input2 ${output1} ${output3} ${output2} ${output4} ILLUMINACLIP:$TM_PATH/adapters/$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >$output.dir/trimmomatic.out 2>$output.dir/trimmomatic.err && md5sum $outputs >$output.dir/trimmomatic.md5","trimmomatic"
            }
        }
    } else {
        produce(products) {
            uses(threads:16) {
                exec "java -jar $TM_JAR SE -threads $threads $input $output ILLUMINACLIP:$TM_PATH/adapters/$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >$output.dir/trimmomatic.out 2>$output.dir/trimmomatic.err && md5sum $output >$output.dir/trimmomatic.md5","trimmomatic"
            }
        }
    }

    check {
        exec "[ -s $output1 ]"
    } otherwise {
        succeed "The Trimmomatic filtering left not reads, stopping this branch ($branch.name)"
    }
}
