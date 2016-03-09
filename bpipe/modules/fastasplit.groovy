// File operation module

fastasplit = {

	doc title: "Splits a FASTA file into chunks, using the CHUNKS variable",
		desc: """
			Splits a fasta-formatted file into chunks of similar
			size using the exonerate tool fastasplit.

			Input:
			A fasta-formatted sequence file

			Requires:
			CHUNKS : defines the number of chunks to split into.

		""",

		constraints: """
			The number of chunks should be chosen so that at least 
			10 sequences are included in each chunk - else the splitting
			may fail.
		""",
		author: "marc.hoeppner@bils.se"

        requires CHUNKS : "Must set variable CHUNKS"

	if (branch.sample_dir) { output.dir = branch.outdir } 	

	def chunkfiles = []

	for ( value in (0..(CHUNKS.toInteger()-1)) ) {
		if (value < 10) {
			chunkfiles.push(input+"_chunk_000000${value}")
		} else if (value < 100) {
			chunkfiles.push(input+"_chunk_00000${value}")
		} else if (value < 1000) {
			chunkfiles.push(input+"_chunk_0000${value}")
		} else {
			chunkfiles.push(input+"_chunk_000${value}")
		}
	}

      	//  def chunkfiles = (0..(CHUNKS.toInteger()-1)).collect{ input+"_chunk_000000${it}" }

	produce(chunkfiles) {
                exec "fastasplit -f $input -o ${output.dir} -c $CHUNKS"
        }

}

