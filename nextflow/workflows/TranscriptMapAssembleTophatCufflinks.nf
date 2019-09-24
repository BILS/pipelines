// A pipeline to create to assemble transcripts for annotation

// about title: "A pipeline to assemble transcripts from RNAseq data based on Cufflinks/Tophat"
//
// inputs "fq.gz" : "Requires compressed FastQ file(s) as input"
//
// load 'pipeline.config'
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	trimmomatic()
	tophat()
	cufflinks()
	cufflinks2maker()
}

// run { ~"(.*)_.+.f*q*[.gz]?" *
// 	[ verify_generic.using(binary:"tophat") + verify_generic.using(binary:"cufflinks") + sample_dir_prepare.using(sample_dir:true) +
// 		trimmomatic +
// 		tophat +
// 		cufflinks.using(cufflinks_j:0.1, cufflinks_F:0.15) +
// 		cufflinks2maker
// 	]
// }
