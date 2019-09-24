// A pipeline to create to assemble transcripts for annotation

// about title: "A pipeline to assemble transcripts from RNAseq data based on Tophat/Stringtie (input ...)"

// load 'pipeline.config'
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	trimmomatic()
	tophat()
	stringtie()
}

// run { ~"(.*)_.+.f*q*[.gz]?" *
// 	[ verify_generic.using(binary:"tophat") + verify_generic.using(binary:"stringtie") + sample_dir_prepare.using(sample_dir:true) +
// 		trimmomatic +
// 		tophat +
// 		stringtie
// 	]
// }
