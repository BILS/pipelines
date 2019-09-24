// about title: "Please add a tile"
//
// inputs "fq.gz" : "Requires read files in zipped FastQ format (fq.gz)"
//
// load 'pipeline.config'
//
nextflow.preview.dsl=2
include './../modules/annotation_modules'

workflow {
	data = Channel.fromPath(params.reads)
	trimmomatic
	tophat
	rseqc_bam_stat
	rseq_read_distribution
	rseqc_juntion_annotation
}

// run {
//     "%_*.fq.gz" * [
//         sample_dir_prepare.using(sample_dir:true)
//         + trimmomatic.using(sample_dir:true,paired:false,threads:8)
//         + tophat.using(sample_dir:true,paired:false,threads:8)
//         + [ rseqc_bam_stat.using(sample_dir:true),
//             rseqc_read_distribution.using(sample_dir:true),
//             rseqc_junction_annotation.using(sample_dir:true) ]
//     ]
// }
