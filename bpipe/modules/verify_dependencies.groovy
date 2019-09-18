verify_generic = {

	var binary : ""

	if("which $binary".execute().waitFor()!=0) {
                fail "The path provided to $binary could not be resolved"
        }

        forward inputs
}
	

// Verifications do not require an execute, so they will run on the head node and
// can be grouped as a segment (else it would mean submitting 4 or more jobs...)

verify_dependencies_rnaseq = segment {
	verify_generic.using(binary:"tophat") + verify_generic.using(binary:"cufflinks") + verify_generic.using(binary:"samtools") + verify_generic.using(binary:"bowtie2") 
}

verify_annotation_preprocess = segment {
	verify_generic.using(binary:"fastasplit") + verify_generic.using(binary:"Rscript") + verify_generic.using(binary:"cpgplot") + verify_generic.using(binary:"bowtie2-build")
}

verify_dependencies_annotation_models = segment {
	verify_generic.using(binary:"gffread") + verify_generic.using(binary:"makeblastdb") + verify_generic.using(binary:"gff2gbSmallDNA.pl")
}


