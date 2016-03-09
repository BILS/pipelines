realignIntervals = {

	// Hard-coded to take 2 known indels files right now

	doc "Identify realign intervals with GATK"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
    	requires REF : "Must provide a path to reference (REF)"
	requires GOLD_STANDARD_INDELS : "Must provide a path to gold standard indels (GOLD_STANDARD_INDELS)"
	requires INDELS_100G : "Must provide a path to indels 100G (INDELS_100G)"
	requires LOG : "Must provide a location for log file (LOG)"

	exec """
        	java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I $input.bam --known $GOLD_STANDARD_INDELS --known $INDELS_100G -log $LOG -o $output.intervals
   	 """
}

realign = {
	
	doc "Realign indels with GATK"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"

    	exec """
        	java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $input.bam -targetIntervals $input.intervals -log $LOG -o $output.bam
    	"""
}


baseQualRecalCount = {

    	doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"
	requires DBSNP : "Must provide a path to known SNP DB (DBSNP)"
    
	exec "java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $input.bam -R $REF --knownSites $DBSNP -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -log $LOG -o $output.counts"
}

baseQualRecalTabulate = {

	doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/align"
        } else {
                output.dir="align"
        }    

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"

	exec "java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -I $input.bam -BQSR $input.counts -R $REF -l INFO -log $LOG -o $output"
}

callSNPs = {
    	doc "Call SNPs/SNVs using GATK Unified Genotyper"
    
	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }
   
	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"
	requires DBSNP : "Must provide path to refence SNP DB (DBSNP)"

 	exec """
            java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
               -nt $threads 
               -R $REF 
               -I $input.bam 
               --dbsnp $DBSNP 
               -stand_call_conf 50.0 -stand_emit_conf 10.0 
               -dcov 1600 
               -l INFO 
               -A AlleleBalance -A DepthOfCoverage -A FisherStrand 
               -glm SNP -log $LOG 
               -o $output.vcf
        """
}

callIndels = {

	doc "Call variants using GATK Unified Genotyper"
    	
	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }

	requires GATK : "Must provide a path to GATK (GATK)"
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"
	requires DBSNP : "Must provide path to refence SNP DB (DBSNP)"

    	exec """
        	java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
             	-nt $threads
             	-R $REF 
             	-I $input.bam 
             	--dbsnp $DBSNP 
             	-stand_call_conf 50.0 -stand_emit_conf 10.0 
             	-dcov 1600 
             	-l INFO 
             	-A AlleleBalance -A DepthOfCoverage -A FisherStrand 
             	-glm INDEL 
             	-log $LOG -o $output.vcf
    	"""
}

@filter("filter")
filterSNPs = {
    	// Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.

	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }

	requires GATK : "Must provide a path to GATK (GATK)"   
        requires LOG : "Must provide a location for log file (LOG)"

    	exec """
        	java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration 
             	-R $REF 
             	--filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' 
             	--filterName 'GATK_MINIMAL_FILTER'
             	--variant $input.vcf 
             	-log $LOG 
             	-o $output.vcf
    	"""
}

@filter("filter")
filterIndels = {
    	doc """
            Filter data using very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
            If you have 10 or more samples GATK also recommends the filter InbreedingCoeff < -0.8
        	"""
	var sample_dir : false
        if (branch.sample_dir) { sample_dir = true }

        if (branch.outdir) {
                output.dir = branch.outdir + "/variants"
        } else {
                output.dir="variants"
        }

	requires GATK : "Must provide a path to GATK (GATK)"   
        requires REF : "Must provide a path to reference (REF)"
        requires LOG : "Must provide a location for log file (LOG)"

    	exec """
        	java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration 
                    -R $REF 
                    --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0' 
                    --filterName 'GATK_MINIMAL_FILTER' -log $LOG 
                    --variant $input.vcf 
                    -o $output.vcf
    	"""
}
