#		pipeline collection
#		###################

List of available pipelines:

1. [AnnotationAnnie](#annotationAnnie)
2. AnnotationBlast2Go
3. AnnotationGffStats
4. AnnotationModels2Training
5. AnnotationPreprocessing
6. AnnotationTranscriptMapAssemble
7. AnnotationTranscriptMapAssemble_fq
8. RNAseqAlignQuantify
9. RNAseqBamQualityControl
10. [RNAseqQualityControl](###rNAseqQualityControl)

###	AnnotationAnnie

Title:		"A pipeline to execute blastp and interproscan searches to create input for the 'Annie' functional annotation pipeline"</br>
Input(s):	gff : Requires an annotation in GFF format as input (.gff)

----------------

### AnnotationBlast2Go

Title:		"A pipeline to execute blastp and interproscan searches and pass them to Blast2Go"</br>
Input(s):	gff : Requires an annotation in GFF format as input (.gff)

----------------

###	AnnotationGffStats

Title:		"Pipeline to perform post-processing/statistical evaluation of GFF3-formatted annotation files"</br>
Input(s):	gff : Genome annotation file in GFF3 format

----------------

###	AnnotationModels2Training

Title:		"Takes a genome anntation in GFF3 format and extracts data for training augustus profile"</br>
Input(s):	gff : A gene annotation file in GFF format

----------------

###	AnnotationPreprocessing

Title:		"A pipeline to generate all pre-annotation production output from a genome sequence"</br>
Input(s):	fa : Requires genome sequence in fasta format

----------------

###	AnnotationTranscriptMapAssemble

Title:		"A pipeline to assemble transcripts from RNAseq data based on Cufflinks/Tophat"</br>
Input(s):	fq.gz : Requires compressed FastQ file(s) as input

----------------

###	AnnotationTranscriptMapAssemble_fq

Title:		"A pipeline to assemble transcripts from RNAseq data based on Cufflinks/Tophat"</br>
Input(s):	fq : Requires compressed FastQ file(s) as input

----------------

###	RNAseqAlignQuantify

Title:		"RNA-seq pipeline to align reads against a reference, clean the alignment and quantify using cufflinks"</br>
Input(s):	fq.gz : RNA-seq reads in gzipped fastq format, paired-end data expected in the format %_*.fq.gz.

----------------

###	RNAseqBamQualityControl
Title:		"Please add a tile"</br>
Input(s):	fq.gz : Requires read files in zipped FastQ format (fq.gz)

----------------

###	RNAseqQualityControl

Title:		"RNA-seq pipeline to perform quality control and trimming on RNA-seq read files"</br>
Input(s):	fq.gz : RNA-seq reads in gzipped fastq format, paired-end data expected in the format %_*.fq.gz.

----------------
