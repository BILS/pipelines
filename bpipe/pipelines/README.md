#		pipeline collection
#		###################

List of available pipelines:

1. [AnnotationAnnie](#annotationannie)
2. [AnnotationBlast2Go](#annotationblast2go)
3. [AnnotationGffStats](#annotationgffstats)
4. [AnnotationModels2Training](#annotationmodels2training)
5. [AnnotationPreprocessing](#annotationpreprocessing)
6. [AnnotationTranscriptMapAssemble](#annotationtranscriptmapAssemble)
7. [AnnotationTranscriptMapAssemble_fq](#annotationtranscriptmapassemble_fq)
8. [RNAseqAlignQuantify](#rnaseqalignquantify)
9. [RNAseqBamQualityControl](#rnaseqbamqualitycontrol)
10. [RNAseqQualityControl](#rnaseqqualitycontrol)

###	AnnotationAnnie

Title:		"A pipeline to execute blastp and interproscan searches to create input for the 'Annie' functional annotation pipeline"</br>
Input(s):	gff : Requires an annotation in GFF format as input (.gff)</br>
Dependencies: emboss, cufflinks, bioperl, blast, exonerate

----------------

### AnnotationBlast2Go

Title:		"A pipeline to execute blastp and interproscan searches and pass them to Blast2Go"</br>
Input(s):	gff : Requires an annotation in GFF format as input (.gff)
Dependencies: interproscan, blastp
----------------

###	AnnotationGffStats

Title:		"Pipeline to perform post-processing/statistical evaluation of GFF3-formatted annotation files"</br>
Input(s):	gff : Genome annotation file in GFF3 format
Dependencies:

----------------

###	AnnotationModels2Training

Title:		"Takes a genome anntation in GFF3 format and extracts data for training augustus profile"</br>
Input(s):	gff : A gene annotation file in GFF format
Dependencies:
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
