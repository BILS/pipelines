#! /usr/bin/env nextflow

// nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "$baseDir/test_data/test.gff"
params.genome = "$baseDir/test_data/genome.fasta"
params.outdir = "results"

params.records_per_file = 1000

params.blast_db_fasta = '/path/to/protein/database.fasta'

params.interproscan_db = 'all'

log.info("""
NBIS
 _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Functional annotation input preparation workflow
 ===================================

 General parameters
     gff_annotation    : ${params.gff_annotation}
     genome            : ${params.genome}
     outdir            : ${params.outdir}


 Parallelisation parameters
     records_per_file  : ${params.records_per_file}

 Blast parameters
     blast_db_fasta    : ${params.blast_db_fasta}

 Interproscan parameters
     interproscan_db   : ${params.interproscan_db}

 """)

// include './../workflows/annotation_workflows' params(params)
//
// workflow {
//
// 	main:
// 	functional_annotation_input_preparation(Channel.fromPath(params.gff_annotation, checkIfExists: true))
//
// 	publish:
// 	functional_annotation_input_preparation.out to: "${params.outdir}"
// }

// workflow functional_annotation_input_preparation {
//
// 	get:
// 		gff_file
// 		genome
//
// 	main:
// 		gff2protein(gff_file,genome.collect())
// 		blastp(gff2protein.out.splitFasta(by: params.records_per_file))
// 		merge_blast_tab(blastp.out.collect())
// 		interpro(gff2protein.out.splitFasta(by: params.records_per_file))
// 		merge_interpro_tsv(interpro.out.collect())
//
// 	emit:
// 		blast_results = merge_blast_tab.out
// 		interpro_tsv = merge_interpro_tsv.out
//
// }

Channel.fromPath(params.gff_annotation, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find gff file matching ${params.gff_annotation}!\n" }
    .into { gff_for_gff2protein; gff_for_functional_merge }
Channel.fromPath(params.genome, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
    .into { genome_for_gene_model; genome_for_gff2protein; genome_for_gff2gbk }
Channel.fromPath("${params.blast_db}{,.p*}", checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find blast database files matching ${params.blast_db}{,.p*}" }
	.into { blastdb_files; blastdb_files_for_gff_merge }

process gff2protein {

    // tag "${gff_file.baseName}"
    label 'GAAS'

    input:
    file gff_file from gff_for_gff2protein
    file genome_fasta from genome_for_gff2protein.collect()

    output:
    file "${gff_file.baseName}_proteins.fasta" into fasta_for_blast, fasta_for_interpro

    script:
    """
    gff3_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
        -p -cfs -cis -ct ${params.codon_table} --gff $gff_file
    """
    // gff3_sp_extract_sequences.pl is a script in the NBIS pipelines repository in bin

}

process blastp {

    // tag "$database"

    input:
    file fasta_file from fasta_for_blast.splitFasta(by: params.records_per_file)
    file blastdb from blastdb_files.collect()

    output:
    file "${fasta_file.baseName}_blast.tsv" into blast_tsvs

    script:
    // database = blastdb[0].toString() - ~/.p\w\w$/
    """
    blastp -query $fasta_file -db ${params.blast_db_fasta} -num_threads ${task.cpus} \\
        -outfmt 6 -out ${fasta_file.baseName}_blast.tsv
    """

}

// process merge_blast_tab {
//
//     // tag "Merge: Blast TSVs"
//     publishDir "${params.outdir}/blast_tsv", mode: 'copy'
//
//     input:
//     file blast_fragments from blast_tsvs.collect()
//
//     output:
//     file 'merged_blast_results.tsv'
//
//     script:
//     """
//     cat $blast_fragments > merged_blast_results.tsv
//     """
//
// }

process interproscan {

    // tag "InterProScan: Protein function classification"

    input:
    file protein_fasta from fasta_for_interpro.splitFasta(by: params.records_per_file)

    output:
    // file '*.gff3' into interpro_gffs
    // file 'results/*.xml' into interpro_xmls
    file 'results/*.tsv' into interpro_tsvs

    script:
    """
    interproscan.sh -appl ${params.interproscan_db} -i $protein_fasta -d results \\
	    --iprlookup --goterms -pa -dp -t p
    """

}

// process merge_interpro_xml {
//
//     tag "Merge: InterProScan XMLs"
//     publishDir "${params.outdir}/interproscan_xml", mode: 'copy'
//
//     input:
//     file xml_files from interpro_xmls.collect()
//
//     output:
//     file 'interpro_search.xml'
//
//     // This code is not robust at all. Need to rewrite (e.g. the -v "xml" already excludes "protein-matches" lines because of the "xmlns" attributes)
//     script:
//     """
//     head -n 2 ${xml_files[0]} > interpro_search.xml
//     for XML in $xml_files; do
//     grep -v "xml" \$XML | grep -v "protein-matches" >> interpro_search.xml
//     done
//     tail -n 1 ${xml_files[0]} >> interpro_search.xml
//     """
// }

// process merge_interpro_tsv {
//
//     // tag "Merge InterProScan TSVs"
//     publishDir "${params.outdir}/interproscan_tsv", mode: 'copy'
//
//     input:
//     file tsv_files from interpro_tsvs.collect()
//
//     output:
//     file 'interpro_search.tsv'
//
//     script:
//     """
//     cat $tsv_files > interpro_search.tsv
//     """
//
// }

process merge_functional_annotation {

	publishDir "${params.outdir}/blast_tsv", mode: 'copy', pattern: 'blast_merged.tsv'
	publishDir "${params.outdir}/interproscan_tsv", mode: 'copy', pattern: 'interproscan_merged.tsv'
	publishDir "${params.outdir}/final_annotation", mode: 'copy', pattern: "${gff_annotation.baseName}_plus-functional-annotation.gff"
	label 'GAAS'

	input:
	file gff_annotation from gff_for_functional_merge
	file merged_blast_results from blast_tsvs.collectFile(name:'blast_merged.tsv').collect()
	file merged_interproscan_results from interpro_tsvs.collectFile(name:'interproscan_merged.tsv').collect()
	file blast_files from blastdb_files_for_gff_merge.collect()

	output:
	file "${gff_annotation.baseName}_plus-functional-annotation.gff"

	script:
	""""
	gff3_sp_manage_functional_annotation.pl -f ${gff_annotation} \\
	    -b ${merged_blast_results} -i ${merged_interproscan_results} \\
		-db ${params.blast_db_fasta} -id BALIM \\
		-o ${gff_annotation.baseName}_plus-functional-annotation.gff
	"""
	// gff3_sp_manage_functional_annotation.pl is a script in the NBIS GAAS repository

}

workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}
