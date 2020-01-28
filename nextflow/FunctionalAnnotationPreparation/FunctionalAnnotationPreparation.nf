#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "/path/to/annotation.gff"
params.genome = "/path/to/genome.fasta"
params.outdir = "results"

params.codon_table = 1

params.records_per_file = 1000

params.blast_db_fasta = '/path/to/protein/database.fasta'

params.interproscan_db = 'all'

params.merge_annotation_identifier = 'ID'

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
     gff_annotation                 : ${params.gff_annotation}
     genome                         : ${params.genome}
     outdir                         : ${params.outdir}

 Parallelisation parameters
     records_per_file               : ${params.records_per_file}

 Gff2Protein parameters
     codon_table                    : ${params.codon_table}

 Blast parameters
     blast_db_fasta                 : ${params.blast_db_fasta}

 Interproscan parameters
     interproscan_db                : ${params.interproscan_db}

 Merge functional annotation parameters
     merge_annotation_identifier    : ${params.merge_annotation_identifier}

 """)

// include './../workflows/annotation_workflows' params(params)
//
workflow {

	main:
    annotation = Channel.fromPath(params.gff_annotation, checkIfExists: true)
        .ifEmpty { exit 1, "Cannot find gff file matching ${params.gff_annotation}!\n" }
    genome = Channel.fromPath(params.genome, checkIfExists: true)
        .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
    blastdb = Channel.fromPath("${params.blast_db_fasta}{,.p*}", checkIfExists: true)
        .ifEmpty { exit 1, "Cannot find blast database files matching ${params.blast_db_fasta}{,.p*}" }
	functional_annotation_input_preparation(annotation,genome,blastdb)

	// publish:
	// functional_annotation_input_preparation.out to: "${params.outdir}"
}

workflow functional_annotation_input_preparation {

	get:
		gff_file
		genome
        blastdb

	main:
		gff2protein(gff_file,genome.collect())
		blastp(gff2protein.out.splitFasta(by: params.records_per_file, file: true),blastdb.collect())
		interproscan(gff2protein.out.splitFasta(by: params.records_per_file, file: true))
        merge_functional_annotation(gff_file,
            blastp.out.collectFile(name:'blast_merged.tsv').collect(),
            interproscan.out.collectFile(name:'interproscan_merged.tsv').collect(),
            blastdb.collect()
            )

	// emit:
	// 	blast_results = merge_blast_tab.out
	// 	interpro_tsv = merge_interpro_tsv.out

}

// Channel.fromPath(params.gff_annotation, checkIfExists: true)
//     .ifEmpty { exit 1, "Cannot find gff file matching ${params.gff_annotation}!\n" }
//     .into { gff_for_gff2protein; gff_for_functional_merge }
// Channel.fromPath(params.genome, checkIfExists: true)
//     .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
//     .into { genome_for_gene_model; genome_for_gff2protein; genome_for_gff2gbk }
// Channel.fromPath("${params.blast_db_fasta}{,.p*}", checkIfExists: true)
//     .ifEmpty { exit 1, "Cannot find blast database files matching ${params.blast_db_fasta}{,.p*}" }
//     .into { blastdb_files; blastdb_files_for_gff_merge }

process gff2protein {

    // tag "${gff_file.baseName}"
    label 'AGAT'

    input:
    path gff_file // from gff_for_gff2protein
    path genome_fasta // from genome_for_gff2protein.collect()

    output:
    path "${gff_file.baseName}_proteins.fasta"  //into fasta_for_blast, fasta_for_interpro

    script:
    """
    agat_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
        -p -cfs -cis -ct ${params.codon_table} --gff $gff_file
    """
    // agat_sp_extract_sequences.pl is a script from AGAT

}

process blastp {

    // tag "$database"

    input:
    path fasta_file // from fasta_for_blast.splitFasta(by: params.records_per_file)
    path blastdb // from blastdb_files.collect()

    output:
    path "${fasta_file.baseName}_blast.tsv" //into blast_tsvs

    script:
    // database = blastdb[0].toString() - ~/.p\w\w$/
    """
    blastp -query $fasta_file -db ${blastdb} -num_threads ${task.cpus} \\
        -outfmt 6 -out ${fasta_file.baseName}_blast.tsv
    """

}

process interproscan {

    // tag "InterProScan: Protein function classification"

    input:
    path protein_fasta // from fasta_for_interpro.splitFasta(by: params.records_per_file)

    output:
    // file '*.gff3' into interpro_gffs
    // file 'results/*.xml' into interpro_xmls
    path '*.tsv' // into interpro_tsvs

    script:
    applications = { params.interproscan_db ? "-appl ${params.interproscan_db}" : '' }
    """
    interproscan.sh ${applications} -i $protein_fasta -o ${protein_fasta.baseName}.tsv \\
        -f TSV --iprlookup --goterms -pa -dp -t p
    """

}

process merge_functional_annotation {

    publishDir "${params.outdir}/blast_tsv", mode: 'copy', pattern: 'blast_merged.tsv'
    publishDir "${params.outdir}/interproscan_tsv", mode: 'copy', pattern: 'interproscan_merged.tsv'
    publishDir "${params.outdir}/final_annotation", mode: 'copy', pattern: "${gff_annotation.baseName}_plus-functional-annotation.gff"
    label 'AGAT'

    input:
    path gff_annotation //from gff_for_functional_merge
    path merged_blast_results //from blast_tsvs.collectFile(name:'blast_merged.tsv').collect()
    path merged_interproscan_results //from interpro_tsvs.collectFile(name:'interproscan_merged.tsv').collect()
    path blast_files //from blastdb_files_for_gff_merge.collect()

    output:
    path "${gff_annotation.baseName}_plus-functional-annotation.gff"

    script:
    """
    agat_sp_manage_functional_annotation.pl -f ${gff_annotation} \\
        -b ${merged_blast_results} -i ${merged_interproscan_results} \\
        -db ${params.blast_db_fasta} -id ${params.merge_annotation_identifier} \\
        -o ${gff_annotation.baseName}_plus-functional-annotation.gff
    """
    // agat_sp_manage_functional_annotation.pl is a script from AGAT

}

workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}
