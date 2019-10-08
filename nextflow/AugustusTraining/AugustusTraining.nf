#! /usr/bin/env nextflow

// nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "$baseDir/test_data/test.gff"
params.genome = "$baseDir/test_data/test.gff"
params.outdir = "results"

params.gff_gene_model_filter_options = '-c -r -d 500 -a 0.3'

params.codon_table = 1

params.test_size = 100
params.flank_region_size = 500

log.info """
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Augustus training dataset workflow
 ===================================

 General Parameters
     gff_annotation                : ${params.gff_annotation}
     genome                        : ${params.genome}
     outdir                        : ${params.outdir}

 Gene Model Filtering parameters
     gff_gene_model_filter_options : ${params.gff_gene_model_filter_options}

 Protein Sequence extraction parameters
     codon_table                   : ${params.codon_table}

 Augustus training parameters
     test_size                     : ${params.test_size}
     flank_region_size             : ${params.flank_region_size}

 """
//
// include './../workflows/annotation_workflows' params(params)
//
// workflow {
//
// 	main:
// 	augustus_training_dataset(Channel.fromPath(gff_annotation, checkIfExists: true))
//
// 	publish:
// 	gbk2augustus.out.dataset to: "${params.outdir}/augustus_training_dataset"
//
// }
//
// workflow augustus_training_dataset {
//
// 	get:
// 		gff_annotation
//
// 	main:
// 		gff_filter_gene_models(gff_annotation)
// 		gff_longest_cds(gff_filter_gene_models.out)
// 		gff2protein(gff_longest_cds.out)
// 		blast_makeblastdb(gff2protein.out)
// 		blast_recursive(gff2protein.out,blast_makeblastdb.out)
// 		gff_filter_by_blast(gff_annotation,blast_recursive.out)
// 		gff2gbk(gff_filter_by_blast.out)
// 		gbk2augustus(gff2gbk.out)
//
// 	emit:
// 		dataset = gbk2augustus.out
//
// }

Channel.fromPath(params.gff_annotation, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find gff file matching ${params.gff_annotation}!\n" }
    .set { gff_for_gene_model_filter }
Channel.fromPath(params.genome, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
    .into { genome_for_gene_model; genome_for_gff2protein; genome_for_gff2gbk }

process gff_filter_gene_models {

    tag "${gff3_file.baseName}"

    input:
    file gff3_file from gff_for_gene_model_filter
    file genome_fasta from genome_for_gene_model.collect()

    output:
    file "${gff3_file.baseName}_model-filtered.gff3" into gff_for_longest_cds

    script:
    """
    filter_sort.pl -f $gff3_file -F $genome_fasta \\
        -o ${gff3_file.baseName}_model-filtered.gff3 ${params.gff_gene_model_filter_options}
    """
    // filter_sort.pl is a script in the NBIS
}

process gff_longest_cds {

    tag "${gff3_file.baseName}"

    input:
    file gff3_file from gff_for_longest_cds

    output:
    file "codingGeneFeatures.filter.longest_cds.gff" into gff_for_gff2protein, gff_for_blast_filter

    script:
    """
    gff3_sp_keep_longest_isoform.pl -f $gff3_file \\
        -o codingGeneFeatures.filter.longest_cds.gff
    """
    // gff3_sp_keep_longest_isoform.pl is a script in the NBIS GAAS repo
}

process gff2protein {

    tag "${gff_file.baseName}"

    input:
    file gff_file from gff_for_gff2protein
    file genome_fasta from genome_for_gff2protein.collect()

    output:
    file "${gff_file.baseName}_proteins.fasta" into fasta_for_blast, fasta_for_blastdb

    script:
    """
    gff3_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
        -p -cfs -cis -ct ${params.codon_table} --gff $gff_file
    """
    // gff3_sp_extract_sequences.pl is a script in the NBIS pipelines repository in bin

}

process blast_makeblastdb {

    tag "${fasta_file.baseName} type: $dbtype"

    input:
    file fasta_file from fasta_for_blastdb

    output:
    file "*.{phr,pin,psq}" into blastdb_files

    script:
    dbtype = "${params.dbtype}" ? "${params.dbtype}" : 'prot'
    """
    makeblastdb -in $fasta_file -dbtype $dbtype
    """

}

process blast_recursive {

    tag "${fasta_file.baseName}"

    input:
    file fasta_file from fasta_for_blast
    file blastdb from blastdb_files.collect()

    output:
    file "${fasta_file.baseName}_blast.tsv" into blast_tsv

    script:
    database = blastdb[0].toString() - ~/.p\w\w$/
    """
    blastp -query $fasta_file -db ${database} -num_threads ${task.cpus} \\
        -outfmt 6 -out ${fasta_file.baseName}_blast.tsv
    """

}

process gff_filter_by_blast {

    tag "${gff_file.baseName}"
    publishDir "${params.outdir}/BlastFilteredGFF", mode: 'copy'

    input:
    file gff_file from gff_for_blast_filter
    file blast_file from blast_tsv

    output:
    file "${gff_file.baseName}_blast-filtered.gff3" into blast_filtered_gff

    script:
    """
    gff_filter_by_mrna_id.pl --gff $gff_file --blast $blast_file \\
        --outfile ${gff_file.baseName}_blast-filtered.gff3
    """
    // gff_filter_by_mrna_id.pl is a script in the NBIS pipelines repository in bin

}

process gff2gbk {

    tag "${gff_file.baseName}"

    input:
    file gff_file from blast_filtered_gff
    file genome_fasta from genome_for_gff2gbk

    output:
    file "${gff_file.baseName}.gbk" into genbank_files

    script:
    """
    gff2gbSmallDNA.pl $gff_file $genome_fasta ${params.flank_region_size} ${gff_file.baseName}.gbk
    """
    // gff2gbSmallDNA.pl is a script in the Augustus package

}

process gbk2augustus {

    tag "Make Augustus training set: ${genbank_file.baseName}"
    publishDir "${params.outdir}/Augustus", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".train") > 0)        "TrainingData/$filename"
            else if (filename.indexOf(".test") > 0)    "TestingData/$filename"
            else if (filename.indexOf(".gbk") > 0)     "GenbankFile/$filename"
            else filename }

    input:
    file genbank_file from genbank_files

    output:
    file "${genbank_file.baseName}.train"
    file "${genbank_file.baseName}.test"
    file "${genbank_file.baseName}"

    script:
    """
    randomSplit.pl $genbank_file ${params.test_size}
    """
    // randomSplit.pl is a script in the Augustus package

}

workflow.onComplete {
    log.info ( workflow.success ? "\nAugustus training dataset complete!\n" : "Oops .. something went wrong\n" )
}
