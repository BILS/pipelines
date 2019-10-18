#! /usr/bin/env nextflow

// nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.maker_evidence_gff = "$baseDir/test_data/test.gff"
params.genome = "$baseDir/test_data/genome.fasta"
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
     maker_evidence_gff            : ${params.maker_evidence_gff}
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

Channel.fromPath(params.maker_evidence_gff, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find gff file matching ${params.maker_evidence_gff}!\n" }
    .set { gff_for_split_maker_evidence }
Channel.fromPath(params.genome, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
    .into { genome_for_gene_model; genome_for_gff2protein; genome_for_gff2gbk }

process split_maker_evidence {

    tag "${maker_evidence.baseName}"
    publishDir "${params.outdir}/maker_results_noAbinitio_clean", mode: 'copy'
    label 'GAAS'

    input:
    file maker_evidence from gff_for_split_maker_evidence

    output:
    file "maker_results_noAbinitio_clean/mrna.gff" into gff_for_model_select_by_AED
    file "maker_results_noAbinitio_clean/*"

    script:
    """
    gff3_sp_splitByLevel2Feature.pl -g ${maker_evidence} -o maker_results_noAbinitio_clean
    """
    // gff3_sp_splitByLevel2Feature.pl is a script in the NBIS GAAS repository
}

process model_selection_by_AED {

    tag "${mrna_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'GAAS'

    input:
    file mrna_gff from gff_for_model_select_by_AED

    output:
    file "codingGeneFeatures.filter.gff" into gff_for_longest_isoform

    script:
    """
    maker_select_models_by_AED_score.pl -f ${mrna_gff} -v 0.3 -t "<=" -o codingGeneFeatures.filter.gff
    """
    // maker_select_models_by_AED_score.pl is a script in the NBIS GAAS repository
}

process retain_longest_isoform {

    tag "${coding_gene_features_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'GAAS'

    input:
    file coding_gene_features_gff from gff_for_longest_isoform

    output:
    file "codingGeneFeatures.filter.longest_cds.gff" into gff_for_incomplete_gene_model_removal

    script:
    """
    gff3_sp_keep_longest_isoform.pl -f ${coding_gene_features_gff} -o codingGeneFeatures.filter.longest_cds.gff
    """
    // gff3_sp_keep_longest_isoform.pl is a script in the NBIS GAAS respository
}

process remove_incomplete_gene_models {

    tag "${coding_gene_features_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'GAAS'

    input:
    file coding_gene_features_gff from gff_for_incomplete_gene_model_removal
    file genome_fasta from genome_for_gene_model.collect()

    output:
    file "codingGeneFeatures.filter.longest_cds.complete.gff" into gff_complete_gene_models

    script:
    """
    gff3_sp_filter_incomplete_gene_coding_models.pl --gff ${coding_gene_features_gff} \
        -f ${genome_fasta} -o codingGeneFeatures.filter.longest_cds.complete.gff
    """
    // gff3_sp_filter_incomplete_gene_coding_models.pl is a script in the NBIS GAAS repository
}

process filter_by_locus_distance {

    tag "${coding_gene_features_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'GAAS'

    input:
    file coding_gene_features_gff from gff_complete_gene_models

    output:
    file "codingGeneFeatures.filter.longest_cds.complete.good_distance.gff" into gff_for_protein_extraction, gff_for_blast_filter

    script:
    """
    gff3_sp_filter_by_locus_distance.pl --gff ${coding_gene_features_gff} -o codingGeneFeatures.filter.longest_cds.complete.good_distance.gff
    """
    // gff3_sp_filter_by_locus_distance.pl is a script in the NBIS GAAS repository
}

// process gff_filter_gene_models {
//
//     tag "${gff3_file.baseName}"
//
//     input:
//     file gff3_file from gff_for_gene_model_filter
//     file genome_fasta from genome_for_gene_model.collect()
//
//     output:
//     file "${gff3_file.baseName}_model-filtered.gff3" into gff_for_longest_cds
//
//     script:
//     """
//     filter_sort.pl -f $gff3_file -F $genome_fasta \\
//         -o ${gff3_file.baseName}_model-filtered.gff3 ${params.gff_gene_model_filter_options}
//     """
//     // filter_sort.pl is a script in the NBIS
// }

// process gff_longest_cds {
//
//     tag "${gff3_file.baseName}"
//
//     input:
//     file gff3_file from gff_for_longest_cds
//
//     output:
//     file "codingGeneFeatures.filter.longest_cds.gff" into gff_for_gff2protein, gff_for_blast_filter
//
//     script:
//     """
//     gff3_sp_keep_longest_isoform.pl -f $gff3_file \\
//         -o codingGeneFeatures.filter.longest_cds.gff
//     """
//     // gff3_sp_keep_longest_isoform.pl is a script in the NBIS GAAS repo
// }

process extract_protein_sequence {

    tag "${gff_file.baseName}"
    label 'GAAS'

    input:
    file gff_file from gff_for_protein_extraction
    file genome_fasta from genome_for_gff2protein.collect()

    output:
    file "${gff_file.baseName}_proteins.fasta" into fasta_for_blast, fasta_for_blastdb

    script:
    """
    gff3_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
        -p -cfs -cis -ct ${params.codon_table} --g $gff_file
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
    """
    makeblastdb -in $fasta_file -dbtype prot
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
    label 'GAAS'

    input:
    file gff_file from gff_for_blast_filter
    file blast_file from blast_tsv.collect()

    output:
    file "${gff_file.baseName}_blast-filtered.gff3" into blast_filtered_gff

    script:
    """
    gff3_sp_filter_by_mrnaBlastValue_bioperl.pl --gff $gff_file --blast $blast_file \\
        --outfile ${gff_file.baseName}_blast-filtered.gff3
    """
    // gff_filter_by_mrna_id.pl is a script in the NBIS pipelines repository in bin

}

process gff2gbk {

    tag "${gff_file.baseName}"

    input:
    file gff_file from blast_filtered_gff
    file genome_fasta from genome_for_gff2gbk.collect()

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
