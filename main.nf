#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/taxon_abundance.nf'
include { fastp_json_to_csv } from './modules/taxon_abundance.nf'
include { kraken2 } from './modules/taxon_abundance.nf'
include { count_unclassified_reads } from './modules/taxon_abundance.nf'
include { bracken } from './modules/taxon_abundance.nf'
include { abundance_top_5 } from './modules/taxon_abundance.nf'

workflow {
  ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  ch_kraken_db = Channel.fromPath( "${params.kraken_db}", type: 'dir')
  ch_bracken_db = Channel.fromPath( "${params.bracken_db}", type: 'dir')

  main:
  fastp(ch_fastq)
  fastp_json_to_csv(fastp.out.fastp_json).map{ it -> it[1] }.collectFile(name:'read_qc.csv', keepHeader: true, sort: { it.text }, storeDir: "${params.outdir}")
  kraken2(fastp.out.reads.combine(ch_kraken_db))
  bracken(kraken2.out.combine(ch_bracken_db))
  abundance_top_5(bracken.out).map{ it -> it[1] }.collectFile(name:'abundances.csv', keepHeader: true, sort: { it.text }, storeDir: "${params.outdir}")
}
