#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/taxon_abundance.nf'
include { kraken2 } from './modules/taxon_abundance.nf'
include { bracken } from './modules/taxon_abundance.nf'
include { abundance_top_5 } from './modules/taxon_abundance.nf'
include { extract_reads } from './modules/taxon_abundance.nf'


workflow {

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  ch_kraken_db = Channel.fromPath( "${params.kraken_db}", type: 'dir')
  ch_bracken_db = Channel.fromPath( "${params.bracken_db}", type: 'dir')

  main:
    fastp(ch_fastq)
    kraken2(fastp.out.reads.combine(ch_kraken_db))
    bracken(kraken2.out.combine(ch_bracken_db))
    abundance_top_5(bracken.out)

    if (params.extract_reads) {
      ch_to_extract = bracken.out.map{ it -> it[1] }.splitCsv(header: true).filter{ it -> Float.parseFloat(it['fraction_total_reads']) > params.extract_reads_threshold }.map{ it -> [it['sample_id'], it['taxonomy_id']] }
      extract_reads(ch_fastq.join(kraken2.out).combine(ch_to_extract, by: 0))
    }
}
