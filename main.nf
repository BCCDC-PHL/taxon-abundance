#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files } from './modules/hash_files.nf'
include { fastp } from './modules/taxon_abundance.nf'
include { kraken2 } from './modules/taxon_abundance.nf'
include { bracken } from './modules/taxon_abundance.nf'
include { abundance_top_5 } from './modules/taxon_abundance.nf'
include { extract_reads } from './modules/taxon_abundance.nf'
include { pipeline_provenance } from './modules/provenance.nf'
include { collect_provenance } from './modules/provenance.nf'


workflow {

  ch_start_time = Channel.of(LocalDateTime.now())
  ch_pipeline_name = Channel.of(workflow.manifest.name)
  ch_pipeline_version = Channel.of(workflow.manifest.version)

  ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  ch_kraken_db = Channel.fromPath( "${params.kraken_db}", type: 'dir')
  ch_bracken_db = Channel.fromPath( "${params.bracken_db}", type: 'dir')

  main:
    hash_files(ch_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq-input")))
    fastp(ch_fastq)
    kraken2(fastp.out.reads.combine(ch_kraken_db))
    bracken(kraken2.out.report.combine(ch_bracken_db))
    abundance_top_5(bracken.out.abundances)

    if (params.extract_reads) {
      ch_to_extract = bracken.out.abundances.map{ it -> it[1] }.splitCsv(header: true).filter{ it -> Float.parseFloat(it['fraction_total_reads']) > params.extract_reads_threshold }.map{ it -> [it['sample_id'], it['taxonomy_id']] }
      extract_reads(ch_fastq.join(kraken2.out.report).combine(ch_to_extract, by: 0))
    }

    ch_provenance = fastp.out.provenance
    ch_provenance = ch_provenance.join(kraken2.out.provenance).map{ it -> [it[0], [it[1], it[2]]] }
    ch_provenance = ch_provenance.join(bracken.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
