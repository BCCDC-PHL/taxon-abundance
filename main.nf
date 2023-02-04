#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files }             from './modules/hash_files.nf'
include { fastp }                  from './modules/taxon_abundance.nf'
include { kraken2 }                from './modules/taxon_abundance.nf'
include { bracken }                from './modules/taxon_abundance.nf'
include { abundance_top_5 }        from './modules/taxon_abundance.nf'
include { abundance_top_5_kraken } from './modules/taxon_abundance.nf'
include { kraken_abundances }      from './modules/taxon_abundance.nf'
include { extract_reads }          from './modules/taxon_abundance.nf'
include { pipeline_provenance }    from './modules/provenance.nf'
include { collect_provenance }     from './modules/provenance.nf'
include { bracken_subspeciation }                from './modules/subspeciation_taxon_abundance.nf'
include { abundance_top_5_subspeciation }        from './modules/subspeciation_taxon_abundance.nf'

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

    if (!params.skip_bracken) {
      bracken(kraken2.out.report.combine(ch_bracken_db))
      abundance_top_5(bracken.out.abundances)
      ch_abundances = bracken.out.abundances
    } else {
      abundance_top_5_kraken(kraken2.out.report)
      ch_abundances = kraken_abundances(kraken2.out.report)
    }

    if (params.extract_reads) {
      ch_to_extract = ch_abundances.map{ it -> it[1] }.splitCsv(header: true).filter{ it -> Float.parseFloat(it['fraction_total_reads']) > params.extract_reads_threshold }.map{ it -> [it['sample_id'], it['taxonomy_id']] }
      extract_reads(ch_fastq.join(kraken2.out.report).combine(ch_to_extract, by: 0))
    }

    if (params.collect_outputs) {
      fastp.out.csv.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_fastp.csv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: { it -> it.readLines()[1] })
      if (!params.skip_bracken) {
        ch_abundances.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_" + params.taxonomic_level + "_bracken_abundances.csv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: { it -> it.readLines()[1].split(',')[0] })
        abundance_top_5.out.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_" + params.taxonomic_level + "_bracken_abundances_top_5.csv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: { it -> it.readLines()[1] })
      } else {
        ch_abundances.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_" + params.taxonomic_level + "_kraken_abundances.csv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: { it -> it.readLines()[1].split(',')[0] })
        abundance_top_5_kraken.out.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_" + params.taxonomic_level + "_kraken_abundances_top_5.csv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: { it -> it.readLines()[1] })
      }
    }

    if (params.subspeciation) {        
      bracken_subspeciation(kraken2.out.report.combine(ch_bracken_db))
      abundance_top_5_subspeciation(bracken_subspeciation.out.abundances)
      ch_abundances = bracken_subspeciation.out.abundances
    }

    ch_provenance = fastp.out.provenance

    ch_provenance = ch_provenance.join(kraken2.out.provenance).map{ it -> [it[0], [it[1], it[2]]] }

    if (!params.skip_bracken) {
      ch_provenance = ch_provenance.join(bracken.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }
    
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
