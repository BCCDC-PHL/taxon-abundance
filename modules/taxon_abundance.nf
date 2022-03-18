process fastp {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_fastp.{json,csv}"

  input:
  tuple val(sample_id), path(reads_1), path(reads_2)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.{json,csv}"), emit: fastp_reports
  tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads

  script:
  """
  fastp -i ${reads_1} -I ${reads_2} -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz
  mv fastp.json ${sample_id}_fastp.json
  fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
  """
}

process kraken2 {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_kraken2.txt"

  input:
  tuple val(sample_id), path(reads_1), path(reads_2), path(kraken2_db)

  output:
  tuple val(sample_id), path("${sample_id}_kraken2.txt")

  script:
  """
  kraken2 --db ${kraken2_db} --threads ${task.cpus} --output "-" --report ${sample_id}_kraken2.txt --paired ${reads_1} ${reads_2}
  """
}

process bracken {

  tag { sample_id + " / " + params.taxonomic_level }
  
  errorStrategy 'ignore'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_bracken_abundances.csv"

  input:
  tuple val(sample_id), path(kraken2_report), path(bracken_db)

  output:
  tuple val(sample_id), path("${sample_id}_${params.taxonomic_level}_bracken_abundances.csv")

  script:
  """
  bracken -d ${bracken_db} \
    -i ${kraken2_report} \
    -w ${sample_id}_${params.taxonomic_level}_bracken.txt \
    -o ${sample_id}_${params.taxonomic_level}_bracken_abundances_unsorted.tsv \
    -r ${params.read_length} \
    -l ${params.taxonomic_level}
  head -n 1 ${sample_id}_${params.taxonomic_level}_bracken_abundances_unsorted.tsv | tr \$'\\t' ',' > bracken_abundances_header.csv
  adjust_bracken_percentages_for_unclassified_reads.py \
    -k ${kraken2_report} \
    -b ${sample_id}_${params.taxonomic_level}_bracken.txt \
    -a ${sample_id}_${params.taxonomic_level}_bracken_abundances_unsorted.tsv \
    > ${sample_id}_${params.taxonomic_level}_bracken_abundances_unsorted_with_unclassified.csv
  tail -n+2 ${sample_id}_${params.taxonomic_level}_bracken_abundances_unsorted_with_unclassified.csv | sort -t ',' -nrk 7,7 > ${sample_id}_${params.taxonomic_level}_bracken_abundances_data.csv
  cat bracken_abundances_header.csv ${sample_id}_${params.taxonomic_level}_bracken_abundances_data.csv > ${sample_id}_${params.taxonomic_level}_bracken_abundances.csv
  """
}

process abundance_top_5 {

  tag { sample_id + " / " + params.taxonomic_level }

  executor 'local'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_top_5.csv"

  input:
  tuple val(sample_id), path(bracken_abundances)

  output:
  tuple val(sample_id), path("${sample_id}_${params.taxonomic_level}_top_5.csv")

  script:
  """
  bracken_top_n_linelist.py ${bracken_abundances} -n 5 -s ${sample_id} > ${sample_id}_${params.taxonomic_level}_top_5.csv
  """
}