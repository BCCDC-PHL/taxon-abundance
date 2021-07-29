process fastp {

  tag { sample_id }

  input:
  tuple val(grouping_key), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
  tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads

  script:
  if (grouping_key =~ '_S[0-9]+_') {
    sample_id = grouping_key.split("_S[0-9]+_")[0]
  } else if (grouping_key =~ '_') {
    sample_id = grouping_key.split("_")[0]
  } else {
    sample_id = grouping_key
  }
  """
  fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz
  mv fastp.json ${sample_id}_fastp.json
  """
}

process fastp_json_to_csv {

  tag { sample_id }

  executor 'local'

  input:
  tuple val(sample_id), path(fastp_json)

  output:
  tuple val(sample_id), path("${sample_id}_read_count.csv")

  script:
  """
  fastp_json_to_csv.py -s ${sample_id} ${fastp_json} > ${sample_id}_read_count.csv
  """
}

process kraken2 {

  tag { sample_id }

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
  head -n 1 ${sample_id}_${params.taxonomic_level}_bracken_abundances_unsorted.tsv > bracken_abundances_header.tsv
  tail -n+2 ${sample_id}_${params.taxonomic_level}_bracken_abundances_unsorted.tsv | sort -t \$'\\t' -nrk 7,7 > ${sample_id}_${params.taxonomic_level}_bracken_abundances_data.tsv
  cat bracken_abundances_header.tsv ${sample_id}_${params.taxonomic_level}_bracken_abundances_data.tsv | sed 's/\\t/,/g' > ${sample_id}_${params.taxonomic_level}_bracken_abundances.csv
  """
}

process abundance_top_5 {

  tag { sample_id + " / " + params.taxonomic_level }

  executor 'local'

  cpus 1

  input:
  tuple val(sample_id), path(bracken_abundances)

  output:
  tuple val(sample_id), path("${sample_id}_${params.taxonomic_level}_top_5.csv")

  script:
  """
  bracken_top_n_linelist.py ${bracken_abundances} -n 5 -s ${sample_id} > ${sample_id}_${params.taxonomic_level}_top_5.csv
  """
}