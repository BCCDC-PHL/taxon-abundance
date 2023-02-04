process bracken_subspeciation {

  tag { sample_id + " / " + params.subspecies_taxonomic_level }
  
  errorStrategy 'ignore'

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_bracken_abundances.csv"

  input:
  tuple val(sample_id), path(kraken2_output), path(kraken2_report), path(bracken_db)

  output:
  tuple val(sample_id), path("${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances.csv"), emit: abundances
  tuple val(sample_id), path("${sample_id}_bracken_provenance.yml"), emit: provenance

  script:
  """
  printf -- "- process_name: bracken\\n" > ${sample_id}_bracken_provenance.yml
  printf -- "  tool_name: bracken\\n  tool_version: 2.6.1\\n" >> ${sample_id}_bracken_provenance.yml
  printf -- "  database_path: \$(readlink -f ${bracken_db})\\n" >> ${sample_id}_bracken_provenance.yml
  bracken -d ${bracken_db} \
    -i ${kraken2_report} \
    -w ${sample_id}_${params.subspecies_taxonomic_level}_bracken.txt \
    -o ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances_unsorted.tsv \
    -r ${params.read_length} \
    -l ${params.subspecies_taxonomic_level}

  paste <(echo "sample_id") <(head -n 1 ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances_unsorted.tsv) | tr \$'\\t' ',' > bracken_abundances_header.csv

  adjust_bracken_percentages_for_unclassified_reads.py \
    -k ${kraken2_report} \
    -b ${sample_id}_${params.subspecies_taxonomic_level}_bracken.txt \
    -a ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances_unsorted.tsv \
    > ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances_unsorted_with_unclassified.csv

  tail -n+2 ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances_unsorted_with_unclassified.csv | \
    sort -t ',' -nrk 7,7 | \
    awk -F ',' 'BEGIN {OFS=FS}; {print "${sample_id}",\$0}' > ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances_data.csv

  cat bracken_abundances_header.csv ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances_data.csv > ${sample_id}_${params.subspecies_taxonomic_level}_bracken_abundances.csv
  """
}

process abundance_top_5_subspeciation {

  tag { sample_id + " / " + params.subspecies_taxonomic_level }

  executor 'local'

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_top_5.csv"

  input:
  tuple val(sample_id), path(bracken_abundances)

  output:
  tuple val(sample_id), path("${sample_id}_${params.subspecies_taxonomic_level}_top_5.csv")

  script:
  """
  bracken_top_n_linelist.py ${bracken_abundances} -n 5 -s ${sample_id} > ${sample_id}_${params.subspecies_taxonomic_level}_top_5.csv
  """
}

