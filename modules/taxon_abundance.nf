process fastp {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_fastp.{json,csv}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: csv
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n"  >> ${sample_id}_fastp_provenance.yml
    printf -- "  tools:\\n"               >> ${sample_id}_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml
    
    fastp \
	-i ${reads_1} \
	-I ${reads_2} \
	-o ${sample_id}_trimmed_R1.fastq.gz \
	-O ${sample_id}_trimmed_R2.fastq.gz
    
    mv fastp.json ${sample_id}_fastp.json
    
    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
    """
}

process kraken2 {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_kraken2_report.txt"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(kraken2_db)

    output:
    tuple val(sample_id), path("${sample_id}_kraken2_output.tsv"), path("${sample_id}_kraken2_report.txt"), emit: report
    tuple val(sample_id), path("${sample_id}_kraken2_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: kraken2\\n"                                                      >> ${sample_id}_kraken2_provenance.yml
    printf -- "  tools:\\n"                                                                     >> ${sample_id}_kraken2_provenance.yml
    printf -- "    - tool_name: kraken2\\n"                                                     >> ${sample_id}_kraken2_provenance.yml
    printf -- "      tool_version: \$(kraken2 --version | grep 'version' | cut -d ' ' -f 3)\\n" >> ${sample_id}_kraken2_provenance.yml
    printf -- "      parameters:\\n"                                                            >> ${sample_id}_kraken2_provenance.yml
    printf -- "        - name: confidence\\n"                                                   >> ${sample_id}_kraken2_provenance.yml
    printf -- "          value: ${params.confidence}\\n"                                        >> ${sample_id}_kraken2_provenance.yml

    get_kraken_db_metadata.py --db ${kraken2_db} --provenance ${sample_id}_kraken2_provenance.yml

    kraken2 \
	--db ${kraken2_db} \
	--threads ${task.cpus} \
	--confidence ${params.confidence} \
	--output ${sample_id}_kraken2_output.tsv \
	--report ${sample_id}_kraken2_report.txt \
	--paired ${reads_1} ${reads_2}
  """
}

process bracken {

    tag { sample_id + " / " + taxonomic_level }
  
    errorStrategy 'ignore'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_bracken_abundances.csv"

    input:
    tuple val(sample_id), path(kraken2_output), path(kraken2_report), path(bracken_db), val(taxonomic_level)

    output:
    tuple val(sample_id), path("${sample_id}_${taxonomic_level}_bracken_abundances.csv"), val(taxonomic_level), emit: abundances
    tuple val(sample_id), path("${sample_id}_bracken_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: bracken\\n"    >> ${sample_id}_bracken_provenance.ym
    printf -- "  tools:\\n"                   >> ${sample_id}_bracken_provenance.yml
    printf -- "    - tool_name: bracken\\n"   >> ${sample_id}_bracken_provenance.yml
    printf -- "      tool_version: 2.6.1\\n"  >> ${sample_id}_bracken_provenance.yml
    printf -- "      parameters:\\n"          >> ${sample_id}_bracken_provenance.yml
    printf -- "        - name: read_length\\n" >> ${sample_id}_bracken_provenance.yml
    printf -- "          value: ${params.read_length}\\n" >> ${sample_id}_bracken_provenance.yml
    printf -- "        - name: taxonomic_level\\n"        >> ${sample_id}_bracken_provenance.yml
    printf -- "          value: ${taxonomic_level}\\n"    >> ${sample_id}_bracken_provenance.yml

    bracken -d ${bracken_db} \
	-i ${kraken2_report} \
	-w ${sample_id}_${taxonomic_level}_bracken.txt \
	-o ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted.tsv \
	-r ${params.read_length} \
	-l ${taxonomic_level}

    paste <(echo "sample_id") <(head -n 1 ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted.tsv) | tr \$'\\t' ',' > bracken_abundances_header.csv

    adjust_bracken_percentages_for_unclassified_reads.py \
	-k ${kraken2_report} \
	-b ${sample_id}_${taxonomic_level}_bracken.txt \
	-a ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted.tsv \
	> ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted_with_unclassified.csv

    tail -n+2 ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted_with_unclassified.csv | \
	sort -t ',' -nrk 7,7 | \
	awk -F ',' 'BEGIN {OFS=FS}; {print "${sample_id}",\$0}' > ${sample_id}_${taxonomic_level}_bracken_abundances_data.csv

    cat bracken_abundances_header.csv ${sample_id}_${taxonomic_level}_bracken_abundances_data.csv > ${sample_id}_${taxonomic_level}_bracken_abundances.csv
    """
}

process abundance_top_5 {

    tag { sample_id + " / " + taxonomic_level }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_top_5.csv"

    input:
    tuple val(sample_id), path(bracken_abundances), val(taxonomic_level)

    output:
    tuple val(sample_id), path("${sample_id}_${taxonomic_level}_top_5.csv"), val(taxonomic_level)

    script:
    """
    bracken_top_n_linelist.py ${bracken_abundances} -n 5 -s ${sample_id} > ${sample_id}_${taxonomic_level}_top_5.csv
    """
}

process abundance_top_5_kraken {

    tag { sample_id + " / " + taxonomic_level }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_top_5.csv"

    input:
    tuple val(sample_id), path(kraken_output), path(kraken_report), val(taxonomic_level)

    output:
    tuple val(sample_id), path("${sample_id}_${taxonomic_level}_top_5.csv"), val(taxonomic_level)

    script:
    """
    kraken_top_n_linelist.py ${kraken_report} -n 5 -s ${sample_id} --taxonomy-lvl ${taxonomic_level} > ${sample_id}_${taxonomic_level}_top_5.csv
    """
}


process kraken_abundances {

    tag { sample_id + " / " + taxonomic_level }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*_kraken2_abundances.csv"

    input:
    tuple val(sample_id), path(kraken_output), path(kraken_report), val(taxonomic_level)

    output:
    tuple val(sample_id), path("${sample_id}_${taxonomic_level}_kraken2_abundances.csv"), val(taxonomic_level)

    script:
    """
    kraken_abundances.py ${kraken_report} -s ${sample_id} --taxonomy-lvl ${taxonomic_level} > ${sample_id}_${taxonomic_level}_kraken2_abundances.csv
    """
}


process extract_reads {

    tag { sample_id + ' / ' + taxid }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output/extracted_reads_by_taxid" : "${params.outdir}/${sample_id}/extracted_reads_by_taxid", mode: 'copy', pattern: "${taxid}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(kraken2_output), path(kraken2_report), val(taxid)

    output:
    tuple val(sample_id), val(taxid), path("${taxid}", type: "dir")

    script:
    """
    mkdir ${taxid}

    extract_kraken_reads.py \
	-k ${kraken2_output} \
	-r ${kraken2_report} \
	-1 ${reads_1} \
	-2 ${reads_2} \
	--taxid ${taxid} \
	--fastq-output \
	-o ${taxid}/${sample_id}-taxid-${taxid}_R1.fastq \
	-o2 ${taxid}/${sample_id}-taxid-${taxid}_R2.fastq \
	--include-children

    gzip ${taxid}/*.fastq
    """
}
