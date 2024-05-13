process collect_provenance {

    tag { sample_id }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_*_provenance.yml", mode: 'copy'

    input:
    tuple val(sample_id), path(provenance_files)

    output:
    tuple val(sample_id), file("${sample_id}_*_provenance.yml")

    script:
    """
    cat ${provenance_files} > ${sample_id}_\$(date +%Y%m%d%H%M%S)_provenance.yml
    """
}

process pipeline_provenance {

    tag { pipeline_name + " / " + pipeline_version }

    executor 'local'

    input:
    tuple val(session_id), val(run_name), val(pipeline_name), val(pipeline_version), val(timestamp_analysis_start)

    output:
    file("pipeline_provenance.yml")

    script:
    """
    printf -- "- pipeline_name: ${pipeline_name}\\n"       >> pipeline_provenance.yml
    printf -- "  pipeline_version: ${pipeline_version}\\n" >> pipeline_provenance.yml
    printf -- "  nextflow_session_id: ${session_id}\\n"    >> pipeline_provenance.yml
    printf -- "  nextflow_run_name: ${run_name}\\n"        >> pipeline_provenance.yml
    printf -- "  timestamp_analysis_start: ${timestamp_analysis_start}\\n" >> pipeline_provenance.yml
    """
}
