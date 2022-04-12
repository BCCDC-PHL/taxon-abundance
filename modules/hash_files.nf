process hash_files {

    tag { sample_id + " / " + file_type }

    input:
    tuple  val(sample_id), path(files_to_hash), val(file_type)

    output:
    tuple  val(sample_id), path("${sample_id}_${file_type}.sha256.csv"), emit: csv
    tuple  val(sample_id), path("${sample_id}_${file_type}_provenance.yml"), emit: provenance

    script:
    """
    shasum -a 256 ${files_to_hash} | tr -s ' ' ',' > ${sample_id}_${file_type}.sha256.csv
    while IFS=',' read -r hash filename; do
      printf -- "- input_filename: \$filename\\n  sha256: \$hash\\n" >> ${sample_id}_${file_type}_provenance.yml
    done < ${sample_id}_${file_type}.sha256.csv
    """
}
