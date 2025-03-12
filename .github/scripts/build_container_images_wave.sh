#!/bin/bash

mkdir -p artifacts
mkdir -p wave_images

for env_yaml in environments/*.yml; do
    image_name=$(head -n 1 $env_yaml | cut -d ' ' -f 2)
    echo "building image ${image_name} from file ${env_yaml}..."
    wave \
	--conda-file ${env_yaml} \
	--singularity \
	--freeze \
	--await \
	--output json \
	| python -m json.tool \
	| tee wave_images/${image_name}.json
    echo "done building image ${image_name}"
    cp wave_images/${image_name}.json artifacts/

done


