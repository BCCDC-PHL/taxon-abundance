#!/usr/bin/env python

import argparse
import glob
import json
import os
import subprocess


def pull_image(image_url, destination_image_file):
    """
    """
    apptainer_pull_cmd = [
        "apptainer",
        "pull",
        destination_image_file,
        image_url,
    ]
    subprocess.run(apptainer_pull_cmd)


def push_image(source_image_file, image_url):
    """
    """

    apptainer_push_cmd = [
        "apptainer",
        "push",
        source_image_file,
        image_url,
    ]
    subprocess.run(apptainer_push_cmd)


def main(args):
    repo_owner = os.environ['GITHUB_REPOSITORY_OWNER'].lower()
    
    wave_jsons = glob.glob(os.path.join(args.wave_jsons_dir, "*.json"))
    for wave_json in wave_jsons:
        with open(wave_json, 'r') as f:
            w = json.load(f)
            pull_image_url = w['containerImage']
            image_name_with_version = pull_image_url.split('/')[-1]
            image_name, image_version = image_name_with_version.split(':')
            pull_destination = os.path.join(args.images_dir, f"{image_name}--{image_version}.img")
            pull_image(pull_image_url, pull_destination)

            push_image_url = f"oras://ghcr.io/{repo_owner}/{image_name}:{image_version}"
            push_image(pull_destination, push_image_url)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--wave-jsons-dir')
    parser.add_argument('--images-dir')
    args = parser.parse_args()
    main(args)
