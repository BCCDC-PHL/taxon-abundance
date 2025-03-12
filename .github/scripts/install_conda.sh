#!/bin/bash
set -eo pipefail

artifacts_dir="artifacts"

echo "Install Miniconda .." >> ${artifacts_dir}/test.log

export PATH=/opt/miniconda3/bin:$PATH

wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh

/bin/bash ~/miniconda.sh -b -p /opt/miniconda3

rm ~/miniconda.sh

echo ". /opt/minconda3/etc/profile.d/conda.sh" >> ~/.bashrc

conda update -n base -c defaults conda

conda install -y -c conda-forge mamba

conda init bash
