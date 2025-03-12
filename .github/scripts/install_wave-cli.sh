#!/bin/bash

set -eo pipefail

wget https://github.com/seqeralabs/wave-cli/releases/download/v1.4.1/wave-1.4.1-linux-x86_64

mv wave-1.4.1-linux-x86_64 wave

chmod +x wave

mkdir -p /opt/wave/bin

mv wave /opt/wave/bin

echo "/opt/wave/bin" >> $GITHUB_PATH
