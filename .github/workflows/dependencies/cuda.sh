#!/usr/bin/bash

set -eu -o pipefail

sudo apt-get -qq update
sudo apt-get install -y cmake g++ wget libblas-dev liblapack-dev

sudo wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub
sudo apt-key add 7fa2af80.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64 /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list

sudo apt-get -qq update
sudo apt-get install -y cuda-compiler-11-2 cuda-command-line-tools-11-2 cuda-nvtx-11-2 cuda-runtime-11-2

sudo ln -s cuda-11.2 /usr/local/cuda
