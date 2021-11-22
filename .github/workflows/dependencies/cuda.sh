#!/usr/bin/bash

set -eu -o pipefail

sudo apt-get -qq update
sudo apt-get install -y cmake g++ wget libblas-dev liblapacke-dev

sudo wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub
sudo apt-key add 7fa2af80.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64 /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list

sudo apt-get -qq update
sudo apt-get install -y cuda-compiler-11-0 cuda-nvtx-11-0 cuda-runtime-11-0 libcusolver-dev-11-0

sudo ln -s cuda-11.0 /usr/local/cuda
