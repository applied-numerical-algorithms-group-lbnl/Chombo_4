#!/usr/bin/bash

set -eu -o pipefail

sudo apt-get -qq update
sudo apt-get install -y cmake g++ wget

sudo wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub
sudo apt-key add 7fa2af80.pub

sudo apt-get -qq update
sudo apt-get install -y cuda-compiler-11-2

sudo ln -s cuda-11.2 /usr/local/cuda
