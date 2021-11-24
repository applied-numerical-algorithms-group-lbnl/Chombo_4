#!/usr/bin/bash

set -eu -o pipefail
sudo apt-get update
sudo apt-get install -y build-essential cmake g++ liblapacke-dev
sudo apt-get clean
