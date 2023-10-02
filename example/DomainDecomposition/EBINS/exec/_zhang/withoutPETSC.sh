#!/bin/bash
#SBATCH --qos=regular
#SBATCH -C cpu
#SBATCH -N 1
#SBATCH --time=01:00:00
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --account=m1411

./main.exe ./_2023/case1.inputs

