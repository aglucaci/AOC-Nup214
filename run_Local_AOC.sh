#!/bin/bash

clear

echo "# Analysis of Orthologous Collections (AOC)."
echo "# @Author: Alexander G. Lucaci"
echo ""

# Set up the pipeline failure expectations.
set -euo pipefail

# Downloading hyphy-analyses.
FOLDER="hyphy-analyses"
URL="https://github.com/veg/hyphy-analyses.git"

if [ ! -d "$FOLDER" ] ; then
    echo "# Downloading hyphy-analyses repository"
    git clone "$URL" "$FOLDER"
fi

if [ ! -d "logs" ] ; then
    echo "# Creating 'logs' directory"
    mkdir -p logs
fi

#echo "Executing HPC Snakemake command"

echo "# Executing Snakemake command"

snakemake \
      -s Snakefile.smk.py \
      --jobs 1 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 120 \
      --conda-frontend conda

# End of file 

echo "# Done"
exit 0





