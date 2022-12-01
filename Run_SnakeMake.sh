#!/bin/bash

SnakeFile=$1
# SnakeFile=~/codes/scRNAseq_Benchmark/Snakemake/Interdata.Snakefile


mkdir bu && mv BMoutput/HieRFIT/*.csv BMoutput/evaluation/Summary/HieRFIT.csv bu/;
SnakeFile=~/codes/scRNAseq_Benchmark/Snakemake/Snakefile

cores=12
. /n/home13/yasinkaymaz/miniconda3/etc/profile.d/conda.sh
#conda env export -n R3.3 > R3.6.yml
conda activate R3.6

snakemake -k \
  -s $SnakeFile \
  --reason \
  --configfile config.yml \
  --use-singularity \
  --use-conda \
  -j ${cores} \
  --singularity-args '--bind /n/home13/yasinkaymaz/LabSpace/data --bind /n/home13/yasinkaymaz/LabSpace/results/' -np
