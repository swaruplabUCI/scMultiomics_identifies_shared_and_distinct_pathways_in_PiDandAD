#!bin/bash

# directory containing cellranger outputs
cellranger_dir=$1

# output
outdir=$2

mkdir $outdir

# copy web reports:
for dir in $cellranger_dir/*
do
  name=$(basename $dir)
  cp $dir"/outs/web_summary.html" $outdir/$name"_web_summary.html"
done
