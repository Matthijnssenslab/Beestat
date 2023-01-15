# Beestat

Accompanying repository for association between viral loads and bee colony health in Belgium.  
All source data can be found under the 'data' folder.
To generate the figures, first create the conda environment:

  > conda create -f conda.yaml -n beestat

and run

  > conda activate beestat
  > snakemake --cores 1 -s generate_figures.smk

