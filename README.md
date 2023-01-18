# Beestat

Accompanying repository for association between viral loads and bee colony health in Belgium.

## generate figures

All source data can be found under the 'data' folder.
To generate the figures, first clone the repository and create the conda environment:

  > git clone git@github.com:Matthijnssenslab/Beestat.git  
  > cd Beestat  
  > conda env create -f conda.yaml -n beestat  

or if you have mamba

  > mamba env create -f conda.yaml -n beestat  

and run

  > conda activate beestat  
  > snakemake --cores 1 -s generate_figures.smk  

If you have multiple threads available you can up the number of cores (though the entire workflow shouldn't take longer then a couple of minutes on a regular pc).
The figures will be generated inside the 'output' folder, some dataframes will be generated inside the 'data' folder, prefixed with 'out_'.

## data

 - gadm41_BEL_4.json

geojson file for Belgium, downloaded from [gadm](https://gadm.org/download_country.html)

  - metadata.tsv

contains spatial information (on commune level) for all samples.

  - metagenomic_coverage.tsv

normalised coverage for all included viruses, as assayed in our [metagenomic study](https://doi.org/10.1101/2020.09.15.298042) 

  - qpcr_abs.tsv

Absolute quantification values (by qRT-PCR) for all seven included viruses.

  - sequences.fasta

Raw (nucleotide) sequences for all seven viruses.

  - trees/*

Different phylogenetic trees as inferred by iqtree (plotted in figure 1).
