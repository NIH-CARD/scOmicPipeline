The scPipeline is a pipeline using scanpy to do single-cell analysis. We hope to make new additions to that soon.

## Installation

1. First, make sure you have miniconda or conda set up. If using Biowulf, first set up conda in your data directory, or use the miniconda modules on biowulf.

2. Then, install create a new environment in this directory using
```
conda env create -p scPipeline_env/ -f env_requirements.yml
```
3. Activate the environment using `conda activate scPipeline_env/` and you're almost ready to go.

## Setup

1. Set the arguments of config/config_snakemake_bash.yaml
These are the name of the projet, which will be prepended to all output files, as well as
whether to convert the final output database to h5ad format when finished.

2. Each of the arguments.txt files in the config folder pertains to one of the pipeline steps.
Go through each and set the arguments as required. 
Within these files, each argument consists of a flag, an equals sign, and a setting. If listing something, 
simply add a space between each item, e.g. `--gene_list=MAP2 FOX4 APOE`
To learn about the possible arguments and what each means, call each command with the --help flag.
For example:
```
python3 bash_pipeline/plot_markers_4.py --help
```

## Running 
Run the snakemake command, followed by how many cores to use. If using a Biowulf node,
consider how many were requested for that node. Do not use more than the maximum number 
of cores available.
```
snakmake -c4
```






