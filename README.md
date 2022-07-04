The scPipeline is a pipeline using scanpy to do single-cell analysis. We hope to make new additions to that soon.

## Installation

1. First, make sure you have miniconda or conda set up. IF using Biowulf, first set up conda in your data directory, or use the miniconda modules on biowulf.

2. Then, install create a new environment in this directory using
```
conda env create -p scPipeline_env/ -f env_requirements.yml
```
3. Activate the environment using `conda activate -p scPipeline_env/` and you're almost ready to go.

## Setup

Go through each arguments.txt file in the config folder and set the arguments as required. 
Within these files, each argument consists of a flag, an equals sign, and a setting. If listing something, 
simply add a space between each item, e.g. `--gene_list=MAP2 FOX4 APOE`
To learn about the possible arguments and what each means, call each command with the --help flag.
For example:
```
python3 bash_pipeline/plot_markers_4.py --help
```

## Running 







