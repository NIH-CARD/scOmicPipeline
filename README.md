[![PyPI version](https://badge.fury.io/py/card-scrnaseq-pipeline.svg)](https://badge.fury.io/py/card-scrnaseq-pipeline)

The scPipeline is a pipeline using scanpy to do single-cell analysis. We hope to make new additions to that soon.

## Installation
There are two installation methods:  
1. Install the full CLI tool:
    1. `pip install card-scrnaseq-pipeline`

2. Install the environment separately and run the tool manually:
    1. First, make sure you have miniconda or conda set up. If using Biowulf, first set up conda in your data directory, or use the miniconda modules on biowulf.
    2. Then, install create a new environment in this directory using
    ```
    conda env create -p scPipeline_env/ -f env_requirements.yml
    ```
    3. Activate the environment using `conda activate scPipeline_env/` and you're almost ready to go.


## Running 
There are several steps to the pipeline. We recommend creating a shell script to run each command consecutively.
Commands can be provided to the command line or in a script,
or they can be listed in a file, specified with @ in the command line. For example, 
```
scrnapipeline qc --input_file myfile.h5ad @qc_arguments.txt
```

The steps are:
- qc -- requires an input_file. There is no default.  
- cluster  
- ranking  
- plot_marker  
- annotate  


Run each of these in order. Provide the project argument for each command. The project name will be prefixed to any output files, including alternative output filenames specified by the user. This is helpful for keeping separate runs of the pipeline separate. Therefore,
we recommend the user always suggest a project name.
For convenience and consistency, each command already has default
output names for all files; if the user changes these names, the input_file name for the following command must be specified.
Apart from the project argument
and the input_file for the qc command, nearly all arguments are optional.

### Running each step
To learn about the possible arguments of each subcaommand and what each means, call each command with the --help flag.
For example:
```
# Using the pip package command
plot_markers --help

# Running the command manually if the pip package was not installed
python3 src/scrnapipeline/scrnapipeline.py plot_marker --help
```  
The pipeline can be run either command by command on the command line, or in a script.
Running using a script has the disadvantage that any error will require the script to be run again (in the current version. In the future, we will all restarting from the most recent complete step). The advantage is a record is created about what was done.
However, that can also be achieved by providing arguments in a file for each subcommand. 
To do that, specify the name of the file with @ at the command line. E.g.:
```
scrnapipeline qc --input_file myfile.h5ad @qc_arguments.txt
``` 
Within argument files, each argument consists of a flag, an equals sign, and a setting. If listing something, 
simply add a space between each item, e.g. `--gene_list=MAP2 FOX4 APOE`

To create a script, simply put all commands and arguments in a file and name it my_scrnaseq_pipeline_script.sh.
At the top, put the following two lines:
```
#!/bin/bash
input_file=$1 # $1 means the first argument if this script is run from the command line
              # You can also place the actual file name here if you prefer.
```             
Run using 
```
bash my_scrnaseq_pipeline_script.sh data/my_data_file.h5ad
```

The following is an example of a full execution script. The same format can be followed to run each command individually:
```
#!/bin/bash

input_file=$1 # $1 means the first argument if this script is run from the command line
              # You can also place the actual file name here if you prefer.
ref_path="data/reference_markers/marker_panglao_brain_dic_update.json"

echo "Processing ${input_file}"

# qc
scrnapipeline qc -i ${input_file} \
                          --project testing \
                          --n_genes_by_counts 9000 \
                          --pct_counts_mt 15 \
                          --exclude_highly_expressed yes \
                          --show no
echo "finished QC"
echo "+++++++++++++++\n"

# cluster
scrnapipeline cluster --project testing \
                               -r 1.0 \
                               -k leiden_1.0 \
                               -C leiden_1.0 \
                               --show no

echo "finished clustering"
echo "+++++++++++++++\n"

# marker genes
scrnapipeline ranking --project testing \
                               -k rank_genes_groups_r1.0 \
                               --groupby leiden_1.0 \
                               --show no
echo "finished ranking"
echo "+++++++++++++++\n"

# plot markers for neurons: 
#python3 $scrnapipeline plot_marker --project testing \
#                                   -S no \
#                                   -p neuron_marker \
#                                   -t umap \
#                                   -g MAP2

scrnapipeline plot_marker --project testing \
                                   -S no \
                                   -t umap \
                                   -g MAP2 
echo "finished plotting\n"
echo "+++++++++++++++\n"

# annotate
scrnapipeline annotate --project testing \
                                -r rank_genes_groups_r1.0 \
                                -k leiden_1.0 \
                                -S no \
                                -m ${ref_path}
echo "finished annotating"
echo "+++++++++++++++\n"
```

### Other info
Reference markers must be provided in json format.

A utility is provided to convert pickle formatted files to h5ad.  
`scrnapipeline convert_pickle my_pickle_file.pickle`

## Cleanup
After running, there may be several pickle files. The only one that is important is the (your_project_name)_after_annotate.pickle file.
When using the pip installable, the default is to automatically delete all pickle files except the after_annotate.pickle file as each becomes unnecessary following the successful completion of each step. When using the python hooks, however, that is not true and must be done manually. This is done for backwards compatibility.
Feel free to delete the rest as they may be large. Depending on space, it may be necessary to run the commands one at a time in order to delete
the previous pickle file.




