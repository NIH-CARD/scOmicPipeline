The scPipeline is a pipeline using scanpy to do single-cell analysis. We hope to make new additions to that soon.

## Installation

1. First, make sure you have miniconda or conda set up. If using Biowulf, first set up conda in your data directory, or use the miniconda modules on biowulf.

2. Then, install create a new environment in this directory using
```
conda env create -p scPipeline_env/ -f env_requirements.yml
```
3. Activate the environment using `conda activate scPipeline_env/` and you're almost ready to go.


## Running 
There are several steps to the pipeline. Ideally, create a shell script that runs each individually.  

The steps are:  
qc -- requires an input_file. There is no default.  
cluster  
ranking  
plot_marker  
annotate  

Run each of these in order. Provide the project argument for each command. The project name will be prefixed to any output files, including alternative output filenames specified by the user. This is helpful for keeping separate runs of the pipeline separate. Therefore,
we recommend the user always suggest a project name.
For convenience and consistency, each command already has default
output names for all files; if the user changes these names, the input_file name for the following command must be specified.
Apart from the project argument
and the input_file for the qc command, nearly all arguments are optional.

### Providing arguments
To learn about the possible arguments and what each means, call each command with the --help flag.
For example:
```
python3 bash_pipeline/plot_markers_4.py --help
```  

Commands can be provided to the command line or in a script,
or they can be listed in a file, specified with @ in the command line. For example, 
```
bash_pipeline/scRNAPipeline qc --input_file myfile.h5ad @qc_arguments.txt
```  
Within argument files, each argument consists of a flag, an equals sign, and a setting. If listing something, 
simply add a space between each item, e.g. `--gene_list=MAP2 FOX4 APOE`

<<<<<<< HEAD
## Running 
There are several steps to the pipeline. Ideally, create a shell script that runs each individually.
Commands can be provided to the command line or in a script,
or they can be listed in a file, specified with @ in the command line. For example, 
bash_pipeline/scRNAPipeline qc --input_file myfile.h5ad @qc_arguments.txt.

The steps are:
qc -- requires an input_file. There is no default.
cluster
ranking
plot_marker
annotate

Run each of these in order. Provide the project argument for each command. The project name will be prefixed to any output files, including alternative output filenames specified by the user. This is helpful for keeping separate runs of the pipeline separate. Therefore,
we recommend the user always suggest a project name.
For convenience and consistency, each command already has default
output names for all files; if the user changes these names, the input_file name for the following command must be specified.
Apart from the project argument
and the input_file for the qc command, nearly all arguments are optional.

An example can be found in the example_run.sh file.
The example can be run by simply specifying an input file.
You may use the script as template to achieve similar results.
## Running 
Run the snakemake command, followed by how many cores to use. If using a Biowulf node,
consider how many were requested for that node. Do not use more than the maximum number 
of cores available.
### Example script
An example can be found in the example_run.sh file.
The example can be run by simply specifying an input file.
You may use the script as a template to achieve similar results.
To run the script as-is, use the following command using your data file.
```
bash example_run.sh data/my_data_file.h5ad
```

### Other info
Reference markers must be provided in json format.

A utility is provided to convert pickle formatted files to h5ad.
Simply run bash_pipeline/scRNAPipeline convert_pickle my_pickle_file.pickle

## Cleanup
After running, there will be several pickle files. The only one that is important is the (your_project_name)_after_annotate.pickle file.
Feel free to delete the rest as they may be large. Depending on space, it may be necessary to run the commands one at a time in order to delete
the previous pickle file.




