#########
# This Snakefile is adapted from bash modules made by Lirong. These
# are available for reference from 
# CARD_iNDI/Data/Discovery_planning phase/Informatics_Lirong/projects_lirong/my_scanpy_modules
# Will also need to run conda install -c conda-forge leidenalg
# #######

from lib import cluster_specific
import h5py
import scanpy as sc
import yaml
import pickle
import anndata

configfile:  "config/config_snakemake_bash.yaml"

shell.prefix(
    'set -euo pipefail; export R_PROFILE_USER=; export TMPDIR={};'
    .format(cluster_specific.tempdir_for_biowulf())
)
shell.executable('/bin/bash')

sc.settings.verbosity = 3
sc.logging.print_versions()
sc.settings.set_figure_params(config, facecolor="white") # frameon=False
sc.logging.print_version_and_date()

# We only do this at the end because there is a bug in the code that 
# does not encode state properly in the h5ad format, messing up 
# the rank_genes_groups method in the marker_gene rule
def get_final_file_format(name):
    if config['convert_to_h5ad']:
        return name + ".h5ad"
    else:
        return name + ".pickle"

rule all:
    input:
        get_final_file_format(config["project"] + "_after_annotate"),
        config["project"] + "_after_plotting.txt",
        config["project"] + "_after_leiden.pickle",
        config['project'] + "_after_ranking_gene.pickle",



rule qc_1:
    input: "data/ngn_d0_28_raw.h5ad"
    output: temp(config["project"] + "_count_after_QC.pickle")
    log: "logs/qc_log.txt"
    params: 
        file="config/qc_arguments.txt"
    shell:
        """
        python3 bash_pipeline/qc_1.py -o {output} \
            --project {config[project]} \
            @{params.file} 2>&1 | \
        tee {log}
        """
        
rule cluster_2:
    input: config["project"] + "_count_after_QC.pickle"
    output: temp(config["project"] + "_after_leiden.pickle")
    log: "logs/cluster_log.txt"
    params:
        file="config/cluster_arguments.txt"
    shell:
        """
        python3 bash_pipeline/cluster_2.py -i {input} \
            --project {config[project]} \
            -o {output} @{params.file} 2>&1 | \
        tee {log}
        """

rule marker_gene_3:
    input: config['project'] + "_after_leiden.pickle"
    output: temp(config['project'] + "_after_ranking_gene.pickle")
    log: "logs/ranking_log.txt"
    params:
        file="config/marker_arguments.txt"
    shell:
        """
        python3 bash_pipeline/marker_gene_3.py -i {input} \
                --project {config[project]} \
                -o {output} @{params.file} 2>&1 | 
        tee {log}
        """

rule plot_markers_4:
    input: config['project'] + "_after_ranking_gene.pickle"
    output: temp(touch(config["project"] + "_after_plotting.txt"))
    log: "logs/plot_markers_log.txt"
    params:
        file="config/plotting_arguments.txt"
    shell:
        """
        python3 bash_pipeline/plot_markers_4.py -i {input} \
                --project {config[project]} \
                @{params.file} 2>&1 | \
        tee {log}
        """

rule annotate_5:
    input: file=config['project'] + "_after_ranking_gene.pickle"
    output: config['project'] + "_after_annotate.pickle"
    log: "logs/annotate_log.txt"
    params:
        file="config/annotation_arguments.txt"
    shell:
        """
        python3 bash_pipeline/annotate_5.py -i {input.file} \
                --project {config[project]} \
                -o {output} @{params.file} 2>&1 | \
        tee {log}
        """

if config['convert_to_h5ad']:
    rule convert_to_h5ad:
        input: config['project'] + "_after_annotate.pickle"
        output: config['project'] + "_after_annotate.h5ad"
        run:
            with open(input, "rb") as f:
                adata = pickle.load(f)
                adata.write(output, compression='gzip')

            


