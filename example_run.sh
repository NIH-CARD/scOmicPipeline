#!/bin/bash

input_file=$1 # $1 means the first argument if this script is run from the command line
              # You can also place the actual file name here if you prefer.
module_path="bash_pipeline"
ref_path="data/reference_markers/marker_panglao_brain_dic_update.json"
echo "Processing ${input_file}"

scRNAPipeline="scrnapipeline"

# qc
$scRNAPipeline qc -i ${input_file} \
                          --project testing \
                          --n_genes_by_counts 9000 \
                          --pct_counts_mt 15 \
                          --exclude_highly_expressed yes \
                          --show no
echo "finished QC"
echo "+++++++++++++++\n"

# cluster
$scRNAPipeline cluster --project testing \
                               -r 1.0 \
                               -k leiden_1.0 \
                               -C leiden_1.0 \
                               --show no

echo "finished clustering"
echo "+++++++++++++++\n"

# marker genes
$scRNAPipeline ranking --project testing \
                               -k rank_genes_groups_r1.0 \
                               --groupby leiden_1.0 \
                               --show no
echo "finished ranking"
echo "+++++++++++++++\n"

# plot markers for neurons: 
#python3 $scRNAPipeline plot_marker --project testing \
#                                   -S no \
#                                   -p neuron_marker \
#                                   -t umap \
#                                   -g MAP2

$scRNAPipeline plot_marker --project testing \
                                   -S no \
                                   -t umap \
                                   -g MAP2 
echo "finished plotting\n"
echo "+++++++++++++++\n"

# annotate
$scRNAPipeline annotate --project testing \
                                -r rank_genes_groups_r1.0 \
                                -k leiden_1.0 \
                                -S no \
                                -m ${ref_path}
echo "finished annotating"
echo "+++++++++++++++\n"
