#!/bin/bash

input_file=$1
module_path="../my_scanpy_modules"
ref_path="../scanpy_markers/reference_markers/marker_panglao_brain_dic_update.json"
echo "Processing ${input_file}"

qc="${module_path}/qc_1.py"	
cluster="${module_path}/cluster_2.py"
rank="${module_path}/marker_gene_3.py"
plot="${module_path}/plot_markers_4.py"
ann="${module_path}/annotate_5.py"

# qc
python $qc -i ${input_file} -n 9000 -m 15 -e yes
echo "finished QC"
echo "+++++++++++++++\n"

# cluster
python $cluster -i after_leiden.h5ad -r 1.0 -k leiden_1.0 -C leiden_1.0

echo "finished clustering"
echo "+++++++++++++++\n"

# marker genes
python $rank -i after_ranking_gene.h5ad -b leiden_1.0 -k rank_genes_groups_r1.0
echo "finished ranking"
echo "+++++++++++++++\n"

# plot markers for neurons: 
python $plot -S no -p neuron_marker -t umap -g MAP2
python $plot -S no -p glutamatergic -t umap -g SLC17A7 SLC17A6 GRIN1 GRIN2B GLS GLUL
echo "finished plotting\n"
echo "+++++++++++++++\n"

# annotate
python $ann -i after_annotated.h5ad -r rank_genes_groups_r1.0 -k leiden_1.0 -m ${ref_path}
echo "finished annotating"
echo "+++++++++++++++\n"
