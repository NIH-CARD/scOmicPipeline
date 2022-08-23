#! /usr/bin/env python

# cell annotation using the reference marker list and a rank_genes_groups at a specific resolution

#import h5py
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import json
import helpers
from helpers import MultiLineArgAndFileParser
from collections import namedtuple

def python_annotate(input, output, **kwargs):
    kwargs["input"] = input
    kwargs["out"] = output
    kwargs.update(kwargs['annotation']) #collapse 
    args = namedtuple("args", kwargs.keys())(*kwargs.values()) # turn kwargs into an object

    # load freom pickle if no anndata object provided
    adata = helpers.choose_adata_src(args.input)    
    # get the reference gene list from a path
    marker_ref = get_ref_list(args.marker_ref_path)
    # get the cell types after comparing 3 methods of overlapping using the function of cell_types()
    cell_types_all = compare_types(adata, marker_ref, rank_key=args.rank_key)
    # generate the final cell type dic
    final_type_dic = final_types(cell_types_all)

    if args.new_cluster_names:
        annotate_2(adata, args.new_cluster_names, key=args.key, out=args.out, project=args.project, figure_type=args.figure_type)
    else:
        annotate(adata, final_type_dic, key=args.key, out=args.out, project=args.project, figure_type=args.figure_type)

    
def get_ref_list(marker_ref_path):
    """
    Read the reference marker (a dictionary saved as pickle)
    """
    with open(marker_ref_path, "r") as f:
      marker_ref = json.load(f)

    return marker_ref

def cell_types(adata, marker_ref, rank_key=None, method=None, normalize=None):
    """
    Generate a list containing the cell type corresponding to the clusters at a specific resolution
    rank_key: the key of the rank_genes_groups at a resolution
    method: one of the 4 methods provided by sc.tl.marker_gene_overlap
    normalize: use when the method is "overlap_count"; options: "data", "reference", "None"
    Choose the type with the highest score to generate the list of cell types corresponding to each cluster
    """
    annotation = sc.tl.marker_gene_overlap(adata, marker_ref, key=rank_key, method=method, normalize=normalize)

    # sns.set(rc={'figure.figsize':(20, 20)})
    # sns.heatmap(annotation, annot=True)
    fig, ax = plt.subplots(figsize=(20,20))
    ax = sns.heatmap(annotation, annot=True)
    plt.savefig(rank_key+method+".pdf")

    cell_types = []
    for col in annotation.columns.to_list(): # columns are cluster names
        # the index column contains the cell type names
        # choosing the cell type (indice) with the highest value
        if annotation[col].max() > 0:
            cell_type = annotation[col].idxmax()
        else:
            cell_type = "Unknown"
        # generate the list of cell types, the sequence matches the cluster names
        cell_types.append(cell_type)

    return cell_types

def compare_types(adata, marker_ref, rank_key=None):
    """
    Try 3 methods ("overlap_count", jaccard", "overlap_coef") and merge them to a datafram
    The option of normalize in the method of "overlap_count" can change the type greatly
    """
    cell_types_1 = cell_types(adata, marker_ref, rank_key=rank_key, method="overlap_count", normalize="reference")
    cell_types_2 = cell_types(adata, marker_ref, rank_key=rank_key, method="jaccard")
    cell_types_3 = cell_types(adata, marker_ref, rank_key=rank_key, method="overlap_coef")
    cell_types_all = pd.DataFrame({"overlap_count": cell_types_1, "jaccard": cell_types_2, "overlap_coef": cell_types_3})

    return cell_types_all

def final_types(cell_types_all):
    """
    Pick the cell types that win the votes from the 3 methods.
    Convert it into a dictionary; also cast the numeric representation of clusters into string
    """
    final_type = cell_types_all.mode(axis=1)
    final_type_dic = final_type.iloc[:,0].to_dict()
    final_type_dic = {str(key): val for key, val in final_type_dic.items()}

    return final_type_dic

# add the annotation to adata
def annotate(adata, final_type_dic, key=None, out=None, project=None, figure_type=None):
    """
    Add the cell type annotation to the adata using a cluster-type dictionary, 
    instead of overwrite the original cluster label.
    """
    adata.obs[key+"_annotation"] = adata.obs[key].replace(final_type_dic)
    counts_per_type = adata.obs[key+"_annotation"].value_counts()
    per_type = 100*counts_per_type/counts_per_type.sum()
    sc.pl.umap(adata, color=key+"_annotation", save=project+"_"+key+"_ann."+figure_type)
    helpers.dump_to_pickle(out, adata)

    return per_type

def annotate_2(adata, new_cluster_names, show=True, key=None, out=None, project=None, figure_type=None):
    """
    This is scanpy default methods in annotating the clusters using self-define new_cluster_names.
    Take a list/array for new cell type names.
    However, the categories must be unique.
    The argument of new_cluster_names must be provided. It must be a list and have the same length as the number of clusters.
    """
    #make a copy of the origanl key
    adata.obs[key+"_orig"] = adata.obs[key]
    adata.rename_categories(key, new_cluster_names)
    sc.pl.umap(adata, color=key, legend_loc='on data', show=show, save=project+"."+figure_type)
    helpers.dump_to_pickle(out, adata)

def main():
    
    parser = MultiLineArgAndFileParser(fromfile_prefix_chars="@", description="Arguments for annotate cell types for clusters")
    # optional argument
    parser.add_argument("-i", "--input_file", type=str, help="path of the input of 'after_ranking_gene.pickle'", default="after_ranking_gene.pickle")
    parser.add_argument("-m", "--marker_ref_path", type=str, help="path of panglao reference markers", \
        default="../scanpy_scripts/reference_markers/marker_panglao_brain_dic_update.json")
    parser.add_argument("-o", "--out", type=str, help="path of the anndata object to be saved", default="after_annotated.pickle")
    parser.add_argument("-d", "--dpi", type=int, help="resolution of the output figure", default=80)
    parser.add_argument("-s", "--figsize", type=float, nargs=2, help="size of output figure, use 2 numbers, e.g., 2 2")
    parser.add_argument("-f", "--figure_type", type=str, help="define the export type of plot_type, e.g., png, pdf, or svg", default="pdf")
    parser.add_argument("-p", "--project", type=str, help="give the project name", default="")
    parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is show=True; provide no, false, or 0 to block print to screen")
    parser.add_argument("-k", "--key", type=str, help="Choose the key of leiden clustering", default="leiden_0.4")
    parser.add_argument("-r", "--rank_key", type=str, help="Choose the key of rank_genes_groups to be compared to marker_ref, \
        e.g., rank_genes_groups", default="rank_genes_groups_r0.6")
    parser.add_argument("-n", "--new_cluster_names", type=str, nargs="+", help="provide the cell type name corresponding to each cluster")

    # parser.set_defaults(func=annotate) #??
    args = parser.parse_args()
    # args.func(args)
    input_file = args.input_file
    marker_ref_path = args.marker_ref_path
    dpi = args.dpi
    out = args.out
    figsize = args.figsize
    figure_type = args.figure_type
    show = args.show
    project = args.project if (args.project == "") else ("_" + args.project)
    rank_key = args.rank_key
    key = args.key
    show = args.show
    new_cluster_names = args.new_cluster_names

    print("\nThe arguments are: ", args)

    # set scanpy parameters
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_version_and_date()
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=dpi, facecolor='white', figsize=figsize)

    # get the h5ad/pickle after ranking genes
    adata = helpers.choose_adata_src(input_file)
    # get the reference gene list from a path
    marker_ref = get_ref_list(marker_ref_path)
    # get the cell types after comparing 3 methods of overlapping using the function of cell_types()
    cell_types_all = compare_types(adata, marker_ref, rank_key=rank_key)
    # generate the final cell type dic
    final_type_dic = final_types(cell_types_all)

    if new_cluster_names:
        annotate_2(adata, new_cluster_names, show=show, key=key, out=out, project=project, figure_type=figure_type)
    else:
        annotate(adata, final_type_dic, key=key, out=out, project=project, figure_type=figure_type)


if __name__ == "__main__":
    main()
