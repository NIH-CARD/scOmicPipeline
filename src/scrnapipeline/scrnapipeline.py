#!/usr/bin/env python3
import argparse
import json
import os
import pickle
import h5py
from collections import namedtuple

#import sys
#sys.path.append(os.getcwd)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from . import helpers 
from .helpers import MultiLineArgAndFileParser


#### QC
def setup(dpi=None, figsize=None, **kwargs):
    """
    Set up the displaying and output parameters
    """
    sc.settings.verbosity = 3 
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=dpi, facecolor="white", figsize=figsize) # frameon=False
    sc.logging.print_version_and_date()


def qc(input_file=None, n_top=None, show=None, project=None, 
figure_type=None,min_genes=None, min_cells=None, n_genes_by_counts=None, pct_counts_mt=None, **kwargs):
    """
    This is quality check for scRNA-seq data and do some filterings.
    """
    if input_file.endswith(".h5"):
        adata = sc.read_10x_h5(input_file)
    elif input_file.endswith(".h5ad"):
        adata = sc.read_h5ad(input_file)
    else:
        # There is a potential bug in saving file later if use cache=True
        # here input is the directory name which contains matrix file
        adata = sc.read_10x_mtx(input_file, var_names='gene_symbols', cache=True)

    adata.var_names_make_unique()

    #check the most abundantly expressed genes
    sc.pl.highest_expr_genes(adata, n_top=n_top, show=show, save=project+"."+figure_type)

    ## Basic filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    ## Calculate the percentage of mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    ## Check sequencing quality
    #choose the threthold of gene numbers to remove, e.g., n_genes = 4500
    #choose the threthold of mitochondial genes to remove, e.g., percent_mito = 0.15
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                jitter=0.4, multi_panel=True, show=show, save=project+"_qc."+figure_type)

    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=show, save="1"+project+"."+figure_type)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=show, save="2"+project+"."+figure_type)

    adata = adata[adata.obs.n_genes_by_counts < n_genes_by_counts, :]
    adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]

    return adata

def normalize_regress(adata, show=None, project=None, figure_type=None, exclude_highly_expressed=None, **kwargs):
    """
    Normalize and logarithmize adata, retrieve highly_variable_genes, and futher scale and remove outliers.
    """
    # Scale and logarithmize the data
    # here can exclude highly expresseddd gene using
    # sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=exclude_highly_expressed)
    sc.pp.log1p(adata)
    # optional: store the unnormalized data in .raw
    adata.raw=adata
    # Choosing highly-variable genes for further analysis
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata, show=show, save=project+"."+figure_type)

    # this below filtering step is unnecessary if only just for the purpos of PCA
    adata = adata[:, adata.var.highly_variable]

    # Regression out effects of total counts per cell and the percentage of mitochondrial genes expressed
    # and the data to unit variance, Clip values exceeding standard deviation 10.
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    return adata

def pca(adata, color_gene=None, show=None, project=None, figure_type=None, out=None, **kwargs):

    # Principal component analysis
    sc.tl.pca(adata, svd_solver='arpack')
    # remember to change the parameter of color = "-"
    sc.pl.pca(adata, color=color_gene, show=show, save=project+"."+figure_type)
    sc.pl.pca_variance_ratio(adata, log=True,show=show, save=project+"."+figure_type)

    # adata.write_h5ad(out) # recommend h5ad type
    # Dump adata out
    helpers.dump_to_pickle(out, adata)

    return adata

def qc_main(args):

    input_file= args.input_file
    out = args.out
    dpi = args.dpi
    figsize = args.figsize
    figure_type = args.figure_type
    show = args.show
    project = args.project if (args.project == "") else ("_" + args.project)

    min_genes = args.min_genes
    min_cells = args.min_cells
    n_genes_by_counts = args.n_genes_by_counts
    pct_counts_mt = args.pct_counts_mt
    n_top = args.n_top
    color_gene = args.color_gene
    exclude_highly_expressed = args.exclude_highly_expressed

    print("The arguments are: ", args)
    print()
    print("exclude_highly_expressed is ", exclude_highly_expressed)
    print()

    setup(dpi=dpi, figsize=figsize)
    adata = qc(input_file=input_file, n_top=n_top, show=show, project=project, figure_type=figure_type, 
    min_genes=min_genes, min_cells=min_cells, n_genes_by_counts=n_genes_by_counts, pct_counts_mt=pct_counts_mt)
    adata = normalize_regress(adata, show=show, project=project, figure_type=figure_type, exclude_highly_expressed=exclude_highly_expressed)
    pca(adata, color_gene=color_gene, show=show, project=project, figure_type=figure_type, out=out)

### Cluster 
"""
Wrapper method for cluster
It also deletes the previous, no longer necessary, pickle of the object.
"""
def python_cluster(input, output, **kwargs):
    kwargs["input"] = input
    kwargs["out"] = output
    args = namedtuple("args", kwargs.keys())(*kwargs.values()) # turn kwargs into an object
    cluster_main(args)

def cluster_main(args):
    """
    Clustering cells after computing pca and neiborhood distances.
    """
    
    print()
    print(f"The arguments are {args}")
    
    input = args.input
    out = args.out
    dpi = args.dpi
    figsize = args.figsize
    figure_type = args.figure_type
    show = args.show
    
    project = args.project if (args.project == "") else ("_" + args.project)

    resolution = args.resolution
    n_neighbors = args.n_neighbors
    n_pcs = args.n_pcs
    #method = args.method
    #metric = args.metric
    color_gene = args.color_gene
    key_added = args.key_added

    # set scanpy parameters
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    # in scanpy version 1.6.1 tutorial: sc.logging.print_header()
    sc.logging.print_version_and_date()
    sc.logging.print_versions()
    # default figsize=None, means it doesn't change the seaborn defined default parameters
    sc.settings.set_figure_params(dpi=dpi, facecolor='white', figsize=figsize)
    
    # load freom pickle if no anndata object provided
    adata = helpers.choose_adata_src(input)    

    ### Computing, embedding, and clustering the neighborhood graph
    # defaults are: n_neighbors= 15, n_pcs=None
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    # plot umap using raw data: normalized and logarithimized but not regressed out
    # sc.pl.umap(adata, color=color, save="_on_raw_"+project+"."+figure_type)
    # plot umap using scaled and corrected gene expression
    # sc.pl.umap(adata, color=color, use_raw=False, save="_"+project+"."+figure_type)

    # cluster using leeiden graph-clustering method
    # default resolution=1.0
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    sc.pl.umap(adata, color=color_gene, show=show, save="_after_leiden"+project+"."+figure_type)

    # Dump adata out
    helpers.dump_to_pickle(out, adata)

    return adata

### Ranking
"""
Python wrapper method for the ranking method
"""
def python_ranking(input, output, **kwargs):
    kwargs["input"] = input
    kwargs["out"] = output
    kwargs.update(kwargs['marker_gene']) #collapse 
    args = namedtuple("args", kwargs.keys())(*kwargs.values()) # turn kwargs into an object
    ranking_main(args)

def ranking_main(args):
    """
    This is finding marker genes for each cell cluster.
    """
    
    print()
    print("The arguments are: ", args)

    input = args.input
    out = args.out
    dpi = args.dpi
    figsize = args.figsize
    figure_type = args.figure_type
    show = args.show
    project = args.project if (args.project == "") else ("_" + args.project)

    groupby = args.groupby
    n_genes = args.n_genes
    export_num = args.export_num
    method = args.method
    corr_method = args.corr_method
    reference = args.reference
    groups = args.groups
    key_added = args.key_added

    # set scanpy parameters
    # reduce the verbosity from 3 to 2 in the setting of logging output
    sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_versions()
    sc.logging.print_version_and_date()
    sc.settings.set_figure_params(dpi=dpi, facecolor='white', figsize=figsize) #frame=False

    adata = helpers.choose_adata_src(input)    

    # one vs rest comparison using Mann-Whitney-U-test (recommend)
    # default: groupby="leiden", groups="all", method="wilcoxon", key_added=None
    # sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    # sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
    sc.tl.rank_genes_groups(adata, groupby=groupby, groups=groups, reference=reference, method=method, corr_method=corr_method, key_added=key_added)
    sc.pl.rank_genes_groups(adata, key=key_added, n_genes=n_genes, show=show, sharey=False, save=project+"."+figure_type)


    # export a dataframe for marker genes
    if key_added is None:
        result = adata.uns['rank_genes_groups']
    else:
        result = adata.uns[key_added]

    groups = result['names'].dtype.names

    top_marker_genes = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(export_num)

    if key_added is None:
        top_marker_genes.to_csv("rank_genes_groups" + project + '.csv')
    else: 
        top_marker_genes.to_csv(key_added + project + '.csv')

    # Dump adata out
    helpers.dump_to_pickle(out, adata)

    return adata



######## Plot Makers
"""
Wrapper method for plot_makers
"""
def python_plot_makers(input, **kwargs):
    kwargs["input"] = input
    kwargs.update(kwargs['plot_makers']) #collapse 
    args = namedtuple("args", kwargs.keys())(*kwargs.values()) # turn kwargs into an object
    plot_makers_main(args)


def plot_makers_main(args):
    """
    Visualize the differential expression of marker genes across clusters.
    """

    print("\nThe arguments are: ", args)

    
    input_file = args.input_file
    figure_type = args.figure_type
    dpi = args.dpi
    figsize = args.figsize
    show = args.show
    project = args.project if (args.project == "") else ("_" + args.project)
    use_raw = args.use_raw

    groupby = args.groupby
    n_genes = args.n_genes
    groups = args.groups
    plot_type = args.plot_type
    gene_list = args.gene_list
    series = args.series
    adata = helpers.choose_adata_src(input_file) 

    setup(dpi, figsize)

    # choose the plotting types
    if len(plot_type) != 0:
        if  "violin" in plot_type:
            print("plotting violin")
            sc.pl.violin(adata, gene_list, groupby=groupby, show=show, use_raw=use_raw, save=helpers.get_save_name(series=series, project=project, figure_type=figure_type))
        if "dotplot" in plot_type:
            print("plotting dotplot")
            sc.pl.dotplot(adata, gene_list, groupby=groupby, show=show, use_raw=use_raw, save=helpers.get_save_name(series=series, project=project, figure_type=figure_type))    
        if "stacked_violin" in plot_type:
            print("plotting stacked_violin")
            sc.pl.stacked_violin(adata, gene_list, groupby=groupby, rotation=90, show=show, use_raw=use_raw, save=helpers.get_save_name(series=series, project=project, figure_type=figure_type))
        if "rank_genes_groups_violin" in plot_type:
            print("plotting rank_genes_groups_violin")
            sc.pl.rank_genes_groups_violin(adata, groups=groups, n_genes=n_genes, show=show, use_raw=use_raw, save=helpers.get_save_name(series=series, project=project, figure_type=figure_type))
        if "umap" in plot_type:
            print("plotting umap")
            sc.pl.umap(adata, color=gene_list, show=show, use_raw=use_raw, save=helpers.get_save_name(prefix="_gene_expr", series=series, project=project, figure_type=figure_type))
    else:
        print("No such type of plot")

## Annotate
# cell annotation using the reference marker list and a rank_genes_groups at a specific resolution


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
def annotate(adata, final_type_dic, key=None, out=None, project=None, figure_type=None, show=True):
    """
    Add the cell type annotation to the adata using a cluster-type dictionary, 
    instead of overwrite the original cluster label.
    """
    adata.obs[key+"_annotation"] = adata.obs[key].replace(final_type_dic)
    counts_per_type = adata.obs[key+"_annotation"].value_counts()
    per_type = 100*counts_per_type/counts_per_type.sum()
    sc.pl.umap(adata, show=show,  color=key+"_annotation", save=project+"_"+key+"_ann."+figure_type)
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

def annotate_main(args):
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
        annotate(adata, final_type_dic, key=key, out=out, project=project, figure_type=figure_type, show=show)


def convert_pickle_main(args):
    input_file = args.input_file
    output_name = os.path.split(input_file)[0] + os.path.splitext(os.path.split(input_file)[1])[0] + ".h5ad"

   
    with open(input_file, "rb") as f:
        adata = pickle.load(f)
        adata.write(output_name, compression='gzip')

def main():
    # top-level parser
    parser = MultiLineArgAndFileParser()
    subparsers = parser.add_subparsers()

    ### QC
    qc_parser = subparsers.add_parser("qc", fromfile_prefix_chars="@", description="Arguments for scRNA-seq QC")
    
    # basic parameters
    qc_parser.add_argument("-p", "--project", type=str, help="a project name", default="")
    qc_parser.add_argument("-i", "--input_file", type=str, help="the path of the input .h5 file or the matrix folder", required=True)
    qc_parser.add_argument("-o", "--out", type=str, help="the path and file name of the output anndata", default="count_after_QC.pickle")
    qc_parser.add_argument("-d", "--dpi", type=int, help="the resolution of the output figure", default=80)
    qc_parser.add_argument("-s", "--figsize", type=float, nargs=2, help="the size of the output figure, use 2 numbers, e.g., 2 2", default=None)
    qc_parser.add_argument("-f", "--figure_type", type=str, help="the type of the output figure, e.g., pdf, png, or svg", default="pdf")
    qc_parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="block output figures on the screen by providing no, false, or 0")
    
    # filtering and plotting parameters
    qc_parser.add_argument("-g", "--min_genes", type=int, help="cell filter: the minimal number of genes which a cell should have", default=200)
    qc_parser.add_argument("-c", "--min_cells", type=int, help="gene filter: the minimal number of cells which a gene should be in", default=3)
    qc_parser.add_argument("-n", "--n_genes_by_counts", type=int, help="the threthold of gene counts", default=8000)
    qc_parser.add_argument("-m", "--pct_counts_mt", type=float, help="the threthold of mitochondrial genes percentage", default=5.0)
    qc_parser.add_argument("-t", "--n_top", type=int, help="the number of genes to plot in the highest expressed gene plot", default=20)
    qc_parser.add_argument("-C", "--color_gene", type=str, help="the gene whoes expression level is to be colored in the pca plot", default="MAP2")
    qc_parser.add_argument("-e", "--exclude_highly_expressed", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="exclude highly expressed genes when do normalization by providing yes, no, or 1 ")
    
    qc_parser.set_defaults(func=qc_main)

    ### CLUSTER
    cluster_parser = subparsers.add_parser("cluster", fromfile_prefix_chars="@", description="Arguments for scRNA-seq Clustering")
    
    # basic parameters
    cluster_parser.add_argument("-i", "--input", type=str, help="the path and name of count_after_QC.pickle file", default="count_after_QC.pickle")
    cluster_parser.add_argument("-d", "--dpi", type=int, help="the resolution of the output figure", default=80)
    cluster_parser.add_argument("-f", "--figure_type", type=str, help="the export type of plots, e.g., png, pdf, or svg", default="pdf")
    cluster_parser.add_argument("-p", "--project", type=str, help="the project name", default="")
    cluster_parser.add_argument("-o", "--out", type=str, help="the path and file name to save the anndata object", default="after_leiden.pickle")
    cluster_parser.add_argument("-s", "--figsize", type=float, nargs=2, help="the size of output figure, use 2 numbers, e.g., 2 2")
    cluster_parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="block output figures on the screen by providing no, false, or 0")
    
    # umap parmeters
    cluster_parser.add_argument("-n", "--n_neighbors", type=int, help="the size of local neiborhood for manifold approximation", default=15)
    cluster_parser.add_argument("-P", "--n_pcs", type=int, help="the number of PCs to use", default=None)
    cluster_parser.add_argument("-m", "--method", type=str, help="the method for neighborhood graph, either ‘umap’, ‘gauss’, ‘rapids’", default="umap")
    cluster_parser.add_argument("-M", "--metric", type=str, help="the metric for neighborhood graph, [‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’], Literal[‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’],", default="euclidean")

    # leiden parameters
    cluster_parser.add_argument("-r", "--resolution", type=float, help="the resolution for leiden", default=1.0)

    # color parameters and key names to be stored in adata
    cluster_parser.add_argument("-C", "--color_gene", type=str, nargs="*", help="define a list of genes (e.g., MAP2 TEME199 TMEM106B), a key of leiden (e.g., 'leiden' or other key_added like 'leiden_0.6'), or both as color hues in umap plot", default="leiden")
    # parser.add_argument("-g", "--gene_list", type=str, nargs="+", action="store", dest="list", help="define a list of genes to show in umap, e.g., MAP2 TEME199 NIL", default=['leiden'])
    cluster_parser.add_argument("-k", "--key_added", type=str, help="the key name of a ledien anaysis to be addeed to anndata", default='leiden')

    cluster_parser.set_defaults(func=cluster_main)
    
    ### RANKING
    
    marker_parser = subparsers.add_parser("ranking", fromfile_prefix_chars="@", description="Arguments for scRNA-seq Ranking")
    # basic parameters
    marker_parser.add_argument("-i", "--input", type=str, help="the path of after_leiden.pickle", default="after_leiden.pickle")
    marker_parser.add_argument("-d", "--dpi", type=int, help="the resolution of the output figure", default=80)
    marker_parser.add_argument("-f", "--figure_type", type=str, help="the type of plots, e.g., png, pdf, or svg", default="pdf")
    marker_parser.add_argument("-p", "--project", type=str, help="the project name", default="")
    marker_parser.add_argument("-o", "--out", type=str, help="the path and name of the output anndata", default="after_ranking_gene.pickle")
    marker_parser.add_argument("-s", "--figsize", type=float, nargs=2, help="the size of the output figure, use 2 numbers, e.g., 2 2")
    marker_parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="block outputing figures on the screen by providing no, false, or 0")

    # rank gene parmeters
    marker_parser.add_argument("-b", "--groupby", type=str, help="the key of the obs grouping to be condisder, e.g., leiden, leiden_0.6", default="leiden")
    marker_parser.add_argument("-n", "--n_genes", type=int, help="the number of top marker genes to plot", default=25)
    marker_parser.add_argument("-N", "--export_num", type=int, help="the number of top marker genes to be exported to csv", default=50)
    marker_parser.add_argument("-m", "--method", type=str, help="the method for differential expression analysis, either t-test, wilcoxon, t-test_overestim_var, logreg", default="wilcoxon")
    marker_parser.add_argument("-c", "--corr_method", type=str, help="the p-value correction method, either benjamini-hochberg or bonferroni", default="benjamini-hochberg")
    marker_parser.add_argument("-R", "--reference", type=str, help="the cell group to be compared with", default="rest")
    marker_parser.add_argument("-g", "--groups", type=str, nargs="+", help="the subset cell groups to compare, using the group names in anndata.obs, e.g. ['g1', 'g2', 'g3'] or ['0', '1', '2']", default="all")
    marker_parser.add_argument("-k", "--key_added", type=str, help="the key in adata.uns indicating where the information to be saved to")

    marker_parser.set_defaults(func=ranking_main)
    

    ### PLOT_MARKER 
    plot_marker_parser = subparsers.add_parser("plot_marker", fromfile_prefix_chars="@", description="Arguments for scRNA-seq plotting")
    
    # basic parameters
    plot_marker_parser.add_argument("-i", "--input_file", type=str, help="path of the input after_leiden.pickle file", default="after_ranking_gene.pickle")
    plot_marker_parser.add_argument("-d", "--dpi", type=int, help="the resolution of the output figure", default=80)
    plot_marker_parser.add_argument("-s", "--figsize", type=float, nargs=2, help="the size of the output figure, use 2 numbers, e.g., 2 2")
    plot_marker_parser.add_argument("-f", "--figure_type", type=str, help="define the export type of plot_type, e.g., png, pdf, or svg", default="pdf")
    plot_marker_parser.add_argument("-p", "--project", type=str, help="give the project name", default="")
    plot_marker_parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is show=True; provide no, false, or 0 to block print to screen")
    plot_marker_parser.add_argument("-r", "--use_raw", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is using raw; provide no, false, or 0 to not use raw data")

    # plotting gene parmeters
    plot_marker_parser.add_argument("-g", "--gene_list", type=str, nargs="*", help="define a list of genes (e.g., MAP2 TEME199 TMEM106B) for plotting", default=["MAP2"])
    plot_marker_parser.add_argument("-t", "--plot_type", type=str, nargs="*", help="define the plotting types, e.g., dotplot, violin, stacked_violin, and rank_genes_groups_violin", default=["violin", "dotplot", "stacked_violin", "rank_genes_groups_violin", "umap"])
    plot_marker_parser.add_argument("-b", "--groupby", type=str, help="the key of the obs grouping to be considered, e.g., leiden, leiden_0.6")
    #plot_marker_parser.add_argument("-b", "--groupby", type=str, help="the key of the obs grouping to be condisder, e.g., leiden, leiden_0.6", default="leiden")
    plot_marker_parser.add_argument("-n", "--n_genes", type=int, help="number of genes used for plotting rank_genes_groups_violin", default=8)
    plot_marker_parser.add_argument("-G", "--groups", type=str, nargs="*", help="choose the subset cell groups for plotting rank_genes_groups_violin, e.g. ['g1', 'g2', 'g3'] or ['0', '1', '2']")
    plot_marker_parser.add_argument("-X", "--series", type=str, help="A name to append after the project name for output files.")

    plot_marker_parser.set_defaults(func=plot_makers_main)


    ### ANNOTATE 
    anno_parser = subparsers.add_parser('annotate', fromfile_prefix_chars="@", description="Arguments for annotate cell types for clusters")
    
    # optional argument
    anno_parser.add_argument("-i", "--input_file", type=str, help="path of the input of 'after_ranking_gene.pickle'", default="after_ranking_gene.pickle")
    anno_parser.add_argument("-m", "--marker_ref_path", type=str, help="path of panglao reference markers", \
        default="../scanpy_scripts/reference_markers/marker_panglao_brain_dic_update.json")
    anno_parser.add_argument("-o", "--out", type=str, help="the path of the anndata object to be saved", default="after_annotated.pickle")
    anno_parser.add_argument("-d", "--dpi", type=int, help="resolution of the output figure", default=80)
    anno_parser.add_argument("-s", "--figsize", type=float, nargs=2, help="size of output figure, use 2 numbers, e.g., 2 2")
    anno_parser.add_argument("-f", "--figure_type", type=str, help="define the export type of plot_type, e.g., png, pdf, or svg", default="pdf")
    anno_parser.add_argument("-p", "--project", type=str, help="give the project name", default="")
    anno_parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is show=True; provide no, false, or 0 to block print to screen")
    anno_parser.add_argument("-k", "--key", type=str, help="Choose the key of leiden clustering", default="leiden_0.4")
    anno_parser.add_argument("-r", "--rank_key", type=str, help="Choose the key of rank_genes_groups to be compared to marker_ref, \
    .g., rank_genes_groups", default="rank_genes_groups_r0.6")
    anno_parser.add_argument("-n", "--new_cluster_names", type=str, nargs="+", help="provide the cell type name corresponding to each cluster")

    anno_parser.set_defaults(func=annotate_main)

    ### CONVERT TO H5AD
    convert_parser = subparsers.add_parser('convert_pickle', description="Convert from pickle to h5ad")
    convert_parser.add_argument("input_file", type=str, help="path of the file to convert.", default="after_annotated.pickle")
    convert_parser.set_defaults(func=convert_pickle_main)

    args = parser.parse_args()
    args.func(args)
    

if __name__ == "__main__":
    main()
