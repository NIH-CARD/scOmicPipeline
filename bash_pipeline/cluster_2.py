import h5py
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import helpers
from collections import namedtuple
from anndata import AnnData
import pickle


"""
Wrapper method for cluster
It also deletes the previous, no longer necessary, pickle of the object.
"""
def python_cluster(input, output, **kwargs):
    kwargs["input"] = input
    kwargs["out"] = output
    args = namedtuple("args", kwargs.keys())(*kwargs.values()) # turn kwargs into an object
    cluster(args)

def cluster(args):
    """
    Clustering cells after computing pca and neiborhood distances.
    """
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

def main():
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@", description="Arguments for scRNA-seq Clustering")
    # basic parameters
    parser.add_argument("-i", "--input", type=str, help="the path of count_after_QC.pickle file", default="count_after_QC.pickle")
    parser.add_argument("-d", "--dpi", type=int, help="the resolution of the output figure", default=80)
    parser.add_argument("-f", "--figure_type", type=str, help="the export type of plots, e.g., png, pdf, or svg", default="pdf")
    parser.add_argument("-p", "--project", type=str, help="the project name", default="")
    parser.add_argument("-o", "--out", type=str, help="the file name to save the anndata object", default="after_leiden.pickle")
    parser.add_argument("-s", "--figsize", type=float, nargs=2, help="the size of output figure, use 2 numbers, e.g., 2 2")
    parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="block output figures on the screen by providing no, false, or 0")
    
    # umap parmeters
    parser.add_argument("-n", "--n_neighbors", type=int, help="the size of local neiborhood for manifold approximation", default=15)
    parser.add_argument("-P", "--n_pcs", type=int, help="the number of PCs to use", default=None)
    parser.add_argument("-m", "--method", type=str, help="the method for neighborhood graph, either ‘umap’, ‘gauss’, ‘rapids’", default="umap")
    parser.add_argument("-M", "--metric", type=str, help="the metric for neighborhood graph, [‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’], Literal[‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’],", default="euclidean")

    # leiden parameters
    parser.add_argument("-r", "--resolution", type=float, help="the resolution for leiden", default=1.0)

    # color parameters and key names to be stored in adata
    parser.add_argument("-C", "--color_gene", type=str, nargs="*", help="define a list of genes (e.g., MAP2 TEME199 TMEM106B), a key of leiden (e.g., 'leiden' or other key_added like 'leiden_0.6'), or both as color hues in umap plot", default="leiden")
    # parser.add_argument("-g", "--gene_list", type=str, nargs="+", action="store", dest="list", help="define a list of genes to show in umap, e.g., MAP2 TEME199 NIL", default=['leiden'])
    parser.add_argument("-k", "--key_added", type=str, help="the key name of a ledien anaysis to be addeed to anndata", default='leiden')

    parser.set_defaults(func=cluster)
    args = parser.parse_args()
    args.func(args)
    print()
    print(f"The arguments are {args}")

if __name__ == "__main__":
    main()
