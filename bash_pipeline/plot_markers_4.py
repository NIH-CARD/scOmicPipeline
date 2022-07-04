#!/usr/bin/env python3

#import h5py
import argparse
import pandas as pd
import scanpy as sc
import helpers
from collections import namedtuple

def setup(args):
    """
    Set up the displaying and output parameters
    """
    sc.settings.verbosity = 3 
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=args.dpi, facecolor="white", figsize=args.figsize) # frameon=False
    sc.logging.print_version_and_date()

"""
Wrapper method for plot_makers
"""
def python_plot_makers(input, **kwargs):
    kwargs["input"] = input
    kwargs.update(kwargs['plot_makers']) #collapse 
    args = namedtuple("args", kwargs.keys())(*kwargs.values()) # turn kwargs into an object
    plot_makers(args)


def plot_makers(args):
    """
    Visualize the differential expression of marker genes across clusters.
    """
    input_file = args.input_file
    dpi = args.dpi
    figsize = args.figsize
    figure_type = args.figure_type
    show = args.show
    project = args.project if (args.project == "") else ("_" + args.project)
    use_raw = args.use_raw

    groupby = args.groupby
    n_genes = args.n_genes
    groups = args.groups
    plot_type = args.plot_type
    gene_list = args.gene_list
    dpi = args.dpi
    series = args.series
    adata = helpers.choose_adata_src(input_file) 

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

def main():
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@", description="Arguments for scRNA-seq plotting")
    # basic parameters
    parser.add_argument("-i", "--input_file", type=str, help="path of the input after_leiden.h5ad file", default="after_ranking_gene.h5ad")
    parser.add_argument("-d", "--dpi", type=int, help="resolution of the output figure", default=None)
    parser.add_argument("-s", "--figsize", type=float, nargs=2, help="size of output figure, use 2 numbers, e.g., 2 2")
    parser.add_argument("-f", "--figure_type", type=str, help="define the export type of plot_type, e.g., png, pdf, or svg", default="pdf")
    parser.add_argument("-p", "--project", type=str, help="give the project name", default="")
    parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is show=True; provide no, false, or 0 to block print to screen")
    parser.add_argument("-r", "--use_raw", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is using raw; provide no, false, or 0 to not use raw data")

    # plotting gene parmeters
    parser.add_argument("-g", "--gene_list", type=str, nargs="*", help="define a list of genes (e.g., MAP2 TEME199 TMEM106B) for plotting", default=["MAP2"])
    parser.add_argument("-t", "--plot_type", type=str, nargs="*", help="define the plotting types, e.g., dotplot, violin, stacked_violin, and rank_genes_groups_violin", default=["violin", "dotplot", "stacked_violin", "rank_genes_groups_violin", "umap"])
    parser.add_argument("-b", "--groupby", type=str, help="the key of the obs grouping to be considered, e.g., leiden, leiden_0.6")
    #parser.add_argument("-b", "--groupby", type=str, help="the key of the obs grouping to be condisder, e.g., leiden, leiden_0.6", default="leiden")
    parser.add_argument("-n", "--n_genes", type=int, help="number of genes used for plotting rank_genes_groups_violin", default=8)
    parser.add_argument("-G", "--groups", type=str, nargs="*", help="choose the subset cell groups for plotting rank_genes_groups_violin, e.g. ['g1', 'g2', 'g3'] or ['0', '1', '2']")
    parser.add_argument("-X", "--series", type=str, help="A name to append after the project name for output files.")


    parser.set_defaults(func=plot_makers)
    args = parser.parse_args()
    args.func(args)



    print("\nThe arguments are: ", args)

    setup(args)



if __name__ == "__main__":
    main()
