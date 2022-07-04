#!/usr/bin/env python

#import h5py
import argparse
import pandas as pd
import scanpy as sc
from collections import namedtuple
from anndata import AnnData
import helpers

"""
Python wrapper method for the markers method
"""
def python_marker(input, output, **kwargs):
    kwargs["input"] = input
    kwargs["out"] = output
    kwargs.update(kwargs['marker_gene']) #collapse 
    args = namedtuple("args", kwargs.keys())(*kwargs.values()) # turn kwargs into an object
    markers(args)

def markers(args):
    """
    This is finding marker genes for each cell cluster.
    """

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
    sc.tl.rank_genes_groups(adata, groupby=groupby, groups=groups, reference=reference, method=method, key_added=key_added)
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


def main():

    parser = argparse.ArgumentParser(fromfile_prefix_chars="@", description="Arguments for scRNA-seq Clustering")
    # basic parameters
    parser.add_argument("-i", "--input", type=str, help="the path of after_leiden.h5ad", default="after_leiden.h5ad")
    parser.add_argument("-d", "--dpi", type=int, help="the resolution of the output figure", default=80)
    parser.add_argument("-f", "--figure_type", type=str, help="the type of plots, e.g., png, pdf, or svg", default="pdf")
    parser.add_argument("-p", "--project", type=str, help="the project name", default="")
    parser.add_argument("-o", "--out", type=str, help="the name of the output anndata", default="after_ranking_gene.h5ad")
    parser.add_argument("-s", "--figsize", type=float, nargs=2, help="the size of the output figure, use 2 numbers, e.g., 2 2")
    parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="block outputing figures on the screen by providing no, false, or 0")

    # rank gene parmeters
    parser.add_argument("-b", "--groupby", type=str, help="the key of the obs grouping to be condisder, e.g., leiden, leiden_0.6", default="leiden")
    parser.add_argument("-n", "--n_genes", type=int, help="the number of top marker genes to plot", default=25)
    parser.add_argument("-N", "--export_num", type=int, help="the number of top marker genes to be exported to csv", default=50)
    parser.add_argument("-m", "--method", type=str, help="the method for differential expression analysis, either t-test, wilcoxon, t-test_overestim_var, logreg", default="wilcoxon")
    parser.add_argument("-c", "--corr_method", type=str, help="the p-value correction method, either benjamini-hochberg or bonferroni", default="benjamini-hochberg")
    parser.add_argument("-R", "--reference", type=str, help="the cell group to be compared with", default="rest")
    parser.add_argument("-g", "--groups", type=str, nargs="+", help="the subset cell groups to compare, using the group names in anndata.obs, e.g. ['g1', 'g2', 'g3'] or ['0', '1', '2']", default="all")
    parser.add_argument("-k", "--key_added", type=str, help="the key in adata.uns indicating where the information to be saved to")

    parser.set_defaults(func=markers)
    args = parser.parse_args()
    args.func(args)
    print()
    print("The arguments are: ", args)
    
if __name__ == "__main__":
    main()
