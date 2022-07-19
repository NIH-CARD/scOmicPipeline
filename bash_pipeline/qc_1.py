#!/usr/bin/env python3

import argparse
import pandas as pd
import scanpy as sc
import helpers
from helpers import MyArgumentParser

#from bash_pipeline import helpers

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

def main():
    parser = MyArgumentParser(fromfile_prefix_chars="@", description="Arguments for scRNA-seq QC")
    # basic parameters
    parser.add_argument("-p", "--project", type=str, help="a project name", default="")
    parser.add_argument("-i", "--input_file", type=str, help="the path of the input .h5 file or the matrix folder", default="filtered_feature_bc_matrix.h5")
    parser.add_argument("-o", "--out", type=str, help="the file name prefix of the output anndata", default="count_after_QC.pickle")
    parser.add_argument("-d", "--dpi", type=int, help="the resolution of the output figure", default=80)
    parser.add_argument("-s", "--figsize", type=float, nargs=2, help="the size of the output figure, use 2 numbers, e.g., 2 2", default=None)
    parser.add_argument("-f", "--figure_type", type=str, help="the type of the output figure, e.g., pdf, png, or svg", default="pdf")
    parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="block output figures on the screen by providing no, false, or 0")
    # filtering and plotting parameters
    parser.add_argument("-g", "--min_genes", type=int, help="cell filter: the minimal number of genes which a cell should have", default=200)
    parser.add_argument("-c", "--min_cells", type=int, help="gene filter: the minimal number of cells which a gene should be in", default=3)
    parser.add_argument("-n", "--n_genes_by_counts", type=int, help="the threthold of gene counts", default=8000)
    parser.add_argument("-m", "--pct_counts_mt", type=float, help="the threthold of mitochondrial genes percentage", default=5.0)
    parser.add_argument("-t", "--n_top", type=int, help="the number of genes to plot in the highest expressed gene plot", default=20)
    parser.add_argument("-C", "--color_gene", type=str, help="the gene whoes expression level is to be colored in the pca plot", default="MAP2")
    parser.add_argument("-e", "--exclude_highly_expressed", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="exclude highly expressed genes when do normalization by providing yes, no, or 1 ")

    args = parser.parse_args()

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

if __name__ == "__main__":
    main()
