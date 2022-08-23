import argparse

from collections import namedtuple

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

import helpers
from helpers import MultiLineArgAndFileParser

def annotate_main(args):
    figsize = args.figsize
    print(figsize)
    print(args.out)
    print('asdfasdfasd')

def qc_main(args):
    print("this is qc")

def main():
    # top-level parser
    parser = MultiLineArgAndFileParser()
 
    subparsers = parser.add_subparsers()

    ## ANNOTATE 
    anno_parser = subparsers.add_parser('annotate', fromfile_prefix_chars="@", description="Arguments for annotate cell types for clusters")
    # optional argument
    anno_parser.add_argument("-i", "--input_file", type=str, help="path of the input of 'after_ranking_gene.pickle'", default="after_ranking_gene.pickle")
    anno_parser.add_argument("-m", "--marker_ref_path", type=str, help="path of panglao reference markers", \
        default="../scanpy_scripts/reference_markers/marker_panglao_brain_dic_update.json")
    anno_parser.add_argument("-o", "--out", type=str, help="path of the anndata object to be saved", default="after_annotated.pickle")
    anno_parser.add_argument("-d", "--dpi", type=int, help="resolution of the output figure", default=80)
    anno_parser.add_argument("-s", "--figsize", type=float, nargs=2, help="size of output figure, use 2 numbers, e.g., 2 2")
    anno_parser.add_argument("-f", "--figure_type", type=str, help="define the export type of plot_type, e.g., png, pdf, or svg", default="pdf")
    anno_parser.add_argument("-p", "--project", type=str, help="give the project name", default="")
    anno_parser.add_argument("-S", "--show", type=lambda x: (str(x).lower() in ['true', "1", "yes"]), help="default is show=True; provide no, false, or 0 to block print to screen")
    anno_parser.add_argument("-k", "--key", type=str, help="Choose the key of leiden clustering", default="leiden_0.4")
    anno_parser.add_argument("-r", "--rank_key", type=str, help="Choose the key of rank_genes_groups to be compared to marker_ref, \
        e.g., rank_genes_groups", default="rank_genes_groups_r0.6")
    anno_parser.add_argument("-n", "--new_cluster_names", type=str, nargs="+", help="provide the cell type name corresponding to each cluster")

    anno_parser.set_defaults(func=annotate_main)
    
    qc_parser = subparsers.add_parser('qc')
    qc_parser.add_argument("--hello", type=str, help="test")
    qc_parser.set_defaults(func=qc_main)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

    
