#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-07-22'
__version__ = '0.0.1'

import argparse
import warnings
import random
from distutils.version import LooseVersion
import os
import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc

# avoid tk.Tk issues
# import matplotlib
# matplotlib.use('Agg')


# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Differential cell type composition.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-cond', '--condition_column',
        action='store',
        dest='condition_column',
        required=True,
        help='Condition column. Should be a binary trait.'
    )

    parser.add_argument(
        '-exper', '--experiment_id_column',
        action='store',
        dest='experiment_id_column',
        default='experiment_id',
        help='Column with sample / experiment id groupings. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-cl', '--cell_label_column',
        action='store',
        dest='cell_label_column',
        default='cluster',
        help='Anndata cell type label name in obs slot. (default: %(default)s)'
    )

    parser.add_argument(
        '-cla', '--cell_label_analyse',
        action='store',
        dest='cell_label_analyse',
        default='',
        help='Only analyse cells with these labels. Either "all" or comma \
            seperated list of cell labels. (default: %(default)s)'
    )

    # parser.add_argument(
    #     '-m', '--method',
    #     action='store',
    #     dest='method',
    #     default='wald',
    #     help='Differential expression method to use. (default: %(default)s)'
    # )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='diffxpy',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    experiment_id_column = options.experiment_id_column
    condition_column = options.condition_column
    cell_label_column = options.cell_label_column
    cell_label_analyse = options.cell_label_analyse
    out_file_base = options.of

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # List to stash continuous variables.
    continuous_variables = []

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    # Sort out the condition we want to test
    adata.obs['condition'] = adata.obs[condition_column]
    # Check to see if the condition is a categorical variable
    # if adata.obs[condition_column].dtype.name != 'category':
    #     # raise Exception('Condition is not a category.')
    #     continuous_variables.append('condition')
    #     warnings.warn('Treating condition as a continuous variable')
    # The below code converts the condition into an int and asserts it binary
    adata.obs['condition'] = adata.obs[condition_column].cat.codes
    if len(np.unique(adata.obs['condition']) != 2):
        raise Exception('There are not exactly 2 conditions.')

    # Select the subset of cells we want for analysis
    if cell_label_analyse != 'all':
        cell_label_analyse = cell_label_analyse.split(',')
        adata = adata[adata.obs[cell_label_column].isin(cell_label_analyse)]
    # clusters = np.sort(adata.obs[cell_label_column].unique())

    # Check to make sure that within each cluster, there are >1 condition
    # values.
    if adata.obs[condition_column].dtype.name == 'category':
        n_cells_condition_cluster = adata.obs.groupby(
            [condition_column, cell_label_column]
        ).size().unstack()
        if n_cells_condition_cluster.values.min() == 0:
            raise Exception(
                'For one cell_label_column there are 0 conditions.'
            )
        if len(np.unique(adata.obs[condition_column].cat.codes)) <= 1:
            raise Exception('There is only 1 condition.')

    print('Continuous varibles: {}'.format(','.join(continuous_variables)))

    # TODO: check for nan in covariates and conditions?

    # Get the number of cell counts for each cluster and experiment
    cell_stats = adata.obs.groupby(
        [experiment_id_column, cell_label_column, 'condition']
    ).size().reset_index(name='number_of_cells')
    # Get the total number of cells per sample
    number_cells_per_experiment = adata.obs.groupby(
        [experiment_id_column]
    ).size()
    # Get the fraction of cells per sample
    cell_stats['fraction_of_cells'] = cell_stats[
        'number_of_cells'
    ] / number_cells_per_experiment.loc[
        cell_stats[experiment_id_column]
    ].values

    # For each cluster:
    # 1. Get the number of cells per sample (experiment) and perform
    #    Fisher’s exact test.
    # 2. Get the fraction of cells per sample (experiment) and perform
    #    Mann-Whitney test.
    #
    # Init what will become results dataframe
    results_list = []
    # Iterate over each cluster
    for group_name, df in cell_stats.groupby(
        [cell_label_column]
    ):
        # Perform Fisher’s exact test on number cells for each experiment
        # split by condition
        # ERROR: python does not supprot > 2x2 tables.
        # oddsratio, pvalue = sp.stats.fisher_exact(
        #     table=[[8, 2], [1, 5]],
        #     alternative='two-sided'
        # )
        # results_list.append({
        #     'cell_label': group_name,
        #     'oddsratio': oddsratio,
        #     'pvalue': pvalue,
        #     'method': 'fisher_exact',
        #     'data_type': 'number_of_cells',
        #     'condition': condition_column,
        #     'condition_effect': '0'
        # })

        # TODO: Perform Mann-Whitney test on the fraction of cells for each
        # experiment split by condition
        foo = sp.stats.mannwhitneyu()

        # Append stats

    # Make a dataframe of the results
    # df_results = test_results.summary()
    # df_results['de_method'] = options.method
    # df_results['condition'] = condition_column
    # df_results['covariates'] = ','.join(covariate_columns)
    # df_results['cell_label_column'] = cell_label_column
    # df_results['cell_label_analysed'] = ','.join(cell_label_analyse)
    # df_results = df_results.sort_values(
    #     by=['pval', 'log2fc', 'mean'],
    #     ascending=[True, False, False]
    # )
    df_results = []
    df_results.to_csv(
        '{}-de_results.tsv.gz'.format(out_file_base),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )

    # Make a stacked barplot of the number of cells per each sample

    # Make a barplot of the fraction of cells per each sample


if __name__ == '__main__':
    main()
