#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-08-17'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import warnings
from patsy import ModelDesc
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import scanpy as sc


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read an h5ad object and write counts matrix (row = sample_id and
            column = cell type). Also write metadata on samples.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5ad', '--h5ad_file',
        action='store',
        dest='h5ad',
        required=True,
        help='H5anndata file.'
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
        default='all',
        help='Only analyse cells with these labels. Either "all" or comma \
            seperated list of cell labels. (default: %(default)s)'
    )

    parser.add_argument(
        '-mdc', '--metadata_columns',
        action='store',
        dest='metadata_columns',
        default='all',
        help='Comma-separated list of metadata columns. \
            If all, write all cols. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-mdf', '--metadata_formula',
        action='store',
        dest='metadata_formula',
        default='',
        help='Right hand side of a formula. Variables will be paresed and \
            added to metadata_columns. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='anndata',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    experiment_id_column = options.experiment_id_column
    cell_label_column = options.cell_label_column
    cell_label_analyse = options.cell_label_analyse
    metadata_columns = options.metadata_columns
    metadata_formula = options.metadata_formula
    out_file_base = options.of

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Load the AnnData file and set matrix to be log data
    adata = sc.read_h5ad(filename=options.h5ad)
    
    # Get the metadata columns
    if metadata_columns != 'all':
        metadata_columns = metadata_columns.split(',')
    else:
        metadata_columns = adata.obs.columns.to_list()

    # Parse the formula columns and add to metadata.
    if metadata_formula != '':
        model_desc = ModelDesc.from_formula(metadata_formula)
        for term in model_desc.rhs_termlist:
            for i in term.factors:
                if i not in metadata_columns:
                    metadata_columns.append(i)

    # Add power to data.
    for i in metadata_columns:
        if '__power' in i:
            term = i.split('__power')[0]
            exponent = float(i.split('__power')[1])
            adata.obs[i] = adata.obs[term] ^ exponent

    # Add in key other columns
    if experiment_id_column not in metadata_columns:
        metadata_columns.append(experiment_id_column)
    # if cell_label_column not in metadata_columns:
    #     metadata_columns.append(cell_label_column)

    # Select the subset of cells we want for analysis
    if cell_label_analyse != 'all':
        cell_label_analyse = cell_label_analyse.split(',')
        adata = adata[adata.obs[cell_label_column].isin(cell_label_analyse)]

    # Now generate the cell type counts matrix
    # Make long dataframe of experiment_id cluster number_of_cells
    df_cell_counts = adata.obs.groupby(
        [experiment_id_column, cell_label_column]
    ).size().reset_index(name='number_of_cells')
    df_cell_counts = df_cell_counts.pivot(
        index=experiment_id_column,
        columns=cell_label_column,
        values='number_of_cells'
    )
    df_cell_counts_cols = df_cell_counts.columns.values
    if is_numeric_dtype(adata.obs[cell_label_column]):
        df_cell_counts.columns = [
            '{}__{}'.format(cell_label_column, i) for i in df_cell_counts_cols
        ]
    # Save the cell type counts matrix
    df_cell_counts.to_csv(
        '{}-counts.tsv.gz'.format(out_file_base),
        sep='\t',
        compression=compression_opts,
        index=True,
        header=True
    )

    # Check for per cell metrics (rather than per sample) and emit
    # warning if any metric is a per cell metric.
    # Summarize by taking the mean per sample.
    n_values = adata.obs.groupby(
        experiment_id_column
    )[metadata_columns].nunique()
    col_val_counts_per_sample = n_values != 1
    col_val_counts = col_val_counts_per_sample.any(axis='rows')
    if col_val_counts.any():
        cols_not_uniq = col_val_counts[col_val_counts].index.tolist()
        warnings.warn(
            'Some columns do not reduce down to a single value:\t{}.'.format(
                 ','.join(cols_not_uniq)
            )
        )
        # Replace these metrics with their mean
        # Calculate mean per experiment_id_column value
        cols_not_uniq_means = adata.obs.groupby(
            experiment_id_column
        )[cols_not_uniq].mean()
        # Expand means to match adata.obs
        new_data = cols_not_uniq_means.loc[adata.obs[experiment_id_column], :]
        # Set index to cellular barcode
        new_data.index = adata.obs.index
        # Delete old values
        adata.obs = adata.obs.drop(cols_not_uniq, axis=1)
        # Write new dataframe
        adata.obs = pd.concat([adata.obs, new_data], axis=1)

        # NOTE: The above operation will not handle the case where there are
        #       non-numeric values that vary within a sample (e.g., cell type
        #       or phase). Drop those variables from metadata_columns.
        a = set(cols_not_uniq)
        b = set(adata.obs.columns.to_list())
        diff = a.difference(b)
        if len(diff) > 0:
            warnings.warn(
                '{}:\t{}.'.format(
                    'Dropping columns that could not be averaged across cells',
                    ','.join(diff)
                )
            )
            for col in diff:
                metadata_columns.remove(col)

    # Now compress the metadata matrix
    # The below code keeps multiple rows per sample
    # df_metadata = adata.obs.set_index(experiment_id_column, drop=True)
    #
    # If there are per cell metrics, these are ignored and only the first
    # value for a sample is saved. Example if score per cell caclulated like
    # villus__bottom_enterocytes... then the score for whatever the first
    # cell in adata.obs is kept, the rest are discarded.
    df_metadata = adata.obs[metadata_columns].drop_duplicates(
        subset=experiment_id_column,
        keep='first'
    )
    df_metadata = df_metadata.set_index(experiment_id_column, drop=True)
    # Save the metadata matrix
    df_metadata.loc[df_cell_counts.index, :].to_csv(
        '{}-metadata.tsv.gz'.format(out_file_base),
        sep='\t',
        compression=compression_opts,
        index=True,
        header=True
    )


if __name__ == '__main__':
    main()
