#!/bin/bash
#############################
# File Name : summarize_perms.py
#
# Purpose : [???]
#
# Creation Date : 28-01-2021
#
# Last Modified : Thu 28 Jan 2021 04:49:38 PM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import numpy as np
import sys
import os

def run(f_path):
    full_path = os.path.join('/scratch/PI/horence/rob/isoform_localizations/perm_by_gene',f_path)
    gene_cell_df = pd.read_csv('../processed_data/20210127_q90_cell_gene_med_norm_ranks.csv')
    gene_name = os.path.basename(full_path).split('_')[-1].split('.')[0]
    
    perm_tab = pd.read_csv(full_path, index_col=0)
    
    gene_summary_info = (
        gene_cell_df[
            gene_cell_df[
                'target_molecule_name'].eq(gene_name)
        ].set_index('cell_id')
    )
    
    perm_summary = pd.DataFrame({
        'num_more_extreme_perms':perm_tab.sub(
            gene_summary_info['normalized_rank'],axis=0
        ).lt(0).sum(axis=1),
        'Q25_perm':perm_tab.quantile(0.25, axis=1),
        'Q50_perm':perm_tab.quantile(0.50, axis=1),
        'Q75_perm':perm_tab.quantile(0.75, axis=1),
    })

    collected_info = pd.concat((perm_summary,gene_summary_info),axis=1)
    collected_info.to_csv('/scratch/PI/horence/rob/isoform_localizations/perm_by_gene/{}_collected.csv'.format(gene_name))


if __name__ == '__main__':
    
    f_path = sys.argv[1]
    run(f_path)

