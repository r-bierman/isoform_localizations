#!/bin/bash
#############################
# File Name : summarize_perms.py
#
# Purpose : [???]
#
# Creation Date : 28-01-2021
#
# Last Modified : Fri 29 Jan 2021 02:07:04 PM PST
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
   
    overall_Q25 = gene_summary_info['normalized_rank'].quantile(0.25)
    overall_Q50 = gene_summary_info['normalized_rank'].quantile(0.50)
    overall_Q75 = gene_summary_info['normalized_rank'].quantile(0.75)

    ret_df = pd.DataFrame({
        'gene':[gene_name],
        'num_cells':[perm_tab.shape[0]],
        'real_Q25':[overall_Q25],
        'real_Q50':[overall_Q50],
        'real_Q75':[overall_Q75],
    })

    ret_df['overall_clusts_lt_count_Q25'] = (perm_tab.quantile(0.25) < overall_Q25).sum()
    ret_df['overall_clusts_lt_count_Q50'] = (perm_tab.quantile(0.50) < overall_Q50).sum()
    ret_df['overall_clusts_lt_count_Q75'] = (perm_tab.quantile(0.75) < overall_Q75).sum()

    ret_df['overall_clusts_Q25'] = perm_tab.quantile(0.25).median()
    ret_df['overall_clusts_Q50'] = perm_tab.quantile(0.50).median()
    ret_df['overall_clusts_Q75'] = perm_tab.quantile(0.75).median()


    for clust_id,g in gene_summary_info.groupby('clust_id'):
        clust_Q25 = g['normalized_rank'].quantile(0.25)
        clust_Q50 = g['normalized_rank'].quantile(0.50)
        clust_Q75 = g['normalized_rank'].quantile(0.75)
        
        ret_df['clust_{}_cell_count'.format(clust_id)] = perm_tab.loc[g.index].shape[0]
        
        ret_df['clust_{}_lt_count_Q25'.format(clust_id)] = (perm_tab.loc[g.index].quantile(0.25) < clust_Q25).sum()
        ret_df['clust_{}_lt_count_Q50'.format(clust_id)] = (perm_tab.loc[g.index].quantile(0.50) < clust_Q50).sum()
        ret_df['clust_{}_lt_count_Q75'.format(clust_id)] = (perm_tab.loc[g.index].quantile(0.75) < clust_Q75).sum()
        
        ret_df['clust_{}_Q25'.format(clust_id)] = perm_tab.loc[g.index].quantile(0.25).median()
        ret_df['clust_{}_Q50'.format(clust_id)] = perm_tab.loc[g.index].quantile(0.50).median()
        ret_df['clust_{}_Q75'.format(clust_id)] = perm_tab.loc[g.index].quantile(0.75).median()
 
    ret_df.to_csv(
        '/scratch/PI/horence/rob/isoform_localizations/perm_by_gene/{}_perms_over_cells.csv'.format(gene_name),
        index=False,
    ) 

if __name__ == '__main__':
    
    f_path = sys.argv[1]
    run(f_path)

