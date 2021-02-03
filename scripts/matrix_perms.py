#!/bin/bash
#############################
# File Name : permutations.py
#
# Purpose : [???]
#
# Creation Date : 27-01-2021
#
# Last Modified : Thu 28 Jan 2021 07:24:07 AM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import numpy as np
import glob
import time
import sys
import os

import matplotlib.pyplot as plt
import seaborn as sns


def matrix_approach(test_df,gene_name):
    cells = test_df[test_df['target_molecule_name'].eq(gene_name)]['cell_id'].unique()

    test_df = test_df[test_df['cell_id'].isin(cells)]

    cell_spot_counts = test_df.groupby('cell_id').size()
    cell_gene_spot_counts = test_df[test_df['target_molecule_name'].eq(gene_name)].groupby(['cell_id']).size()
    
    num_cells = cells.size
    num_its = 1000

    ret_df = np.zeros((num_cells, num_its))
    ret_df[:,:] = np.NaN
    
    for i,cell in enumerate(cells):
        
        rand_ranks = np.random.choice(
            cell_spot_counts[cell],
            (num_its, cell_gene_spot_counts[cell])
        )/cell_spot_counts[cell]
        
        ret_df[i,:] = np.median(rand_ranks, axis=1)
            
    ret_df = pd.DataFrame(ret_df)
    ret_df.index = cells
    return ret_df


if __name__ == '__main__':
    gene = sys.argv[1]
    print('|'+gene+'|')

    #Read in all the min_periph_dists that are pre-calculated
    q90_df = pd.read_csv('../processed_data/q90_min_periph_dists.csv')

    start = time.time()
    x = matrix_approach(q90_df,gene)
    print(time.time()-start)
    x.to_csv('/scratch/PI/horence/rob/isoform_localizations/perm_by_gene/test_perms_matrix_{}.csv'.format(gene))


