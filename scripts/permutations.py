#!/bin/bash
#############################
# File Name : permutations.py
#
# Purpose : [???]
#
# Creation Date : 27-01-2021
#
# Last Modified : Wed 27 Jan 2021 11:26:09 AM PST
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


def permutations(sub_df, its=100):
    def gene_perm_median_scores(sub_df):
        for cell_id,cell_g in sub_df.groupby('cell_id'):
            num_spots = cell_g.shape[0]
            cell_g['norm_ranks'] = np.random.permutation(np.arange(num_spots))/num_spots
            gene_medians = cell_g.groupby('target_molecule_name')['norm_ranks'].median()

            yield gene_medians.rename(cell_id)

    for i in range(its):
        single_perm_all_cells = pd.concat(gene_perm_median_scores(sub_df),axis=1)
        single_perm_medians = single_perm_all_cells.median(axis=1)
        yield single_perm_medians.rename(i)



#Read in all the min_periph_dists that are pre-calculated
df = pd.read_csv('../processed_data/min_periph_dists.csv')
cells = pd.read_csv('../processed_data/q90_cells.csv')

#Limit to just run on cells that were in the Q90 with the most spots
q90_df = df[df['cell_id'].isin(cells['cell_id'].unique())]

start = time.time()
x = pd.concat(permutations(q90_df, its=1),axis=1)
print(time.time()-start)
x.to_csv('gene_test_perms.csv')


