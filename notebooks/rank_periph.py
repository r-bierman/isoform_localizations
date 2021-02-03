#!/bin/bash
#############################
# File Name : rank_periph.py
#
# Purpose : [???]
#
# Creation Date : 20-01-2021
#
# Last Modified : Wed 20 Jan 2021 04:08:34 PM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import numpy as np
import glob
import sys
import os

import matplotlib.pyplot as plt
import seaborn as sns

root_dir = '/oak/stanford/groups/horence/rob/isoform_localizations/'

sys.path.append(os.path.join(root_dir,'scripts'))
import spatial_utils
import plot_utils

spots = pd.read_csv('../processed_data/min_periph_dists.csv') #normal spots table with added info

#Use the previously subsamples cells in the q90
cells = pd.read_csv('../processed_data/q90_cells.csv')


#Filter the spots down to just the ones that match the q90 cells
shared_cell_ids = np.intersect1d(spots['cell_id'], cells['cell_id'])
spots = spots[spots['cell_id'].isin(shared_cell_ids)]
cells = cells[cells['cell_id'].isin(shared_cell_ids)]

#add the min dists into the spots df and call it "raw_metric"
spots = spots.rename(columns={'min_boundary_dist':'raw_metric'})

#Calculate num_cell_spots and num_gene_spots
spots = spatial_utils.calculate_spots_per_cell(spots)
print('Calculated spots per cell')


#Calculate normalized spot ranks
spots['raw_metric'] = (
    spots.groupby('cell_id')['raw_metric']
    .transform(
        lambda x: [(x <= v).sum()/len(x) for v in x]
    )
)
print('Calculated normalized spot rank')

periphery_df = (
    spots.groupby(
        ['cell_id','target_molecule_name','num_cell_spots','num_gene_spots']
    )['raw_metric']
    .quantile(0.5) #currently just the median
    .reset_index()
)
print('Calculated median gene normalized spot rank')

periphery_df['metric_name'] = 'periphery_ranks'

#Get the output into the standard metric table
cols = [
   'cell_id', 'target_molecule_name',
   'num_cell_spots', 'num_gene_spots',
   'metric_name', 'raw_metric',
]
metric_df = periphery_df[cols]
metric_df.to_csv('../processed_data/rank_periph.csv',index=False)

