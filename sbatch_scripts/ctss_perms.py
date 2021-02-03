#!/bin/bash
#############################
# File Name : ctss_perms.py
#
# Purpose : [???]
#
# Creation Date : 26-01-2021
#
# Last Modified : Tue 26 Jan 2021 08:49:08 PM PST
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

root_dir = '/oak/stanford/groups/horence/rob/isoform_localizations/'

sys.path.append(os.path.join(root_dir,'scripts'))
import spatial_utils
import plot_utils

import time

start = time.time()
df = pd.read_csv('../processed_data/min_periph_dists.csv')
print('Done with data read in:',time.time()-start)

start = time.time()
gene = 'Ctss'
total_spots = df.groupby('cell_id').size().astype(int)
gene_spots = df[df['target_molecule_name'].eq(gene)].groupby('cell_id').size().astype(int)

gene_counts = pd.concat((gene_spots,total_spots),axis=1).dropna()
gene_counts.columns = ['num_gene_spots','num_total_spots']
gene_counts['num_gene_spots'] = gene_counts['num_gene_spots'].astype(int)
print('Done with counts calculation:',time.time()-start)


start = time.time()
its = 5000
for i in range(its):
    gene_counts['it_{}'.format(i)] = gene_counts.apply(
        lambda r: np.median(
            np.random.choice(
                np.arange(int(r['num_total_spots']))/r['num_total_spots'], 
                int(r['num_gene_spots']))
        ),
        axis=1,
    )

print('Done with perms after:',time.time()-start)

gene_counts.to_csv('ctss_perms.csv')

