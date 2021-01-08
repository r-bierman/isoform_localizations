#!/bin/bash
#############################
# File Name : all_centroids.py
#
# Purpose : [???]
#
# Creation Date : 06-01-2021
#
# Last Modified : Wed 06 Jan 2021 10:17:50 AM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import sys
import os

import matplotlib.pyplot as plt
import seaborn as sns

root_dir = '/oak/stanford/groups/horence/rob/MERFISH_spatial_data'

sys.path.append(os.path.join(root_dir,'scripts'))
import spatial_utils
import plot_utils

spots = pd.read_csv(os.path.join(root_dir,'processed_data','all_passing_rna_spots.csv'))
cells = pd.read_csv(os.path.join(root_dir,'processed_data','segmented_cell_shapes.csv'))

y = spatial_utils.spatial_zscore_centroid_dist(spots, cells)
y.to_csv('somehow_finished_centroids.csv',index=False)

