#!/bin/bash
#############################
# File Name : calculate_periph_dists.py
#
# Purpose : [???]
#
# Creation Date : 20-01-2021
#
# Last Modified : Wed 20 Jan 2021 10:51:56 AM PST
#
# Created By : Rob Bierman
#
##############################

import spatial_utils
import pandas as pd

print('Reading in spots and cells')
spots = pd.read_csv('/oak/stanford/groups/horence/rob/isoform_localizations/processed_data/all_passing_rna_spots.csv')
cells = pd.read_csv('/oak/stanford/groups/horence/rob/isoform_localizations/processed_data/segmented_cell_shapes.csv')

print('Running distance calculation')
res = spatial_utils.calculate_spot_to_boundary_min_dist(spots, cells)

print('Saving results')
res.to_csv('/oak/stanford/groups/horence/rob/isoform_localizations/processed_data/min_periph_dists.csv',index=False)

