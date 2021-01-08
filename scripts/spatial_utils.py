#!/bin/bash
#############################
# File Name : spatial_utils.py
#
# Purpose : [???]
#
# Creation Date : 06-01-2021
#
# Last Modified : Wed 06 Jan 2021 09:08:27 AM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

def spatial_zscore_polarity(spots,cells):
    #Calculate gene centroid within each cell
    gene_centroids = spots.groupby(['cell_id','target_molecule_name'])[['global_x','global_y']].mean().reset_index()
    gene_centroids = gene_centroids.rename(columns={'global_x':'gene_centroidX','global_y':'gene_centroidY'})

    #Calculate distance of gene_centroid to cell_centroid
    gene_centroids = add_cell_centroids_to_spots(cells, gene_centroids)
    gene_centroids['gene_centroid_dist'] = np.sqrt(
        ((gene_centroids['cell_centroidX']-gene_centroids['gene_centroidX'])**2) +
        ((gene_centroids['cell_centroidY']-gene_centroids['gene_centroidY'])**2)
    )
 
    gene_centroids['z_score_polarity'] = zscore_normalize_column(gene_centroids, 'gene_centroid_dist')

    return gene_centroids
    

def spatial_zscore_centroid_dist(spots,cells):

    #For each spot, calculate distance to cell centroid
    spots = add_cell_centroids_to_spots(cells, spots)
    spots['centroid_dist'] = np.sqrt(
        ((spots['cell_centroidX']-spots['global_x'])**2) +
        ((spots['cell_centroidY']-spots['global_y'])**2)
    )
   
    #Calculate mean centroid_dist for each gene in each cell 
    gene_scores = spots.groupby(['cell_id','target_molecule_name'])['centroid_dist'].agg(['count','mean']).reset_index()
    gene_scores = gene_scores.rename(columns={'count':'num_cell_spots','mean':'mean_centroid_dist'})

    #Center and normalize scores 
    gene_scores['z_score_centroid_dist'] = zscore_normalize_column(gene_scores, 'mean_centroid_dist')
    
    return gene_scores


def calculate_cell_centroids(cells):
    """
    Add columns for 'centroidX' and 'centroidY' to a cells dataframe

    Requires having columns boundaryX and boundaryY

    Just averages all the X boundaries, and then all the Y boundaries
    """
    if 'centroidX' in cells.columns and 'centroidY' in cells.columns:
        return cells

    cells['centroidX'] = cells['boundaryX'].str.split(', ').apply(lambda xs: sum(float(x) for x in xs)/len(xs))
    cells['centroidY'] = cells['boundaryY'].str.split(', ').apply(lambda ys: sum(float(y) for y in ys)/len(ys))

    return cells

def add_cell_centroids_to_spots(cells, spots):
    #Add cell centroid info to the spots table
    cells = calculate_cell_centroids(cells)
    cell_id_to_centroidX = dict(cells[['cell_id','centroidX']].values)
    cell_id_to_centroidY = dict(cells[['cell_id','centroidY']].values)
    spots['cell_centroidX'] = spots['cell_id'].map(cell_id_to_centroidX)
    spots['cell_centroidY'] = spots['cell_id'].map(cell_id_to_centroidY)

    return spots

def zscore_normalize_column(df, col_name):
    m = df[col_name].mean()
    var = df[col_name].std()**2
    return (df[col_name]-m)/var
 

