#!/bin/bash
#############################
# File Name : spatial_utils.py
#
# Purpose : [???]
#
# Creation Date : 06-01-2021
#
# Last Modified : Wed 13 Jan 2021 07:48:01 AM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import numpy as np
import argparse

pd.options.mode.chained_assignment = None  # default='warn'

import functools

def metric_summaries(func):
    """
    Decorator for metric functions
    
    Calculates cell_zscore, overall_zscore, etc
    The decorated function must have the following columns
    * cell_id
    * raw_metric
    * target_molecule_name
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        df = func(*args, **kwargs)
        
        cell_metric_means = df.groupby('cell_id')['raw_metric'].transform('mean')
        cell_metric_vars = df.groupby('cell_id')['raw_metric'].transform('std')**2
        df['cell_zscore'] = (df['raw_metric']-cell_metric_means)/cell_metric_vars

        df['gene_mean_cell_zscore'] = df.groupby('target_molecule_name')['cell_zscore'].transform('mean')
        df['gene_var_cell_zscore'] = df.groupby('target_molecule_name')['cell_zscore'].transform('std')**2
        
        return df
    
    return wrapper
        

@metric_summaries
def spatial_metric_polarity(spots,cells):
    """
    Spatial metric: polarity

    Calculates the distance between the gene centroid and the cell centroid
    Intuition is that genes that don't have polar distributions will have a gene centroid near the cell centroid

    Arguments:
        - spots: dataframe object where each row is a MERFISH spot with the following columns
            'cell_id': the id of the cell that contains this spot as a string
            'target_molecule_name': gene name
            'global_x': x-coordinate of the spot
            'global_y': y-coordinate of the spot
            

        - cells: dataframe object where each row is a MERFISH cell with the following columns
            'cell_id': the ide of the cell as a string
            'boundaryX': string of x-coords for the cell boundary like
                "2097.1362953431903, 2097.1362953431903, 2097.1362953431903, ...."
            'boundaryY': string of y-coords for the cell boundary like
                "2097.1362953431903, 2097.1362953431903, 2097.1362953431903, ...."


    Outputs:
        - metric_df: dataframe object where each row is a unique gene/cell combination. columns:
            'cell_id': same as from inputs
            'target_molecule_name': same as from inputs
            'num_cell_spots': number of total spots in this cell
            'num_gene_spots': number of spots of the specific gene in this cell
            'metric_name': always 'polarity' for this function
            'raw_metric': the raw metric value of gene-centroid dist to cell centroid
            any/all processed metrics from the metric_summaries wrapper
    """
    #Calculate num_cell_spots and num_gene_spots
    spots = calculate_spots_per_cell(spots)

    #Calculate gene centroid within each cell
    gene_centroids = spots.groupby(['cell_id','target_molecule_name','num_cell_spots','num_gene_spots'])[['global_x','global_y']].mean().reset_index()
    gene_centroids = gene_centroids.rename(columns={'global_x':'gene_centroidX','global_y':'gene_centroidY'})

    #Calculate distance of gene_centroid to cell_centroid
    gene_centroids = add_cell_centroids_to_spots(cells, gene_centroids)
    gene_centroids['raw_metric'] = np.sqrt(
        ((gene_centroids['cell_centroidX']-gene_centroids['gene_centroidX'])**2) +
        ((gene_centroids['cell_centroidY']-gene_centroids['gene_centroidY'])**2)
    )
    gene_centroids['metric_name'] = 'polarity'

    #Get the output into the standard metric table
    cols = [
       'cell_id', 'target_molecule_name',
       'num_cell_spots', 'num_gene_spots',
       'metric_name', 'raw_metric',
    ]
    metric_df = gene_centroids[cols]

    return metric_df
    
@metric_summaries
def spatial_metric_centrality(spots,cells):
    """
    Spatial metric: centrality

    For each gene in each cell, calculates the mean distance between the RNA spots of that gene from the cell centroid
    Intuition is that genes that have low mean distance from cell centroid are centrally located

    Arguments:
        - spots: dataframe object where each row is a MERFISH spot with the following columns
            'cell_id': the id of the cell that contains this spot as a string
            'target_molecule_name': gene name
            'global_x': x-coordinate of the spot
            'global_y': y-coordinate of the spot
            

        - cells: dataframe object where each row is a MERFISH cell with the following columns
            'cell_id': the ide of the cell as a string
            'boundaryX': string of x-coords for the cell boundary like
                "2097.1362953431903, 2097.1362953431903, 2097.1362953431903, ...."
            'boundaryY': string of y-coords for the cell boundary like
                "2097.1362953431903, 2097.1362953431903, 2097.1362953431903, ...."

    Outputs:
        - metric_df: dataframe object where each row is a unique gene/cell combination. columns:
            'cell_id': same as from inputs
            'target_molecule_name': same as from inputs
            'num_cell_spots': number of total spots in this cell
            'num_gene_spots': number of spots of the specific gene in this cell
            'metric_name': always 'centrality' for this function
            'raw_metric': the raw metric value of gene-centroid dist to cell centroid
            any/all processed metrics from the metric_summaries wrapper
    """
    #For each spot, calculate distance to cell centroid
    spots = add_cell_centroids_to_spots(cells, spots)
    spots['raw_metric'] = np.sqrt(
        ((spots['cell_centroidX']-spots['global_x'])**2) +
        ((spots['cell_centroidY']-spots['global_y'])**2)
    )
   
    #Calculate mean centroid_dist for each gene in each cell 
    centrality_df = spots.groupby(['cell_id','target_molecule_name','num_cell_spots','num_gene_spots'])['raw_metric'].mean().reset_index()
    centrality_df['metric_name'] = 'centrality'

    #Get the output into the standard metric table
    cols = [
       'cell_id', 'target_molecule_name',
       'num_cell_spots', 'num_gene_spots',
       'metric_name', 'raw_metric',
    ]
    metric_df = centrality_df[cols]

    return metric_df


@metric_summaries
def spatial_metric_periphery(spots,cells):
    pass 

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
 
def calculate_spots_per_cell(spots):
    if 'num_cell_spots' not in spots.columns:
        spots['num_cell_spots'] = spots.groupby('cell_id')['cell_id'].transform('count')
        
    if 'num_gene_spots' not in spots.columns:
        spots['num_gene_spots'] = spots.groupby(['cell_id','target_molecule_name'])['cell_id'].transform('count')

    return spots


metrics = {
    'polarity': spatial_metric_polarity,
    'centrality': spatial_metric_centrality,
    'periphery': spatial_metric_periphery,
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run spatial metrics on MERFISH spots and cells')
    parser.add_argument('--cells',dest='cells_path', help='path to cell info', required=True)
    parser.add_argument('--spots',dest='spots_path', help='path to spot info', required=True)
    parser.add_argument('--out',dest='out_path', help='path to save metric output', required=True)
    parser.add_argument('--metric',dest='metric', help='which metric to run', choices=metrics.keys(), required=True)

    args = parser.parse_args()

    metric_func = metrics[args.metric]

    spots = pd.read_csv(args.spots_path)
    cells = pd.read_csv(args.cells_path)

    out_df = metric_func(spots, cells)

    out_df.to_csv(args.out_path, index=False)




