import pandas as pd
import numpy as np
import argparse
import os

import functools
import hashlib
import shutil

import shapely

pd.options.mode.chained_assignment = None  # default='warn'


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
        
        cell_metric_medians = df.groupby('cell_id')['raw_metric'].transform('median')
        cell_metric_vars = df.groupby('cell_id')['raw_metric'].transform('std')**2
        df['cell_zscore'] = (df['raw_metric']-cell_metric_medians)/cell_metric_vars

        df['gene_median_cell_zscore'] = df.groupby('target_molecule_name')['cell_zscore'].transform('median')
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
    gene_centroids = spots.groupby(['cell_id','target_molecule_name','num_cell_spots','num_gene_spots'])[['global_x','global_y']].median().reset_index()
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

    For each gene in each cell, calculates the median distance between the RNA spots of that gene from the cell centroid
    Intuition is that genes that have low median distance from cell centroid are centrally located

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
    #Calculate num_cell_spots and num_gene_spots
    spots = calculate_spots_per_cell(spots)

    #For each spot, calculate distance to cell centroid
    spots = add_cell_centroids_to_spots(cells, spots)
    spots['raw_metric'] = np.sqrt(
        ((spots['cell_centroidX']-spots['global_x'])**2) +
        ((spots['cell_centroidY']-spots['global_y'])**2)
    )

    #Calculate median centroid_dist for each gene in each cell 
    centrality_df = spots.groupby(['cell_id','target_molecule_name','num_cell_spots','num_gene_spots'])['raw_metric'].median().reset_index()
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
    """
    Spatial metric: periphery

    Calculates the median distance between the RNA spot and the nearest cell boundary for each gene
    Intuition is that genes closer to the periphery will have shorter median distances

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
            'metric_name': always 'periphery' for this function
            'raw_metric': the raw metric value of gene-centroid dist to cell centroid
            any/all processed metrics from the metric_summaries wrapper
    """
    #Remove cells that have too few boundary points
    #Having 3 commas means 4 points which is min for shapely
    cells = cells[cells['boundaryX'].str.count(',').gt(3)]

    #Limit the spots to just the filtered cells above
    spots = spots[spots['cell_id'].isin(cells['cell_id'])]

    #Calculate num_cell_spots and num_gene_spots
    spots = calculate_spots_per_cell(spots)

    #Create a shapely.poly object for each cell
    cells['poly'] = cells.apply(lambda r: create_shapely_boundary(r['boundaryX'],r['boundaryY']), axis=1)

    #Add the poly object to each row in spots
    cell_id_to_poly = cells.set_index('cell_id')['poly']
    spots['poly'] = spots['cell_id'].map(cell_id_to_poly)

    #Calculate min boundary dists for each spot
    def calculate_min_boundary_dist(poly,px,py):
        pt = shapely.wkt.loads('POINT({} {})'.format(px,py))
        return poly.boundary.distance(pt)

    spots['raw_metric'] = (
        spots.apply(
            lambda r: calculate_min_boundary_dist(
                r['poly'], 
                r['global_x'],
                r['global_y']),
            axis=1)
    )

    #Calculate median centroid_dist for each gene in each cell 
    periphery_df = spots.groupby(['cell_id','target_molecule_name','num_cell_spots','num_gene_spots'])['raw_metric'].median().reset_index()
    periphery_df['metric_name'] = 'periphery'

    #Get the output into the standard metric table
    cols = [
       'cell_id', 'target_molecule_name',
       'num_cell_spots', 'num_gene_spots',
       'metric_name', 'raw_metric',
    ]
    metric_df = periphery_df[cols]


    return metric_df


def create_shapely_boundary(boundaryX,boundaryY):
    """
    Create a shapely.poly object from cell str boundaires
    """
    xs = [float(x) for x in boundaryX.split(', ')]
    ys = [float(y) for y in boundaryY.split(', ')]
    xs = xs+[xs[0]] #last point must be the same as the first point
    ys = ys+[ys[0]] #last point must be the same as the first point
    str_bounds = ', '.join([str(x)+' '+str(y) for x,y in zip(xs,ys)])

    poly = shapely.wkt.loads('POLYGON(({}))'.format(str_bounds))    
    return poly 


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


def self_hash():
    """
    For versioning
    Read the entire text of this script and hash it
    Copy this script to a dir of frozen hashed scripts if not already there
    Return the hash

    Any change to this file, even comments, changes the hash
    """
    frozen_dir = 'frozen_hashed_scripts'

    self_path = os.path.realpath(__file__)
    with open(self_path,'r') as f_in:
        all_text = f_in.read()

    h = hashlib.sha256(all_text.encode()).hexdigest()

    if not os.path.exists(frozen_dir):
        os.mkdir(frozen_dir)

    frozen_path = os.path.join(frozen_dir,h+'.py')
    if not os.path.exists(frozen_path):
        shutil.copy(self_path, frozen_path)

    return h


metrics = {
    'polarity': spatial_metric_polarity,
    'centrality': spatial_metric_centrality,
    'periphery': spatial_metric_periphery,
}

if __name__ == '__main__':

    #Gets a hash for this script
    code_version = self_hash()

    parser = argparse.ArgumentParser(description='Run spatial metrics on MERFISH spots and cells')
    parser.add_argument('--cells',dest='cells_path', help='path to cell info', required=True)
    parser.add_argument('--spots',dest='spots_path', help='path to spot info', required=True)
    parser.add_argument('--out',dest='out_path', help='path to save metric output', required=True)
    parser.add_argument('--metric',dest='metric', help='which metric to run', choices=metrics.keys(), required=True)
    parser.add_argument('--sample',dest='sample', help='which sample to restrict to, if any', required=False)
    parser.add_argument('--zslice',dest='zslice', help='which zslice to restrict to, if any', required=False)

    args = parser.parse_args()
    print(args)

    metric_func = metrics[args.metric]

    spots = pd.read_csv(args.spots_path)
    cells = pd.read_csv(args.cells_path)

    if args.sample:
        spots = spots[spots['sample'].eq(args.sample)]
    if args.zslice:
        spots = spots[spots['global_z'].eq(float(args.zslice))]

    #Limit to cells that have at least one spot
    cells = cells[cells['cell_id'].isin(spots['cell_id'].unique())]

    print(spots.shape)
    print(cells.shape)

    out_df = metric_func(spots, cells)
    out_df['spatial_utils_code_version'] = code_version

    out_df.to_csv(args.out_path, index=False)




