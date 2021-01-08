import pandas as pd
import numpy as np
import time
import sys

from matplotlib.path import Path

def assign_points_to_cells(cells_path, spots_path, out_path):
    #Read in and clean cells df
    cells = pd.read_csv(cells_path)
    cells = cells.rename(columns={'Unnamed: 0':'cell_id'})

    cells['x_points'] = cells['boundaryX'].str.split(', ').apply(lambda xs: [float(x) for x in xs]).values
    cells['y_points'] = cells['boundaryY'].str.split(', ').apply(lambda ys: [float(y) for y in ys]).values
    cells['boundary'] = cells.apply(lambda r: np.vstack((r['x_points'],r['y_points'])).T, axis=1)
    cells['path'] = cells['boundary'].apply(Path)
    
    points = pd.read_csv(spots_path).drop(columns='Unnamed: 0')
    
    #Loop through each cell and 
    ps = points[['global_x','global_y']].values

    points['num_assigned_cells'] = 0
    
    for i,r in cells.iterrows():
        path = r['path']


        containing = path.contains_points(ps)
        points.loc[containing, 'cell_id'] = r['cell_id']
        points.loc[containing, 'num_assigned_cells'] += 1 #just used to see if there is overwriting occuring

    #Merge the cell info into the points on the cell_id
    points = pd.merge(
        left = points,
        right = cells[['cell_id','boundaryX','boundaryY','slice_id']],
        on = 'cell_id',
        how = 'left',
    )
    
    points.to_csv(out_path, index=False)
   

def assign_points_to_cells_faster(cells_path, spots_path, out_path):
    #Read in and clean cells df
    cells = pd.read_csv(cells_path)
    cells = cells.rename(columns={'Unnamed: 0':'cell_id'})

    cells['x_points'] = cells['boundaryX'].str.split(', ').apply(lambda xs: [float(x) for x in xs]).values
    cells['y_points'] = cells['boundaryY'].str.split(', ').apply(lambda ys: [float(y) for y in ys]).values
    cells['boundary'] = cells.apply(lambda r: np.vstack((r['x_points'],r['y_points'])).T, axis=1)
    cells['path'] = cells['boundary'].apply(Path)
    
    points = pd.read_csv(spots_path)
    print(points.shape)
    sys.stdout.flush()
    
    #Loop through each cell and 
    ps = points[['global_x','global_y']].values
    
    points['num_assigned_cells'] = 0

    start = time.time()
    print(cells.shape)
    for i,r in cells.iterrows():
        if i%500 == 0:
            print(i)
            print(time.time()-start)
            sys.stdout.flush()
            start = time.time()    

        path = r['path']


        containing = path.contains_points(ps)
        points.loc[containing, 'cell_id'] = r['cell_id']
        points.loc[containing, 'num_assigned_cells'] += 1 #just used to see if there is overwriting occuring
        points.loc[containing, 'boundaryX'] = r['boundaryX']
        points.loc[containing, 'boundaryY'] = r['boundaryY']
        points.loc[containing, 'slice_id'] = r['slice_id']

    points.to_csv(out_path, index=False)


     
