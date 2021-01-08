import numpy as np
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import seaborn as sns

def plot_spot_cells(spots, cells, spot_colors='k', spot_size=0.05):
    #Limit to just spots that have matching cells in the cell list
    spots = spots[spots['cell_id'].isin(cells['cell_id'])]

    for cell_id,spot_g in spots.groupby('cell_id'):
        boundaryX = cells.loc[cells['cell_id'].eq(cell_id),'boundaryX'].values[0]
        boundaryY = cells.loc[cells['cell_id'].eq(cell_id),'boundaryY'].values[0]

        fig,ax = _plot_cell(boundaryX,boundaryY)

        ax = _add_spots(ax, spot_g, color=spot_colors, size=spot_size)

        plt.title(cell_id)
        plt.show()
        plt.close()

        

def _add_spots(ax, spots, color='r', size=0.2):
    """
    Draws spots on top of cell polygons
    
    Input is a dataframe with global_x and global_y columns
    """
    if color == 'by_gene':
        genes = spots['target_molecule_name'].unique()
        color_palette = sns.color_palette("hls", genes.size)
        gene_colors = {g:c for g,c in zip(genes,color_palette)}
        colors = spots['target_molecule_name'].map(gene_colors).values
    elif type(color) == dict:
        genes = spots['target_molecule_name'].unique()
        gene_colors = {g: color[g] if g in color else 'k' for g in genes}
        colors = spots['target_molecule_name'].map(gene_colors).values
    elif type(color) == str:
        colors = [color]*spots.shape[0]
    else:
        print('color {} not recognized'.format(color))

    circs = []

    for i,(ind,s) in enumerate(spots.iterrows()):
        c = Circle(
            (s['global_x'],s['global_y']),
            size,
            color=colors[i],
        )
        circs.append(c)
    
    p = PatchCollection(circs, match_original=True)
    ax.add_collection(p)
    return ax

def _plot_cell(boundaryX,boundaryY):
    xs = [float(x) for x in boundaryX.split(', ')]
    ys = [float(y) for y in boundaryY.split(', ')]

    points = np.array([xs,ys]).T
    xmin,ymin = points.min(axis=0)
    xmax,ymax = points.max(axis=0)

    polygon = Polygon(points, fill=None)
    p = PatchCollection([polygon], match_original=True)

    fig, ax = plt.subplots(figsize=(8,8))
    ax.add_collection(p)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.axis('equal')

    return fig,ax

   


def plot_cells(cells):
    """
    cells is a dataframe from a segmented_cells_xxx.csv file
    """
    patches = []

    xmin = None
    xmax = None
    ymin = None
    ymax = None
    
    cells = cells.drop_duplicates(subset=['boundaryX','boundaryY'])
    cells = cells.dropna(subset=['boundaryX','boundaryY'])

    for i in range(cells.shape[0]):
        xs = [float(x) for x in cells.iloc[i]['boundaryX'].split(', ')]
        ys = [float(y) for y in cells.iloc[i]['boundaryY'].split(', ')]

        points = np.array([xs,ys]).T

        polygon = Polygon(points)
        patches.append(polygon)

        sub_xmin,sub_ymin = points.min(axis=0)
        sub_xmax,sub_ymax = points.max(axis=0)

        if xmin is None or xmin > sub_xmin:
            xmin = sub_xmin
        if xmax is None or xmax < sub_xmax:
            xmax = sub_xmax
        if ymin is None or ymin > sub_ymin:
            ymin = sub_ymin
        if ymax is None or ymax < sub_ymax:
            ymax = sub_ymax



    p = PatchCollection(patches)

    fig, ax = plt.subplots()
    ax.add_collection(p)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.axis('equal')

    return fig,ax


