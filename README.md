
Identifying subcellular spatial patterning in MERFISH data from: 
"Molecular, spatial and projection diversity of neurons in primary motor cortex revealed by in situ single-cell transcriptomics"

https://www.biorxiv.org/content/10.1101/2020.06.04.105700v1

Downloaded data on 12/21/2020 from ftp://download.brainimagelibrary.org:8811/02/26/02265ddb0dae51de/

pre-processed data to:
* Assign RNA spots into cells by combining info from segmented* and spot* files
* Combining all cell boundaries into a single file



Directory layout:

notebooks:
* variety of jupyter notebooks

processed_data
* counts.h5ad
    Normalized gene counts per cell
* segmented_cell_shapes.csv
    All cell boundaries from all samples
* all_passing_rna_spots.csv
    All rna spot cell assignments for non-null spots in "kept cells" in the counts.h5ad table

scripts:
* plot_utils.py
    utilities for plotting
* utils.py
    general utilities

unprocessed_downloads.tar.gz:
* dl_scripts
    scripts used to download the MERFISH data
* segmented_cells
    cell boundary info for separate samples
    this info is consolidated in sgemented_cell_shapes.csv
* spot_assigning
    code and output to assign spots to cells
* spot_assignments
    per sample spot assignment into cells
    includes unassigned spots and spot assignments to cells not in the counts.h5ad file
