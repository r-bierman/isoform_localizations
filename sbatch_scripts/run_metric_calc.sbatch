#!/bin/bash
#############################
# File Name : run_metric_calc.sbatch
#
# Purpose : [???]
#
# Creation Date : 11-01-2021
#
# Last Modified : Fri 15 Jan 2021 11:59:41 AM PST
#
# Created By : Rob Bierman
#
##############################

#########################################################
#                                                       #
#      Parameters that should be set every time         #
#                                                       #
#########################################################
#Single commented out SBATCH are active, more than that is inactive
#SBATCH -p owners # which node to run on
#SBATCH --time=16:00:00                     # how much time to run
#SBATCH --mem=45000                         # how much mem in MBs
#SBATCH --job-name=periph


#########################################################
#                                                       #
#               More optional parameters                #
#                                                       #
#########################################################
##SBATCH --output=slurm.out                 # slurm.out file
##SBATCH --error=slurm.err                  # slurm.err file
##SBATCH --mail-type=END,FAIL               # notifications for job done & fail
##SBATCH --mail-user=rbierman@stanford.edu  # email to send it to
##SBATCH --dependency=afterany:???      # wait to submit until this job num finishes
##SBATCH --array=1-35                       # to submit an array of jobs
##SBATCH --reservation r_covid19            # Uncomment this line and next line to get on covid19 queue
##SBATCH --account a_covid19                # Uncomment this line and previous line to get on covid19 queue
#echo $SLURM_ARRAY_TASK_ID                  # this variable takes on the int of the current array

python3 ../scripts/spatial_utils.py \
    --cells ../processed_data/segmented_cell_shapes.csv \
    --spots ../processed_data/all_passing_rna_spots.csv \
    --out 20210115_periphery \
    --metric periphery \



