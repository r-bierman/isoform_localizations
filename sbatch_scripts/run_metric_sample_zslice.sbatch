#!/bin/bash
#############################
# File Name : run_metric_calc.sbatch
#
# Purpose : [???]
#
# Creation Date : 11-01-2021
#
# Last Modified : Thu 14 Jan 2021 03:15:07 PM PST
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
#SBATCH --time=4:00:00                     # how much time to run
#SBATCH --mem=15000                         # how much mem in MBs
#SBATCH --job-name=metric


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
#SBATCH --array=1-77                       # to submit an array of jobs
##SBATCH --reservation r_covid19            # Uncomment this line and next line to get on covid19 queue
##SBATCH --account a_covid19                # Uncomment this line and previous line to get on covid19 queue
echo $SLURM_ARRAY_TASK_ID                  # this variable takes on the int of the current array

#Pull sample and zslice from a list
sample=`head -n $SLURM_ARRAY_TASK_ID ../processed_data/sample_z_slice_flags.txt | tail -n 1 | cut -f1 -d_`
zslice=`head -n $SLURM_ARRAY_TASK_ID ../processed_data/sample_z_slice_flags.txt | tail -n 1 | cut -f2 -d_`

out="20210114_polarity_batched_outs/20210114_polarity_${sample}_${zslice}.csv"

echo $sample
echo $zslice
echo $out

python3 ../scripts/spatial_utils.py \
    --cells ../processed_data/segmented_cell_shapes.csv \
    --spots ../processed_data/all_passing_rna_spots.csv \
    --out $out \
    --metric polarity \
    --sample $sample \
    --zslice $zslice \


