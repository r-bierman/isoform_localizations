#!/bin/bash
#############################
# File Name : run_perms.sbatch
#
# Purpose : [???]
#
# Creation Date : 27-01-2021
#
# Last Modified : Thu 28 Jan 2021 07:26:25 AM PST
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
#SBATCH -p quake,owners                          # which node to run on
#SBATCH --time=1:00:00                     # how much time to run
#SBATCH --mem=10000                         # how much mem in MBs
#SBATCH --job-name=mat_perms


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
#SBATCH --array=1-252%50                       # to submit an array of jobs
##SBATCH --reservation r_covid19            # Uncomment this line and next line to get on covid19 queue
##SBATCH --account a_covid19                # Uncomment this line and previous line to get on covid19 queue
echo $SLURM_ARRAY_TASK_ID                  # this variable takes on the int of the current array

GENE=`head -n $SLURM_ARRAY_TASK_ID gene_names.txt | tail -n 1`
echo $GENE

python3 matrix_perms.py $GENE


