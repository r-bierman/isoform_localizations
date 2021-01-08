#!/bin/bash
#############################
# File Name : run_all_centroids.sbatch
#
# Purpose : [???]
#
# Creation Date : 06-01-2021
#
# Last Modified : Wed 06 Jan 2021 10:18:36 AM PST
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
#SBATCH -p horence                          # which node to run on
#SBATCH --time=48:00:00                     # how much time to run
#SBATCH --mem=60000                         # how much mem in MBs
#SBATCH --job-name=centroids


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

python3 all_centroids.py

