#!/bin/bash
#############################
# File Name : meta_downloads.sbatch
#
# Purpose : Download metadata from yao et al.
#
# Creation Date : 03-02-2021
#
# Last Modified : Wed 03 Feb 2021 02:13:01 PM PST
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
#SBATCH -p quake,owners,horence                          # which node to run on
#SBATCH --time=8:00:00                     # how much time to run
#SBATCH --mem=5000                         # how much mem in MBs
#SBATCH --job-name=meta_dls


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



wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/QC.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/barcode.tsv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/cluster.annotation.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/cluster.membership.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/features.tsv.gz
wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/matrix.mtx.gz
wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/sample_metadata.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/umi_counts.h5
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/QC.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/barcode.tsv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/cluster.annotation.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/cluster.membership.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/features.tsv.gz
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/matrix.mtx.gz
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/sample_metadata.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v3_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/umi_counts.h5
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/QC.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/barcode.tsv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/cluster.annotation.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/cluster.membership.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/features.tsv.gz
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/matrix.mtx.gz
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/sample_metadata.csv
wget -P Zeng_Transcriptomic_Allen_10X_cells_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/umi_counts.h5
wget -P Zeng_Transcriptomic_Allen_SMARTer_cells_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/QC.csv
wget -P Zeng_Transcriptomic_Allen_SMARTer_cells_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/cluster.annotation.csv
wget -P Zeng_Transcriptomic_Allen_SMARTer_cells_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/cluster.membership.csv
wget -P Zeng_Transcriptomic_Allen_SMARTer_cells_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/exon.counts.csv.gz
wget -P Zeng_Transcriptomic_Allen_SMARTer_cells_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/intron.counts.csv.gz
wget -P Zeng_Transcriptomic_Allen_SMARTer_cells_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/sample_metadata.csv.gz
wget -P Zeng_Transcriptomic_Allen_SMARTer_cells_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/tsne.df.csv
wget -P Zeng_Transcriptomic_Allen_SMARTer_nuclei_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/QC.csv
wget -P Zeng_Transcriptomic_Allen_SMARTer_nuclei_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/cluster.annotation.csv
wget -P Zeng_Transcriptomic_Allen_SMARTer_nuclei_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/cluster.membership.csv
wget -P Zeng_Transcriptomic_Allen_SMARTer_nuclei_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/exon.counts.csv.gz
wget -P Zeng_Transcriptomic_Allen_SMARTer_nuclei_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/intron.counts.csv.gz
wget -P Zeng_Transcriptomic_Allen_SMARTer_nuclei_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/sample_metadata.csv.gz
wget -P Zeng_Transcriptomic_Allen_SMARTer_nuclei_MOp/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/tsne.df.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/QC.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/barcode.tsv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/cluster.annotation.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/cluster.membership.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/features.tsv.gz
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/matrix.mtx.gz
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/sample_metadata.csv
wget -P Zeng_Transcriptomic_Allen_10X_nuclei_v2_AIBS/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/umi_counts.h5
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/QC.csv
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/barcodes.csv
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/cluster.annotation.csv
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/cluster.membership.csv
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/features.csv
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/matrix.mtx.gz
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/sample_metadata.csv
wget -P BICCN_MOp_MiniAtlas_10X_nuclei_v3_Broad/ http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/umi_counts.h5

