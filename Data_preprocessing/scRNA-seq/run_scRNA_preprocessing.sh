#!/bin/bash

#SBATCH --job-name=A2T
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --error=%j.err
#SBATCH --output=%j.out

cellranger count --id A2T_run --transcriptome /mnt/hpc/users/Ref_Genome/Crab_eating_macaque/cellranger/MacaqueT2T --fastqs /mnt/hpc/users/01_raw_data/scRNA_testis/A2 --sample A2
