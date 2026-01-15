#!/bin/bash

#SBATCH --job-name=A2TA4  
#SBATCH --partition=cpuPartition   
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=20
#SBATCH --error=%j.err
#SBATCH --output=%j.out

umi_tools extract --extract-method=regex --bc-pattern2="^(?P<umi_1>.)(?P<discard_1>ATCGGCGTACGACT){s<=1}.{8}(?P<discard_2>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=1}.{8}(?P<discard_3>CCCATGATCGTCCGATGCAGTCGTGCCATGAGATGTGTATAAGAGACAG){s<=2,i<=1,d<=1}.*$" -I /mnt/hpc/users/01_raw_data/ATAC/A2TA4/A2TA4_R1.fq.gz -S A2TA4_1.extract.fq.gz --read2-in=/mnt/hpc/users/01_raw_data/ATAC/A2TA4/A2TA4_R2.fq.gz --read2-out=A2TA4_2.extract.fq.gz -L ./extract.log

zcat A2TA4_2.extract.fq.gz | awk 'NR%4==2{print substr($1,1,16)}' > A2TA4_2.extract.barcode
mkdir fastq_test
mv A2TA4_2.extract.barcode extract_barcode

zcat A2TA4_2.extract.fq.gz | seqkit subseq -r 1:16 | gzip -c - > A2TA4_2.barcode.extract.fq.gz
zcat A2TA4_2.extract.fq.gz | seqkit subseq -r 17:-1 | gzip -c - > A2TA4_2.read2.extract.fq.gz

mv A2TA4_1.extract.fq.gz A2TA4_S1_L001_R1_001.fastq.gz
mv A2TA4_2.barcode.extract.fq.gz A2TA4_S1_L001_R2_001.fastq.gz
mv A2TA4_2.read2.extract.fq.gz A2TA4_S1_L001_R3_001.fastq.gz
mkdir fastq
mv A2TA4_S1_* fastq

cellranger-atac count --id=A2TA4 --fastqs=./fastq --sample=A2TA4 \
--reference=/mnt/hpc/users/data/3download/Macaca_T2T \
--jobmode=local --localcores=30 --localmem=80

