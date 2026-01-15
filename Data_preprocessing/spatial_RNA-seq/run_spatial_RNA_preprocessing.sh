#!/bin/bash

#SBATCH --job-name=A2TR4
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --error=%j.err
#SBATCH --output=%j.out

# 1
core=5
mem=90
MAP=/mnt/hpc/users/data/3download/UCSC_Maf/STAR_index
ANN=/mnt/hpc/users/data/3download/UCSC_Maf/T2T-MFA8v1.0_genomic.chrname.transcript_as_exon.gtf
CONT=/mnt/hpc/users/data/3download/Macaca_rDNA/STAR_index
R_path=/mnt/hpc/users/anaconda3/envs/SRT_env/bin/Rscript

# 2 Variable parameter
batch=exp1_0123
slice=A2TR4
directory=/mnt/hpc/users/data/2results/T2T/A2TR4
# 3 Variable parameter
ID_File=/mnt/hpc/users/data/Spatial_barcodes/combined_96_barcode.txt
#ID_File=/mnt/hpc/users/data/Spatial_barcodes/combined_196_barcode.txt
RNA_R1=/mnt/hpc/users/data/1rawdata/A2TR4/A2TR4_R1.fq.gz
RNA_R2=/mnt/hpc/users/data/1rawdata/A2TR4/A2TR4_R2.fq.gz

# Fixed combination
name=RNA_${batch}_${slice}
out_path=${directory}/scRNA_${slice}

# 4
/mnt/hpc/users/data/2results/T2T/A2TR4/extract_barcode.sh $name $RNA_R1 $RNA_R2 $out_path/fastq
# 5
$R_path /mnt/hpc/users/data/2results/T2T/A2TR4/check_barcode.R -b $ID_File -p $out_path/fastq -n ${name}_2.extract.barcode
rm $out_path/fastq/${name}_2.extract.barcode

# 6
#conda activate spatail_env
/mnt/hpc/users/data/2results/T2T/A2TR4/run_st_pipeline.soft_clipping.sh $core $mem $name $out_path \
${out_path}/fastq/${name}_2.extract.fq.gz ${out_path}/fastq/${name}_1.extract.fq.gz $MAP $ANN $CONT $ID_File
# 7
rm ${out_path}/1_stpipeline/tmp/R2_quality_trimmed.bam ${out_path}/1_stpipeline/tmp/contaminated_clean.bam ${out_path}/1_stpipeline/tmp/demultiplexed_matched.bam

