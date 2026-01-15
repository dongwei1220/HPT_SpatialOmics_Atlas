#!/bin/bash
if [ $# -lt 1 ]; then
    echo "error.. need args"
    exit 1
fi

ppn=$1 # 12
mem=$2 # 60 
name=$3 # A2TR4
path=$4 # /mnt/hpc/users/data/2results/T2T/A2TR4/scRNA_A2TR4
FW=$5 # /mnt/hpc/users/data/1rawdata/A2TR4/A2TR4_R1.fq.gz
RV=$6 # /mnt/hpc/users/data/1rawdata/A2TR4/A2TR4_R2.fq.gz
MAP=$7 # /mnt/hpc/users/data/3download/UCSC_Maf/STAR_index
ANN=$8 # /mnt/hpc/users/data/3download/UCSC_Maf/T2T-MFA8v1.0_genomic.chrname.transcript_as_exon.gtf
CONT=$9 # /mnt/hpc/users/data/3download/Macaca_rDNA/STAR_index
ID_File=${10} # /mnt/hpc/users/data/Spatial_barcodes/combined_96_barcode.txt
#ID_File=${10} # /mnt/hpc/users/data/Spatial_barcodes/combined_192_barcode.txt
OUTPUT=1_stpipeline
TMP=1_stpipeline/tmp

cd $path
mkdir -p $TMP

st_pipeline_run.py \
  --output-folder $OUTPUT \
  --temp-folder $TMP \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --ids $ID_File \
  --ref-map $MAP \
  --ref-annotation $ANN \
  --expName $name \
  --htseq-no-ambiguous \
  --verbose \
  --threads $ppn \
  --log-file $OUTPUT/${name}_log.txt \
  --star-two-pass-mode \
  --no-clean-up \
  --contaminant-index $CONT \
  --min-length-qual-trimming 30 \
  --star-sort-mem-limit ${mem}000000000  \
  $FW $RV

