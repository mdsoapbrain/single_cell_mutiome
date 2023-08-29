#!/bin/bash

#
#PBS  -N annot_files_olg_nn_diffdiff_peaks
#
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=80:mem=400gb

## Create conda env (DO THIS STEP BEFORE ATTEMPTING TO RUN LDSC) 
#conda env create --file environment.yml

module load anaconda3/personal
source activate ldsc
# Input parameters
BED_DIR="/rds/general/user/hs2722/home/lifted_bed_files_hg/"
OUTPUT_DIR="/rds/general/user/hs2722/home/LDSC_results/"
BIM_DIR="/rds/general/user/hs2722/home/REF"

# Other parameters
RESULTS_DIR="/rds/general/user/hs2722/home/tmp_3/"
sumstats_dir="/rds/general/user/hs2722/home/LDSC_inputs_3/"
weights_dir="/rds/general/user/hs2722/home/REF/weights_hm3_no_hla/"

# Get list of all bed files
BED_FILES=( ${BED_DIR}*.bed )

# Generate annotation files and calculate LD Scores for each chromosome
for BED_FILE in ${BED_FILES[@]}; do
  samp=$(basename $BED_FILE .bed)
  
  for CHR in {1..22}; do
    BIMFILE=${BIM_DIR}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}.bim

    python /rds/general/user/hs2722/home/ldsc/make_annot.py --bed-file ${BED_FILE} --bimfile ${BIMFILE} --annot-file ${OUTPUT_DIR}/${samp}.${CHR}.annot.gz 

    python /rds/general/user/hs2722/home/ldsc/ldsc.py --l2 --bfile ${BIM_DIR}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR} --ld-wind-cm 1 --annot ${OUTPUT_DIR}/${samp}.${CHR}.annot.gz --thin-annot --out ${OUTPUT_DIR}/${samp}.${CHR} \
    --print-snps /rds/general/user/hs2722/home/REF/1000G_EUR_Phase3_baseline/print_snps.txt &> ${OUTPUT_DIR}${samp}.${CHR}.ldsclog
  done

  # Run LDSC analysis for each sumstats file
  sumstats_files=( "${sumstats_dir}"*.sumstats )
  for sumstats_file in "${sumstats_files[@]}"; do
    # Extract the filename without the directory path and extension
    filename=$(basename "$sumstats_file")
    filename_without_ext="${filename%.*}"

    # Output file path
    output_file="${RESULTS_DIR}${samp}_${filename_without_ext}_h2"

    # Run LDSC analysis
    python /rds/general/user/hs2722/home/ldsc/ldsc.py \
    --h2 "$sumstats_file" \
    --ref-ld-chr ${OUTPUT_DIR}${samp}.,${BIM_DIR}/1000G_EUR_Phase3_baseline/baseline. \
    --w-ld-chr "${weights_dir}weights." \
    --out "$output_file" \
    --overlap-annot \
    --frqfile-chr /rds/general/user/hs2722/home/REF/1000G_Phase3_frq/1000G.EUR.QC. \
    --print-coefficients

    # Print a newline for clarity
    echo
  done
done
