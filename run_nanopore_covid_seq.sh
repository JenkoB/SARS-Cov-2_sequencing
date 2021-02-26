#!/usr/bin/env bash

#===============================================================================
#
#          FILE:  run_nanopore_covid_seq.sh
#		https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
#
#       VERSION:  1.0.0_SuperMicro
#       CREATED:  09/02/2021
#===============================================================================

print_usage() {
  printf "Usage: 
		conda activate artic-ncov2019 
		run_nanopore_covid_seq.sh [Options] Path_to_fast5_dir\n
        Options: 
        -n [N]  number of threads used in bcbio [Default: 20]
        -o [SEQ_FOLDER] name/path of/to output dir - run name [Default: MinION_seq]"	
}

# Default values
THREADS=200
SEQ_FOLDER="MinION_seq"

# Options
while getopts 'n:o:h' OPTION; do
  case "$OPTION" in
    n) THREADS="$OPTARG" ;;
    o) SEQ_FOLDER="$OPTARG" ;;
    h) print_usage
    exit 1 ;;
  esac
done

shift "$(($OPTIND -1))"

if [ -z "$*" ]; then print_usage && exit 1; fi
FAST5=("$@")

PRIMER_SCHEMES_v3="/home/kisld/bioinformatic_tools/artic-ncov2019/primer_schemes"

mkdir -p "${SEQ_FOLDER}"/"demultiplex"
mkdir -p "${SEQ_FOLDER}"/"fastq_files"

guppy_basecaller \
    -i "${FAST5}" \
    -s "${SEQ_FOLDER}"/"fastq_files" \
    -c dna_r9.4.1_450bps_hac.cfg \
    -x auto

 
echo "${THREADS}"

guppy_barcoder \
    --trim_barcodes \
    -i "${SEQ_FOLDER}"/"fastq_files" \
    -s "${SEQ_FOLDER}"/demultiplex \
    --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg" \
    -x auto \
    -t "${THREADS}" \
    --num_barcoding_buffers 20

# source activate /home/kisld/anaconda3/envs/artic-ncov2019

for sample in "${SEQ_FOLDER}"/demultiplex/*; do
    
    echo "${sample}"
    [ -d "${sample}" ] || continue 
    FILENAME=`basename "${sample}" `
    echo -e "\e[33m"${FILENAME}"\033[0m"

    mkdir -p "${SEQ_FOLDER}"/"results"/"${FILENAME}"
    echo "${FILENAME}"
    
    artic guppyplex \
	    --min-length 400 \
	    --max-length 1000 \
	    --directory "${SEQ_FOLDER}"/demultiplex/"${FILENAME}" \
	    --output "${SEQ_FOLDER}"/"results"/"${FILENAME}"/"${FILENAME}".fastq

    
    artic minion \
	    --normalise 200 \
	    --threads "${THREADS}" \
            --scheme-directory "${PRIMER_SCHEMES_v3}" \
            --read-file "${SEQ_FOLDER}"/"results"/"${FILENAME}"/"${FILENAME}".fastq \
            --fast5-directory "${FAST5}" \
            --sequencing-summary "${SEQ_FOLDER}"/"fastq_files"/sequencing_summary.txt \
            nCoV-2019/V3 \
            "${SEQ_FOLDER}"/"results"/"${FILENAME}"/"${FILENAME}"

        
done

cat "${SEQ_FOLDER}"/"results"/*/*.consensus.fasta > "${SEQ_FOLDER}"/"results"/"${SEQ_FOLDER}".consensus_all.fasta
