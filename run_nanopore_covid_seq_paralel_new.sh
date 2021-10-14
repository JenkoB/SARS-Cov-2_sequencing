#!/usr/bin/env bash

#===============================================================================
#
#          FILE:  run_nanopore_covid_seq.sh
#	ORIGINAL INSTRUCTIONS:  https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
#
#       VERSION:  1.0.0
#       CREATED:  26/02/2021
#===============================================================================

print_usage() {
  printf "Usage: 
		conda activate artic-ncov2019 
		run_nanopore_covid_seq.sh [Options] Path_to_fast5_dir\n
        Options: 
        -n [N]  number of threads used [Default: 200]
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

analysis () {

    sample="$1"
    THREADS="$2"
    BWA_INDEX="$3"
    PRIMER_SCHEMES_v3="$4"

    echo "${sample}"
    
    NAME=`basename "${sample}" `
    echo -e "\e[33m$NAME\033[0m"

    bwa mem -t ${THREADS} "${BWA_INDEX}" "${sample}"/*_R1_001.fastq.gz  "${sample}"/*_R2_001.fastq.gz | samtools sort | samtools view -F 4 -o "${sample}"/"${NAME}".sorted.bam
    
    ivar trim -e -i "${sample}"/"${NAME}".sorted.bam -b "${PRIMER_SCHEMES_v3}" -p "${sample}"/"${NAME}".primertrim
    samtools sort "${sample}"/"${NAME}".primertrim.bam -o "${sample}"/"${NAME}".primertrim.sorted.bam

    samtools mpileup -A -d 1000 -B -Q 0 --reference "${BWA_INDEX}" "${sample}"/"${NAME}".primertrim.sorted.bam | ivar consensus -q 13 -m 5 -p "${sample}"/"${NAME}".consensus -n N
    samtools index "${sample}"/"${NAME}".primertrim.sorted.bam
    mosdepth --by "${PRIMER_SCHEMES_v3}" "${sample}"/"${NAME}" "${sample}"/"${NAME}".primertrim.sorted.bam
    sed -n 2p "${sample}"/"${NAME}".mosdepth.summary.txt > "${sample}"/"${NAME}".coverage_report.txt
    awk '{print $0, "'${NAME}'"}' "${sample}"/"${NAME}".coverage_report.txt > "${sample}"/"${NAME}".coverage_report_name.txt

}



for sample in "${SEQ_FOLDER}"/demultiplex/*; do
    
    echo "${sample}"
    [ -d "${sample}" ] || continue 
    FILENAME=`basename "${sample}" `
    echo -e "\e[33m"${FILENAME}"\033[0m"

    mkdir -p "${SEQ_FOLDER}"/"results"/"${FILENAME}"
    echo "${FILENAME}"
    
    artic guppyplex \
	    --skip-quality-check \
	    --min-length 400 \
	    --max-length 1000 \
	    --directory "${SEQ_FOLDER}"/demultiplex/"${FILENAME}" \
	    --prefix "results"/"${FILENAME}"/"${FILENAME}"

    
    artic minion \
	    --normalise 200 \
	    --threads 150 \
        --scheme-directory "${PRIMER_SCHEMES_v3}" \
        --read-file "results"/"${FILENAME}"/"${FILENAME}"_"${FILENAME}".fastq \
        --fast5-directory "${FAST5}" \
        --sequencing-summary "${SEQ_FOLDER}"/"fastq_files"/sequencing_summary.txt \
        nCoV-2019/V3 \
        "results"/"${FILENAME}"/"${FILENAME}"

    #naredi "barcode2_barcode2.fastq" ime in nanopolish posledicno ni prepoznal zato sem popravil v vrstici zgoraj
    
    #artic minion --normalise 200 --threads 150 --scheme-directory "${PRIMER_SCHEMES_v3}" --read-file "${RESULT_DIR}"/"${FILENAME}"/"${FILENAME}.fastq --fast5-directory "${FAST5}" --sequencing-summary "${SEQ_FOLDER}"/sequencing_summary.txt nCoV-2019/V3 "${RESULT_DIR}"/"${FILENAME}"/"${FILENAME}"


export -f analysis

find "${SEQ_FOLDER}"/* -type d | parallel -j ${PARALLEL} analysis {} ${THREADS} "${BWA_INDEX}" "${PRIMER_SCHEMES_v3}"

cat "${SEQ_FOLDER}"/*/*.consensus.fa > "${SEQ_FOLDER}"/"${OUTPUT_NAME}".consensus_all.fasta
cat "${SEQ_FOLDER}"/*/*.coverage_report_name.txt > "${SEQ_FOLDER}"/"${OUTPUT_NAME}".coverage_report_all_tmp.txt
(echo "chrom    length  bases   mean    min max samplename" && cat "${SEQ_FOLDER}"/"${OUTPUT_NAME}".coverage_report_all_tmp.txt) >> "${SEQ_FOLDER}"/"${OUTPUT_NAME}".coverage_report_all.txt

echo -e "\e[33mResulting consensus sequence is in ${SEQ_FOLDER}/${NAME}.consensus_all.fasta\033[0m"

cat "${SEQ_FOLDER}"/"results"/*/*.consensus.fasta > "${SEQ_FOLDER}"/"results"/"${SEQ_FOLDER}".consensus_all.fasta
