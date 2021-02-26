#!/usr/bin/env bash

#===============================================================================
#
#          FILE:  run_illumina_covid_seq.sh
#
#         USAGE: conda enviroment with ivar needs to be activated		
#
#       VERSION:  1.0.0
#       CREATED:  26/02/2021
#===============================================================================

print_usage() {
  printf "Usage: 
        conda activate covid_illumina
        run_illumina_covid_seq_parallel.sh [Options] Path_to_folder_with sequences\n
        Options: 
        -n [THREADS] number of threads used [Default: 20]
        -j [PARALLEL] number of parallel runs [Default: 12]
        -o [OUTPUT_NAME] specify the name of resulting consencuc sequence [Default: Result_seq]"
}

# Default values
THREADS=20
OUTPUT_NAME="Result_seq"
PARALLEL=10

# Options
while getopts 'n:j:o:h' OPTION; do
  case "$OPTION" in
    n) THREADS="$OPTARG" ;;
    j) PARALLEL="$OPTARG" ;;
    o) OUTPUT_NAME="$OPTARG" ;;
    h) print_usage
    exit 1 ;;
  esac
done

shift "$(($OPTIND -1))"

if [ -z "$*" ]; then print_usage && exit 1; fi
SEQ_FOLDER=("$@")

BWA_INDEX="/home/kisld/bioinformatic_tools/artic-ncov2019/covid19_bwa_index/nCoV-2019.reference.fasta"
PRIMER_SCHEMES_v3="/home/kisld/bioinformatic_tools/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.primer.bed"

# source activate /home/kisld/anaconda3/envs/covid_illumina
# conda activate covid_illumina

for file in "${SEQ_FOLDER}"/*; do

    [ -f "${file}" ] || continue

    fastq_name=`basename "${file}" `
    echo "${fastq_name}"
    sample_name=$(echo $"${fastq_name}" | cut -d '_' -f 1)
    mkdir -p "${SEQ_FOLDER}"/"${sample_name}"
    mv "${file}" "${SEQ_FOLDER}"/"${sample_name}"

done

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

export -f analysis

find "${SEQ_FOLDER}"/* -type d | parallel -j ${PARALLEL} analysis {} ${THREADS} "${BWA_INDEX}" "${PRIMER_SCHEMES_v3}"

cat "${SEQ_FOLDER}"/*/*.consensus.fa > "${SEQ_FOLDER}"/"${OUTPUT_NAME}".consensus_all.fasta
cat "${SEQ_FOLDER}"/*/*.coverage_report_name.txt > "${SEQ_FOLDER}"/"${OUTPUT_NAME}".coverage_report_all_tmp.txt
(echo "chrom    length  bases   mean    min max samplename" && cat "${SEQ_FOLDER}"/"${OUTPUT_NAME}".coverage_report_all_tmp.txt) >> "${SEQ_FOLDER}"/"${OUTPUT_NAME}".coverage_report_all.txt

echo -e "\e[33mResulting consensus sequence is in ${SEQ_FOLDER}/${NAME}.consensus_all.fasta\033[0m"
