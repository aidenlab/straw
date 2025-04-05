#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <path_to_straw_executable>"
    exit 1
fi

STRAW_EXE="$1"
LOG_FILE="straw_test_results.log"

echo "Starting tests with straw executable: $STRAW_EXE" | tee "$LOG_FILE"
echo "Timestamp: $(date)" | tee -a "$LOG_FILE"
echo "----------------------------------------" | tee -a "$LOG_FILE"

norms=("NONE" "VC")
files=('/Users/muhammad/Desktop/hg38/GM12878_30.hic')
      #'https://www.encodeproject.org/files/ENCFF757POS/@@download/ENCFF757POS.hic')
#
chroms=("chr2" "chr5")
mtypes=("BP" "MATRIX")
for mtype in "${mtypes[@]}"
do
for chrom in "${chroms[@]}"
do
for norm in "${norms[@]}"
do
for file in "${files[@]}"
do
    echo "Testing configuration:" | tee -a "$LOG_FILE"
    echo "  File: $file" | tee -a "$LOG_FILE"
    echo "  Norm: $norm" | tee -a "$LOG_FILE"
    echo "  Chrom: $chrom" | tee -a "$LOG_FILE"
    echo "  Type: $mtype" | tee -a "$LOG_FILE"

    # Whole chromosome test
    echo "Whole chromosome test:" | tee -a "$LOG_FILE"
    "$STRAW_EXE" observed "$norm" "$file" "$chrom" "$chrom" "$mtype" 1000000 > output.bin
    MD5=$(md5sum output.bin | cut -d' ' -f1)
    echo "  MD5: $MD5" | tee -a "$LOG_FILE"

    # On-diagonal region test
    echo "On-diagonal region test:" | tee -a "$LOG_FILE"
    "$STRAW_EXE" observed "$norm" "$file" "${chrom}:10000000:50000000" "${chrom}:10000000:50000000" "$mtype" 100000 > output.bin
    MD5=$(md5sum output.bin | cut -d' ' -f1)
    echo "  MD5: $MD5" | tee -a "$LOG_FILE"

    # Off-diagonal region test
    echo "Off-diagonal region test:" | tee -a "$LOG_FILE"
    "$STRAW_EXE" observed "$norm" "$file" "${chrom}:10000000:50000000" "${chrom}:40000000:80000000" "$mtype" 100000 > output.bin
    MD5=$(md5sum output.bin | cut -d' ' -f1)
    echo "  MD5: $MD5" | tee -a "$LOG_FILE"

    # Inter-chromosomal test
    echo "Inter-chromosomal test:" | tee -a "$LOG_FILE"
    "$STRAW_EXE" observed "$norm" "$file" "${chrom}:10000000:50000000" "chr10:10000000:30000000" "$mtype" 100000 > output.bin
    MD5=$(md5sum output.bin | cut -d' ' -f1)
    echo "  MD5: $MD5" | tee -a "$LOG_FILE"

    echo "----------------------------------------" | tee -a "$LOG_FILE"
done
done
done
done

rm output.bin
echo "Test results have been logged to $LOG_FILE"
