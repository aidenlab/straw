
echo "test 1, no norms"

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
    echo "."
    echo " file: $file"
    echo -n " norm: $norm"
    echo -n " chrom: $chrom"
    echo -n " mtype: $mtype"
    test_type="whole chrom data dump test"

    ./cmake-build-debug/straw observed ${norm} ${file} ${chrom} ${chrom} ${mtype} 1000000 > new.bin
    ./straw observed ${norm} ${file} ${chrom} ${chrom} ${mtype} 1000000 > prev.bin

    echo -n "."
    cmp --silent new.bin prev.bin && echo -n '.' || echo '### ERROR: Straw versions pulling different results! ### '${test_type}

    test_type="on diag region dump test"

    ./cmake-build-debug/straw observed ${norm} ${file} ${chrom}:10000000:50000000 ${chrom}:10000000:50000000 ${mtype} 100000 > new.bin
    ./straw observed ${norm} ${file} ${chrom}:10000000:50000000 ${chrom}:10000000:50000000 ${mtype} 100000 > prev.bin

    echo -n "."
    cmp --silent new.bin prev.bin && echo -n '.' || echo '### ERROR: Straw versions pulling different results! ### '${test_type}

    test_type="off diag region dump test"

    ./cmake-build-debug/straw observed ${norm} ${file} ${chrom}:10000000:50000000 ${chrom}:40000000:80000000 ${mtype} 100000 > new.bin
    ./straw observed ${norm} ${file} ${chrom}:10000000:50000000 ${chrom}:40000000:80000000 ${mtype} 100000 > prev.bin

    echo -n "."
    cmp --silent new.bin prev.bin && echo -n '.' || echo '### ERROR: Straw versions pulling different results! ### '${test_type}

    test_type="inter chrom diag region vs chr10 dump test"

    ./cmake-build-debug/straw observed ${norm} ${file} ${chrom}:10000000:50000000 chr10:10000000:30000000 ${mtype} 100000 > new.bin
    ./straw observed ${norm} ${file} ${chrom}:10000000:50000000 chr10:10000000:30000000 ${mtype} 100000 > prev.bin

    echo -n "."
    cmp --silent new.bin prev.bin && echo -n '.' || echo '### ERROR: Straw versions pulling different results! ### '${test_type}

done
done
done
done
