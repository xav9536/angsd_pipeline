#!/bin/bash

#To run on VALERIA:
#cat ./02_info/regions_number.txt | parallel -j20 srun -p ibis_medium -c 1 --mem=20G --time=7-00:00 -o log_%j 01_scripts/03C3_LDpruning.sh {} 200000 0.1

#module to load in Valeria
module load StdEnv/2020 gcc/9.3.0 python/3.8 graph-tool/2.45

#module to load on Manitou
#module load miniconda/3-py3.10

REGION_NUM="$1"
DIST=$2
WEIGHT=$3

#Skipping header and pruning graph on LD
zcat 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".ld.gz | tail -n +2 | grep -v -i nan  > 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".ld.nh

~/program_valeria/ngsLD/scripts/prune_ngsLD.py \
    --input 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".ld.nh \
    --max_dist $DIST --min_weight $WEIGHT --field_dist 3 --field_weight 7 --weight_precision 4 \
    --output 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".pruned

rm 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".ld.nh

#Replace ':' by tab in output
sed -i -e 's/:/\t/g' 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".pruned

#Subset beagle files from 03B
python 01_scripts/utility_scripts/beagle_extract_wanted_snp.py \
    03B_gl_maf_canonical/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".beagle.gz \
    03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".pruned \
    03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".pruned.beagle.gz