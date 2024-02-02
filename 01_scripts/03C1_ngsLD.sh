#!/bin/bash

cd $SLURM_SUBMIT_DIR

#To run in parallel
#parallel -a ./02_info/regions_number.txt -j 10 srun -p ibis_medium -c 5 --mem=20G --time=7-00:00  -o 05_ngsLD_per_chr_{}_%j.log 01_scripts/03C1_ngsLD.sh {} &

#module to load on Valeria
module load gsl StdEnv/2018.3 gcc/7.3.0 rstudio-server/1.2.1335

#module to load on Manitou
#module load miniconda/3-py3.10


#Important ngsLD parameters:

#--probs: specification of whether the input is genotype probabilities (likelihoods or posteriors)?
#--n_ind INT: sample size (number of individuals).
#--n_sites INT: total number of sites.
#--max_kb_dist DOUBLE: maximum distance between SNPs (in Kb) to calculate LD. Set to 0(zero) to disable filter. [100]
#--max_snp_dist INT: maximum distance between SNPs (in number of SNPs) to calculate LD. Set to 0 (zero) to disable filter. [0]
#--rnd_sample= sample x% of comparison between snp
#--n_threads INT: number of threads to use. [1]
#--out FILE: output file name. [stdout]
#
REGION_NUM="$1" 
REGION=$(head -n $REGION_NUM 02_info/regions.txt | tail -n 1)
NGSLD=01_scripts/ngsLD
BAMLIST="02_info/bam.filelist"
NCPU=5

mkdir 03C_LDpruned

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l "$BAMLIST" | cut -d " " -f 1)

NSITES=$(wc -l 02_info/sites_by_chr/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"_canonical)


sed -e 's/ /\t/g' 02_info/sites_by_chr/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"_canonical > 02_info/sites_by_chr/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"_canonical.tab

~/program_valeria/ngsLD/ngsLD \
    --geno 03B_gl_maf_canonical/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".beagle.gz \
	--probs \
	--n_ind "$NIND" \
	--n_sites "$NSITES" \
	--min_maf "$MIN_MAF" \
    --pos 02_info/sites_by_chr/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"_canonical.tab \
	--max_kb_dist 500 \
	--rnd_sample 0.5 \
	--n_threads "$NCPU" \
	--out 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".ld 
    
gzip 03C_LDpruned/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".ld 

