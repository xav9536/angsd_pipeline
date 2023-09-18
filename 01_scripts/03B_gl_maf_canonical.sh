#!/bin/bash

### This script will work on all bamfiles and use ANGSD to call SNPs from the previously create list of canonical SNPs,
### and output genotype likelihood (beagle.gz) and MAF (maf.gz)

### This process was parallelized by chromosomes for considerable gains in efficiency
### For example, to run 10 chromosomes at the same time, use:
### cat 02_info/regions_number.txt | parallel -j10 srun -c 4 --mem 20G -p ibis_small -o log_%j --time 1-00:00 ./01_scripts/03B_gl_maf_canonical.sh {}
### Adjust -j and -c to fit your available ressources

#maybe edit
NB_CPU=4 #change accordingly in SLURM header

REGION_NUM="$1" 
REGION=$(head -n $REGION_NUM 02_info/regions.txt | tail -n 1)
BAMLIST="02_info/bam.filelist"

module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

echo " Calculate the MAF and GL for all individuals listed in 02_info/bam.filelist" on canonical SNP sites
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####Calculate the MAF and GL
angsd -P $NB_CPU -nQueueSize 50 \
-domaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 \
-anc 02_info/genome.fasta \
-sites 02_info/sites_by_chr/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"_canonical \
-remove_bads 1 -minMapQ 30 -minQ 20 -skipTriallelic 1 \
-uniqueOnly 1 -only_proper_pairs 1 \
-r "$REGION" \
-b "$BAMLIST" \
-out 03B_gl_maf_canonical/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"

#main features
#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 2 GATK method - export GL in beagle format  -doGLF2) 
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case -fold 1 (car on utilise la ref comme ancestral
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

