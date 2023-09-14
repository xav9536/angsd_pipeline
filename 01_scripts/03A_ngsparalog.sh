#!/bin/bash
#SBATCH -J "03A_ngsparalog"
#SBATCH -o log_%j
#SBATCH -c 4 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=20G

### This script will work on all bamfiles and use ANGSD to call SNPs that pass coverage and MAF filters,
### then calculate the likelihoods that reads wer mismapped at the position of snps using ngsparalog
### to produce list of canonical and deviant SNPs
### This process was parallelized by chromosomes for considerable gains in efficiency

#maybe edit
NB_CPU=4 #change accordingly in SLURM header, 
REGIONS="-rf 02_info/regions_25kb_100snp.txt" #optional edit with your region selected file
REGIONS="" # to remove the options to focus on a limited number of regions

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load angsd/0.931
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

echo " Calculate the SAF, MAF and GL for all individuals listed in 02_info/bam.filelist"
echo "keep loci with at least $MIN_DEPTH read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####Calculate the MAF and HWE
angsd -P $NB_CPU -nQueueSize 50 \
-doMaf 1 -doHWE 1 -GL 2 -doMajorMinor 1 -doCounts 1 \
-remove_bads 1 -minMapQ 30 -minQ 20 -skipTriallelic 1 -uniqueOnly 1\
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
-b 02_info/bam.filelist \
$REGIONS -out 03A_ngsparalog/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

#main features
#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
# -doMaf 1 (allele frequencies) -GL (Genotype likelihood 2 GATK method)
# -doHWE 1 (calculate deviation from Hardy-Weinberg Equilibrium)
# -doMajorMinor 1 use the most frequent allele as major
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ (minimum quality of reads?)
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 50%
#filter on allele frequency -minMaf, set to 0.05 

#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
#order the sites file by chromosome names 
#makes a region file matching the sites files and with same order
#index sites file
echo "from the maf file, extract a list of SNP chr, positoin, major all, minor all"
zcat 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".mafs.gz | cut -f1-4 | tail -n+2 | tee > 02_info/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND" 02_info/regions_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

OUTFILE_sites=02_info/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
OUTFILE_regions=02_info/regions_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

angsd sites index $OUTFILE_sites
