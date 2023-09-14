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

REGION_NUM="$1" 
REGION=$(head -n $REGION_NUM 02_info/regions.txt | tail -n 1)

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load angsd/0.931
module load samtools
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
-remove_bads 1 -minMapQ 30 -minQ 20 -skipTriallelic 1 \
-uniqueOnly 1 -only_proper_pairs 1 \
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
-b 02_info/bam.filelist \
-r $REGION -out 03A_ngsparalog/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"

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
#filter on reads with several hits -uniqueOnly
#filter on pairs of reads not properly mapped -only_proper_pairs 1 (by default in ANGSD)



#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
#output
#index sites file
echo "from the maf file, extract a list of SNP chr, positoin, major all, minor all"
mkdir 02_info/sites_all_by_chr/
gunzip 03A_ngsparalog/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".mafs.gz 

INFILE=03A_ngsparalog/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".mafs
OUTFILE_sites=02_info/sites_all_by_chr/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM"
BEDFILE=02_info/sites_all_by_chr/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".bed

Rscript 01_scripts/Rscripts/make_sites_list_maxdepth_simple.R "$INFILE" "$OUTFILE_sites"

angsd sites index $OUTFILE_sites

# convert site file to bed for ngsparalog
awk '{print $1"\t"$2-1"\t"$2}' $OUTFILE_sites > $BED_FILE

### mpileup and ngsParalog without intermidate files
samtools mpileup -b 02_info/bam.filelist -l $BED_FILE -r $REGION -q 0 -Q 0 --ff UNMAP,DUP |
~/programs/ngsParalog/ngsParalog calcLR \
    -infile - \
    -outfile 03A_ngsparalog/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_chr"$REGION_NUM".ngsparalog \
    -minQ 20 -minind $MIN_IND -mincov $MIN_DEPTH -allow_overwrite 1

### Convert ngsparalog output in list of canonical and deviant SNPs based on p-value threshold


