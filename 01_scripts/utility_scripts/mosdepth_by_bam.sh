#!/bin/bash
#SBATCH -J "mosdepth"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xavier.dallaire.2@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=10G

FIRST_CHR=$(head -n1 02_info/regions.txt)

mkdir 02_info/depth_chr01

## Run mosdepth in fast mode on the first chromosome for each
for i in $(seq $(wc -l 02_info/bam.filelist | cut -d' ' -f1))
do
ID=$(grep -o [A-Z][A-Z][A-Z]s_[0-9][0-9][0-9]-[0-9][0-9] ../02_info/bam.filelist | head -n $i | tail -n 1) ## Edit the pattern to fit your sample ID
FILE=$(head -n $i 02_info/bam.filelist | tail -n 1)
mosdepth -t 4 -n -x \
         -c $FIRST_CHR \ #only on first chromosome
         -i 2 \ #include only reads in proper pairs (by default in ANGSD)
         02_info/depth_chr01/$ID $FILE
done

## Create a summary of all bamfiles
for FILE in $(ls -1 02_info/depth_chr01/*.mosdepth.summary.txt)
do
NAME=`echo "$FILE" | cut -d'.' -f1`
DEPTH=`head -n2 $FILE` | tail -n1 | cut -f4`
echo -e $NAME"\t"$DEPTH >> 02_info/depth_chr01/summary_depth_chr01.txt
done
