#!/bin/bash
#SBATCH -J "mask_deviant"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=X
#SBATCH --time=1-00:00
#SBATCH --mem=10G

#Only in valeria
module load r

mkdir 02_info/mask_by_chr
MASK_LENGTH=150

Rscript 01_scripts/Rscript/deviant_masking.R $MASK_LENGTH
cat 02_info/mask_by_region/mask_deviant_chr*.bed > 02_info/mask_by_region_"$MASK_LENGTH"bp.bed

module load bedtools
module load samtools

#Mask reference genome
bedtools maskfasta -fi 02_info/genome.fasta \
    -bed 02_info/mask_by_region_"$MASK_LENGTH"bp.bed \
    -fo 02_info/genome_maskdev.fasta

samtools faidx 02_info/genome_maskdev.fasta

#Mask ancestral genome (optional)
#bedtools maskfasta -fi 02_info/ancestral_genome.fasta \
#    -bed 02_info/mask_by_region_"$MASK_LENGTH"bp.bed \
#    -fo 02_info/ancestral_genome_maskdev.fasta
#
#samtools faidx 02_info/ancestral_genome_maskdev.fasta