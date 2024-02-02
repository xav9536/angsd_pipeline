#!/bin/bash
#SBATCH -J "05_ld_decay"
#SBATCH -o LD_decay
#SBATCH -c 1
#SBATCH -p ibis_medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=7-00:00
#SBATCH --mem=200G

#Module to load for Valeria
module load gsl StdEnv/2018.3 gcc/7.3.0 rstudio-server/1.2.1335

#module to load on Manitou
#module load miniconda/3-py3.10

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
BAMLIST="02_info/bam.filelist"
N_IND=$(wc -l "$BAMLIST" | cut -d " " -f 1)

### Listing the ld.gz files. If you have many chromosome, considering plotting only the first 10
ls -1 03C_LDpruned/*.ld.gz | head -n 10 > 03C_LDpruned/list.ldfiles

Rscript ~/program_valeria/ngsLD/scripts/fit_LDdecay.R \
    --ld_files 03C_LDpruned/list.ldfiles \
    --ld r2 --n_ind N_IND \
    --max_kb_dist 500 \
    --fit_level 100 --fit_boot 100 \
    --plot_data --plot_no_legend \
    -o 03C_LDpruned/ld_decay_10chrs.pdf
