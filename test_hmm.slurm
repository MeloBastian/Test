#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3gb
#SBATCH -o 44443at40674_filtred_ali.fna_hmm.out


singularity exec /beegfs/home/tricou/singularity_image/hmm.sif HmmCleaner.pl 44443at40674_filtred_ali.fna