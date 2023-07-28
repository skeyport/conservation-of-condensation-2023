#!/bin/bash
#SBATCH --job-name=SK10
#SBATCH --output=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_SK10_%j.out
#SBATCH --error=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_SK10_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=pi-dadrummond
#SBATCH --mem=32G

module load python/anaconda-2022.05
source activate /home/jbard/beagle3-dadrummond/jbard/envs/snake38
cd /home/jbard/midway3-scratch/SK1_SK12/
snakemake --config sample=SK10 \
fastq1="SK10_1.fq.gz" \
fastq2="SK10_2.fq.gz" \
genome=GCA_001417885.1_Kmar_1.0_genomic \
gtf=GCA_001417885.1_Kmar_1.0_genomic.gtf \
--cores all -p --use-conda
