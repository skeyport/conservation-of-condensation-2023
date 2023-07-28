#!/bin/bash
#SBATCH --job-name=SK3
#SBATCH --output=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_SK3_%j.out
#SBATCH --error=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_SK3_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=pi-dadrummond
#SBATCH --mem=32G

module load python/anaconda-2022.05
source activate /home/jbard/beagle3-dadrummond/jbard/envs/snake38
cd /home/jbard/midway3-scratch/SK1_SK12/
snakemake --config sample=SK3 \
fastq1="SK3_1.fq.gz" \
fastq2="SK3_2.fq.gz" \
genome=GCA_947243785.1_Skud-ZP591_genomic \
gtf=GCA_947243785.1_Skud-ZP591_genomic.gtf \
--cores all -p --use-conda
