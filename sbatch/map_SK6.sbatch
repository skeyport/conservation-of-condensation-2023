#!/bin/bash
#SBATCH --job-name=SK6
#SBATCH --output=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_SK6_%j.out
#SBATCH --error=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_SK6_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=pi-dadrummond
#SBATCH --mem=32G

module load python/anaconda-2022.05
source activate /home/jbard/beagle3-dadrummond/jbard/envs/snake38
cd /home/jbard/midway3-scratch/SK1_SK12/
snakemake --config sample=SK6 \
fastq1="SK6_1.fq.gz" \
fastq2="SK6_2.fq.gz" \
genome=GCF_000146045.2_R64_genomic \
gtf=GCF_000146045.2_R64_genomic.gtf \
--cores all -p --use-conda
