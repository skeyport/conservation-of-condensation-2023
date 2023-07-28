#!/bin/bash
#SBATCH --job-name={SAMPLE}
#SBATCH --output=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_{SAMPLE}_%j.out
#SBATCH --error=/home/jbard/midway3-scratch/SK1_SK12/sbatch/map_{SAMPLE}_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=pi-dadrummond
#SBATCH --mem=32G

module load python/anaconda-2022.05
source activate /home/jbard/beagle3-dadrummond/jbard/envs/snake38
cd /home/jbard/midway3-scratch/SK1_SK12/
snakemake --config sample={SAMPLE} \
fastq1="{FASTQ1}" \
fastq2="{FASTQ2}" \
genome={GENOME} \
gtf={GTF} \
--cores all -p --use-conda
