#!/bin/bash
#projdir="."
projdir="/home/jbard/midway3-scratch/SK1_SK12"

declare -a Skud=("SK1" "SK2" "SK3" "SK4")
declare -a Scer=("SK5" "SK6" "SK7" "SK8")
declare -a Kmarx=("SK9" "SK10" "SK11" "SK12")
Skud_genome="GCA_947243785.1_Skud-ZP591_genomic"
Skud_gtf="GCA_947243785.1_Skud-ZP591_genomic.gtf"
Scer_genome="GCF_000146045.2_R64_genomic"
Scer_gtf="GCF_000146045.2_R64_genomic.gtf"
Kmarx_genome="GCA_001417885.1_Kmar_1.0_genomic"
Kmarx_gtf="GCA_001417885.1_Kmar_1.0_genomic.gtf"

template="snake_map_template.sbatch"

## now loop through the above arrays
for sample in "${Skud[@]}"
do
    fastq1="${sample}_1.fq.gz"
    fastq2="${sample}_2.fq.gz"
    genome=$Skud_genome
    gtf=$Skud_gtf

    sed s@'{SAMPLE}'@"${sample}"@g ${projdir}/${template} > ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{FASTQ1}'@"${fastq1}"@g ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{FASTQ2}'@"${fastq2}"@g ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{GENOME}'@"${genome}"@g ${projdir}/sbatch/map_${sample}.sbatch
	sed -i s@'{GTF}'@"${gtf}"@g ${projdir}/sbatch/map_${sample}.sbatch
done

for sample in "${Scer[@]}"
do
    fastq1="${sample}_1.fq.gz"
    fastq2="${sample}_2.fq.gz"
    genome=$Scer_genome
    gtf=$Scer_gtf

    sed s@'{SAMPLE}'@"${sample}"@g ${projdir}/${template} > ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{FASTQ1}'@"${fastq1}"@g ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{FASTQ2}'@"${fastq2}"@g ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{GENOME}'@"${genome}"@g ${projdir}/sbatch/map_${sample}.sbatch
	sed -i s@'{GTF}'@"${gtf}"@g ${projdir}/sbatch/map_${sample}.sbatch
done

for sample in "${Kmarx[@]}"
do
    fastq1="${sample}_1.fq.gz"
    fastq2="${sample}_2.fq.gz"
    genome=$Kmarx_genome
    gtf=$Kmarx_gtf

    sed s@'{SAMPLE}'@"${sample}"@g ${projdir}/${template} > ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{FASTQ1}'@"${fastq1}"@g ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{FASTQ2}'@"${fastq2}"@g ${projdir}/sbatch/map_${sample}.sbatch
    sed -i s@'{GENOME}'@"${genome}"@g ${projdir}/sbatch/map_${sample}.sbatch
	sed -i s@'{GTF}'@"${gtf}"@g ${projdir}/sbatch/map_${sample}.sbatch
done