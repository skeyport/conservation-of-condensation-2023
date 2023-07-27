#!/bin/bash
#
declare -a Skud=("SK1" "SK2" "SK3" "SK4")
declare -a Scer=("SK5" "SK6" "SK7" "SK8")
declare -a Kmarx=("SK9" "SK10" "SK11" "SK12")
for sample in "${Skud[@]}"; do
	sbatch sbatch/map_${sample}.sbatch
done
for sample in "${Scer[@]}"; do
	sbatch sbatch/map_${sample}.sbatch
done
for sample in "${Kmarx[@]}"; do
	sbatch sbatch/map_${sample}.sbatch
done

