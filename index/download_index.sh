#!/bin/bash
Skud_fasta="GCA_947243785.1_Skud-ZP591_genomic.fna.gz"
curl -C - -o $Skud_fasta https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Saccharomyces_kudriavzevii/latest_assembly_versions/GCA_947243785.1_Skud-ZP591/${Skud_fasta}
gunzip -c $Skud_fasta > "GCA_947243785.1_Skud-ZP591_genomic.fasta"
Skud_gtf="GCA_947243785.1_Skud-ZP591_genomic.gtf.gz"
curl -C - -o $Skud_gtf https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Saccharomyces_kudriavzevii/latest_assembly_versions/GCA_947243785.1_Skud-ZP591/${Skud_gtf}
gunzip $Skud_gtf


Kmarx_fasta="GCA_001417885.1_Kmar_1.0_genomic.fna.gz"
curl -C - -o $Kmarx_fasta https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Kluyveromyces_marxianus/latest_assembly_versions/GCA_001417885.1_Kmar_1.0/${Kmarx_fasta}
gunzip -c $Kmarx_fasta > "GCA_001417885.1_Kmar_1.0_genomic.fasta"
Kmarx_gtf="GCA_001417885.1_Kmar_1.0_genomic.gtf.gz"
curl -C - -o $Kmarx_gtf https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Kluyveromyces_marxianus/latest_assembly_versions/GCA_001417885.1_Kmar_1.0/${Kmarx_gtf}
gunzip $Kmarx_gtf

Scer_fasta="GCF_000146045.2_R64_genomic.fna.gz"
curl -C - -o $Scer_fasta https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/${Scer_fasta}
gunzip -c $Scer_fasta > GCF_000146045.2_R64_genomic.fasta
Scer_gtf="GCF_000146045.2_R64_genomic.gtf.gz"
curl -C - -o $Scer_gtf https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/${Scer_gtf}
gunzip $Scer_gtf