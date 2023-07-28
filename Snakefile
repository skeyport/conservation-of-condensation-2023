# Jared Bard 230310
# first trime fastq files using trim-galore (removes sequencing adapters)
# then map poly(A) unstranded rna-seeq data using STAR aligner
# finally count reads using htseq-count, using CDS features and grouped by gene_id

sample = config['sample']
fastq1 = "fastq/"+config['fastq1']
fastq2 = "fastq/"+config['fastq2']
genome = config['genome']
gtf = "index/"+config['gtf']

rule all:
    input:
        "mapped_reads/"+sample+"/"+sample+"_Aligned_Sorted.out.bam.bai",
        "counts/"+sample+"_counts.tsv"
        
rule trim_galore:
    input:
        read1=fastq1,
        read2=fastq2,
    output:
        read1="trimmed/"+sample+"_1_val_1.fq.gz",
        read2="trimmed/"+sample+"_2_val_2.fq.gz",
    log:
        "logs/trim_galore/"+sample+"_trim_galore.log"
    conda:
        "envs/star_plus.yaml"
    threads: workflow.cores
    shell:
        "trim_galore --illumina --fastqc_args '--outdir fastqc/' "
        "--paired -o trimmed -j {threads} "
        "{input.read1} {input.read2} &> {log}"

rule STAR_index:
    input:
        gtf=gtf,
        fasta="index/"+genome+".fasta"
    output:
# clever trick to create an artifical output that is required as an input for mapping
# this forces snakemake to run the index first
        touch("index/index_"+genome+".done")
    log:
        "logs/STAR_index/"+genome+"_star_index.log"
    params:
        genome=genome,
        sjdb_params="--sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentGene gene_id"
    conda:
        "envs/star_plus.yaml"
    shell:
        r"""
        (STAR --runMode genomeGenerate --sjdbOverhang 99 {params.sjdb_params} --genomeSAindexNbases 10 \
        --runThreadN {threads} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} \
        --genomeDir index/{params.genome}) >& {log}
        """
        
rule STAR_map_paired:
    input:
        gtf=gtf,
        idxdone = "index/index_"+genome+".done",
        read1="trimmed/"+sample+"_1_val_1.fq.gz",
        read2="trimmed/"+sample+"_2_val_2.fq.gz",
    output:
        temp("mapped_reads/"+sample+"/"+sample+"_Aligned.out.bam"),
        "mapped_reads/"+sample+"/"+sample+"_Log.final.out"
    params:
        out_prefix="mapped_reads/"+sample+"/"+sample+"_",
        genome = "index/"+config["genome"],
        zip_params = "--readFilesCommand gunzip -c ",
        sjdb_params="--sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentGene gene_id"
    log:
        "logs/STAR_map/"+sample+"_star_map.log"
    threads:
        workflow.cores
    conda:
        "envs/star_plus.yaml"
    shell:
        r"""
        STAR --outSAMtype BAM Unsorted \
        {params.zip_params} \
        --sjdbGTFfile {input.gtf} {params.sjdb_params} \
        --runThreadN {threads} --alignMatesGapMax 20000 --limitBAMsortRAM 1445804817 \
        --genomeDir {params.genome} --outFileNamePrefix {params.out_prefix} \
        --readFilesIn {input.read1} {input.read2} >& {log}
        """

rule samtools_sort:
    input:
        "mapped_reads/"+sample+"/"+sample+"_Aligned.out.bam"
    output:
        "mapped_reads/"+sample+"/"+sample+"_Aligned_Sorted.out.bam"
    log:
        "logs/samtools_sort/"+sample+"_samtools_sort.log"
    threads:
        workflow.cores
    conda:
        "envs/star_plus.yaml"
    shell:
        "samtools sort -@ {threads} {input} -o {output} >& {log}"

rule samtools_index:
    input:
        "mapped_reads/"+sample+"/"+sample+"_Aligned_Sorted.out.bam"
    output:
        "mapped_reads/"+sample+"/"+sample+"_Aligned_Sorted.out.bam.bai"
    threads:
        workflow.cores
    conda:
        "envs/star_plus.yaml"
    shell:
        "samtools index -@ {threads} {input}"

rule htseq:
    input:
        bam="mapped_reads/"+sample+"/"+sample+"_Aligned.out.bam",
        gtf=gtf
    output:
        "counts/"+sample+"_counts.tsv"
    log:
        "logs/htseq/"+sample+"/"+sample+"_htseq.log"
    params:
        "--stranded=no --type=CDS --idattr=gene_id" 
    conda:
        "envs/star_plus.yaml"
    shell:
        "htseq-count -c {output} {params} {input.bam} {input.gtf} >& {log}"

        
