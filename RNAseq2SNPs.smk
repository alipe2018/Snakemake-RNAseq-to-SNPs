import yaml

##========= Globals ==========
configfile: 'config.yaml'

## Set samples
FILES = yaml.load(open(config['SampleYaml']))
SAMPLES = sorted(FILES.keys())

## set output Dir 
WORKDIR = config["WorkDir"]

## Set reference file
DNA = config["dna"]
GTF = config["gtf"]
# HISAT2_INDEX_PREFIX = config["GenomeHisat2Index"]

##======== Rules ============
rule all:
    input:
#       expand( WORKDIR + "Step0.Prepare/{sample}.R1.fq.gz", sample=SAMPLES ),
#       expand( WORKDIR + "Step0.Prepare/{sample}.R2.fq.gz", sample=SAMPLES ),
#       expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R1.fq.gz", sample=SAMPLES ),
#       expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R2.fq.gz", sample=SAMPLES ),
#       expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.json", sample=SAMPLES ),
#       expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.html", sample=SAMPLES ),
#       expand( WORKDIR + "Step2.hisat2Align/{sample}.bam", sample=SAMPLES),
        expand( WORKDIR + "Step2.hisat2Align/{sample}.xls", sample=SAMPLES),
        expand( WORKDIR + "Step2.sorted/{sample}.sorted.bam", sample=SAMPLES),
        expand( WORKDIR + "Step2.sorted/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand( WORKDIR + "Step2.sorted/{sample}.sorted_rmdup.bam", sample=SAMPLES),
        expand( WORKDIR + "Step2.sorted/{sample}.sorted_rmdup.mtx", sample=SAMPLES),
        WORKDIR + "Step3.SNP_Calling/bam.list",
        WORKDIR + "Step3.SNP_Calling/calls.bcf"
 

## ======== Step 0  Prepare Rename the fastq file ======== 
rule ReName:
    input:
        lambda wildcards:FILES[wildcards.sample]
    output:
       R1 =  WORKDIR + "Step0.Prepare/{sample}.R1.fq.gz",
       R2 =  WORKDIR + "Step0.Prepare/{sample}.R2.fq.gz"
    run:
        shell("ln -sf {input[0]} {output.R1}")
        shell("ln -sf {input[1]} {output.R2}")


##======== Step1 raw fastq filter ========
rule FastqFilter:
    input:
        R1 =  WORKDIR + "Step0.Prepare/{sample}.R1.fq.gz",
        R2 =  WORKDIR + "Step0.Prepare/{sample}.R2.fq.gz"
    output:
        R1 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R1.fq.gz",
        R2 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R2.fq.gz",
        json = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.json",
        html = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.html"
    log:
        WORKDIR + "logs/Step1.fastqFilter/{sample}.fastqFilter.log"
    benchmark:
        WORKDIR + "benchmark/Step1.fastqFilter/{sample}.benchmark"
    params:
        "--detect_adapter_for_pe --qualified_quality_phred 20 --unqualified_percent_limit 5 --n_base_limit 5 --length_required 50 --correction"
    threads:
        8
    run:
        shell("/MaizeShen/liupeng/Test/soft/fastp {params} -w {threads} -i {input.R1} -I {input.R2} -o {output.R1} "
        "-O {output.R2} -j {output.json} -h {output.html} 2> {log}")


## ======== Step2 hisat alignment ========
rule hisat2Align:
    input:
        R1 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R1.fq.gz",
        R2 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R2.fq.gz"
    output:
        bam = WORKDIR + "Step2.hisat2Align/{sample}.bam",
        sum = WORKDIR + "Step2.hisat2Align/{sample}.xls"
    log:
        WORKDIR + "logs/Step2.hisat2Align/{sample}.align.logs"
    benchmark:
        WORKDIR + "benchmark/Step2.hisat2Align/{sample}.benchmark"
    threads:
        8
    params:
        "--rg-id {sample} --rg {sample}"
    shell:
         " hisat2 {params} -p {threads} -x {HISAT2_INDEX_PREFIX} -1 {input.R1} -2 {input.R2} "
         " --summary-file {output.sum} 2> {log} | samtools view -Sb -@ {threads} -m 5G -o {output.bam} "

## ======== Step 2.2 sort bam ========
rule sortBam:
    input:
        bam = WORKDIR + "Step2.hisat2Align/{sample}.bam"
    output:
        # protected(WORKDIR + "Step2.sorted/{sample}.sorted.bam")
        WORKDIR + "Step2.sorted/{sample}.sorted.bam"
    threads:
        8
    resources: 
        mem_mb=5000
    shell:
        "samtools sort -@ {threads} -m 5G -o {output} {input.bam}"

## ======== Step 2.3 bam index ========
rule bamIndex:
    input:
        WORKDIR + "Step2.sorted/{sample}.sorted.bam"
    output:
        WORKDIR + "Step2.sorted/{sample}.sorted.bam.bai"
    threads:
        8
    shell:
        "samtools index -b -@ {threads} {input} {output}"

## ======== Step 2.4 Mark Duplicates ========
rule MarkDuplicates:
    input:
        WORKDIR + "Step2.sorted/{sample}.sorted.bam"
    output:
        bam = WORKDIR + "Step2.sorted/{sample}.sorted_rmdup.bam",
        mtx = WORKDIR + "Step2.sorted/{sample}.sorted_rmdup.mtx"
    threads:
        2
    resources:
        mem_mb=10000
    params:
        "ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true"
    shell:
        "picard MarkDuplicates {params} INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.mtx}"

## ======== Step 3.1 make bam list =======
rule MakeBamList:
    input:
        expand(WORKDIR + "Step2.sorted/{sample}.sorted_rmdup.bam", sample = SAMPLES)
    output:
        WORKDIR + "Step3.SNP_Calling/bam.list"
    run:
        with open output[0] as f:
            for bam in {input}:
                print(bam, file=f)

## ======= Step 3.2 SNP Calling =======
rule CallSNP:
    input:
        list = WORKDIR + "Step3.SNP_Calling/bam.list"
    output:
        bcf = WORKDIR + "Step3.SNP_Calling/calls.bcf"
    threads:
        16
    params:
        mp = "-C 50",
        ca = "-mv -Ob"
    shell:
        "bcftools mpileup {params.mp} -b {input.list} -f {DNA} --threads {threads} | bcftools call {params.ca} -o {output.bcf}"
