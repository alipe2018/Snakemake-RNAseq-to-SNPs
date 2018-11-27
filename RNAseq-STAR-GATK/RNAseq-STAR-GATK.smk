import yaml
import os
import re

## ========= Globals ==========
configfile: 'config.yaml'

FILES = yaml.load(open(config['SampleYaml']))
SAMPLES = sorted(FILES.keys())
WORKDIR = config["WorkDir"]
DNA = config["dna"]
GTF = config["gtf"]


## ======== Begin Snakemake Rules ============
rule all:
    input:
        # expand(WORKDIR + "Step1.FastqFilter/{sample}.R1.fq.gz", sample=SAMPLES),
        # expand(WORKDIR + "Step1.FastqFilter/{sample}.R2.fq.gz", sample=SAMPLES),
        # expand(WORKDIR + "Step1.FastqFilter/{sample}.html", sample=SAMPLES),
        # expand(WORKDIR + "Step1.FastqFilter/{sample}.json", sample=SAMPLES),
        WORKDIR + "Step2.StarAlign/star-1-index/sjdbList.out.tab",
        expand(WORKDIR + "Step2.StarAlign/star-1-pass/{sample}SJ.out.tab", sample=SAMPLES),
        WORKDIR + "Step2.StarAlign/SJ.Pass-1-Merged.tab",
        WORKDIR + "Step2.StarAlign/star-2-index/sjdbList.out.tab",
        expand(WORKDIR + "Step2.StarAlign/star-2-pass/{sample}SJ.out.tab", sample=SAMPLES),
        expand(WORKDIR + "Step2.StarAlign/star-2-pass/{sample}Aligned.sortedByCoord.out.bam", sample=SAMPLES)


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

## ======== Step 1 Fastq Filter ========
rule FastqFilter:
    input:
        R1 =  WORKDIR + "Step0.Prepare/{sample}.R1.fq.gz",
        R2 =  WORKDIR + "Step0.Prepare/{sample}.R2.fq.gz"
    output:
        R1 = WORKDIR + "Step1.FastqFilter/{sample}.R1.fq.gz",
        R2 = WORKDIR + "Step1.FastqFilter/{sample}.R2.fq.gz",
        json = WORKDIR + "Step1.FastqFilter/{sample}.json",
        html = WORKDIR + "Step1.FastqFilter/{sample}.html"
    log:
        WORKDIR + "logs/Step1.FastqFilter/{sample}.fastqFilter.log"
    benchmark:
        WORKDIR + "benchmark/Step1.FastqFilter/{sample}.benchmark.tsv"
    params:
        " --detect_adapter_for_pe --qualified_quality_phred 20 "
        " --unqualified_percent_limit 5 --n_base_limit 5 "
        " --length_required 50 --correction "
    threads:
        8
    run:
        shell("/MaizeShen/liupeng/Test/soft/fastp {params} -w {threads} "
              "-i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} "
              "-j {output.json} -h {output.html} 2> {log}"
            )

##======== Step 2.1 Build STAR 1-Index ========
rule StarIndex1:
    input:
        dna = DNA,
        gtf = GTF
    output:
        WORKDIR + "Step2.StarAlign/star-1-index/sjdbList.out.tab"
    log:
        WORKDIR + "logs/Step2.StarAlign/build.star-1-index.log"
    benchmark:
        WORKDIR + "benchmark/Step2.StarAlign/star-1-index.benchmark.tsv"
    threads:
        16
    resources:
        mem = 40000
    params:
        "--runMode genomeGenerate --sjdbOverhang 149"
    run:
        GenomeDir = WORKDIR + "Step2.StarAlign/star-1-index/"

        shell("ln -sf {input.dna} {GenomeDir}/genome.fa ")
        shell("STAR {params} --runThreadN {threads} --sjdbGTFfile {input.gtf} "
              "--genomeDir {GenomeDir} --genomeFastaFiles {GenomeDir}/genome.fa "
              "2>{log} "
            )
 
## ======== Step 2.2 STAR aliament 1-pass ========
rule StarPass1:
    input:
        R1 = WORKDIR + "Step1.FastqFilter/{sample}.R1.fq.gz",
        R2 = WORKDIR + "Step1.FastqFilter/{sample}.R2.fq.gz"
    output:
        WORKDIR + "Step2.StarAlign/star-1-pass/{sample}SJ.out.tab"
    log:
        WORKDIR + "logs/Step2.StarAlign/{sample}.star-1-pass.log"
    benchmark:
        WORKDIR + "benchmark/Step2.StarAlign/{sample}.star-1-pass.benchmark.tsv"
    threads:
        8
    resources:
        mem = 20000
    params:
        opt1 = "--runMode alignReads --readFilesCommand zcat ",
        opt2 = "--outSAMmapqUnique 255 --outSAMtype BAM SortedByCoordinate ",
        opt3 = "--outFilterMismatchNmax 10 --outFilterMultimapNmax 20 ",
        opt4 = "--alignIntronMin 21 --alignIntronMax 0 ",
        opt5 = "--alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignMatesGapMax 0 ",
        pfx = WORKDIR + "Step2.StarAlign/star-1-pass/{sample}"
    run:
        GenomeDir = WORKDIR + "Step2.StarAlign/star-1-index/"

        shell(
              "STAR {params.opt1} {params.opt2} {params.opt3} {params.opt4} {params.opt5} "
              "--runThreadN {threads} --genomeDir {GenomeDir} --outFileNamePrefix {params.pfx} "
              "--sjdbGTFfile {GTF} --outBAMsortingThreadN {threads} "
              "--readFilesIn {input.R1} {input.R2} 2>{log} "
            )

## ======== Step 2.3 merge splice junctions ======== 
rule MergeSpliceJunctions:
    input:
        sjs = expand(WORKDIR + "Step2.StarAlign/star-1-pass/{sample}SJ.out.tab", sample = SAMPLES)
    output:
        tab = WORKDIR + "Step2.StarAlign/SJ.Pass-1-Merged.tab"
    log:
        WORKDIR + "logs/Step2.StarAlign/MergeSpliceJunctions.log"
    benchmark:
        WORKDIR + "benchmark/Step2.StarAlign/MergeSpliceJunctions.benchmark.tsv"
    run:
        # Retain splice junctions with at least 3 uniquely mapped fragments per sample.
        shell("cat {input.sjs} | awk '$7 >= 3' | cut -f1-4 | sort -u > {output.tab}")

## ======== Step 2.4 Star 2 index ========
rule StarIndex2:
    input:
        tab = WORKDIR + "Step2.StarAlign/SJ.Pass-1-Merged.tab"
    output:
        WORKDIR + "Step2.StarAlign/star-2-index/sjdbList.out.tab"
    params:
        "--runMode genomeGenerate --sjdbOverhang 149"
    threads:
        16
    log:
        WORKDIR + "logs/Step2.StarAlign/star-2-index.log"
    run:
        GenomeDir = WORKDIR + "Step2.StarAlign/star-2-index/"

        shell(
            "STAR {params} --runThreadN {threads} --genomeDir {GenomeDir} "
            "--genomeFastaFiles {DNA} --sjdbFileChrStartEnd {input.tab} "
            "--sjdbGTFfile {GTF} 2>{log} "
            )

## ======= Step 2.5: Star 2 pass Alignment ========
rule StarPass2:
    input:
        R1 = WORKDIR + "Step1.FastqFilter/{sample}.R1.fq.gz",
        R2 = WORKDIR + "Step1.FastqFilter/{sample}.R2.fq.gz",
        tab = WORKDIR + "Step2.StarAlign/SJ.Pass-1-Merged.tab"
    output:
        bam = WORKDIR + "Step2.StarAlign/star-2-pass/{sample}Aligned.sortedByCoord.out.bam",
        sjs = WORKDIR + "Step2.StarAlign/star-2-pass/{sample}SJ.out.tab"
    log:
        WORKDIR + "logs/Step2.StarAlign/{sample}.star-2-pass.log"
    benchmark:
        WORKDIR + "benchmark/Step2.StarAlign/{sample}.star-2-pass.tsv"
    threads:
        8
    resources:
        mem = 20000
    params:
        opt1 = "--runMode alignReads --readFilesCommand zcat ",
        opt2 = "--outSAMmapqUnique 255 --outSAMtype BAM SortedByCoordinate ",
        opt3 = "--outFilterMismatchNmax 10 --outFilterMultimapNmax 20 ",
        opt4 = "--alignIntronMin 21 --alignIntronMax 0 ",
        opt5 = "--alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignMatesGapMax 0 ",
        pfx = WORKDIR + "Step2.StarAlign/star-2-pass/{sample}"
    run:
        GenomeDir = WORKDIR + "Step2.StarAlign/star-2-index/"

        shell(
            "STAR {params.opt1} {params.opt2} {params.opt3} {params.opt4} {params.opt5} "
              "--runThreadN {threads} --genomeDir {GenomeDir} --outFileNamePrefix {params.pfx} "
              "--sjdbGTFfile {GTF} --outBAMsortingThreadN {threads} "
              "--sjdbFileChrStartEnd {input.tab} --readFilesIn {input.R1} {input.R2} 2>{log} "
            )
            