configfile: "config/config.yaml"

rule haplotypecaller:
    input:
        bam="results/final_bams/{sample_name}_sorted_bqsr.bam",
        idx="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
        ref=config["reference"],
        known=config["known"],
        interval=config["regions"], # optional
    output:
        vcf="results/vcf/haplotypecaller/{sample_name}_haplotypecaller.vcf.gz",
        idx="results/vcf/haplotypecaller/{sample_name}_haplotypecaller.vcf.gz.tbi"
    log:
        "results/logs/gatk/haplotypecaller/{sample_name}.log",
    params:
        extra=lambda wildcards, input: f"-L {input.interval} --max-reads-per-alignment-start 0 --create-output-variant-index",
    threads: config["threads"]
    resources:
        mem_mb=8192,
    wrapper:
        "v5.8.3/bio/gatk/haplotypecaller"

rule freebayes_bed:
    input:
        alns="results/final_bams/{sample_name}_sorted_bqsr.bam",
        idxs="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
        ref=config["reference"],
        regions=config["regions"],
    output:
        vcf = "results/vcf/freebayes/{sample_name}_freebayes.vcf.gz",
    log:
        "results/logs/freebayes/{sample_name}.log",
    params:
    
    threads: config["threads"]
    resources:
        mem_mb=8192,
    wrapper:
        "v5.9.0/bio/freebayes"

rule lofreq:
    input:
        bam="results/final_bams/{sample_name}_sorted_bqsr.bam",
        bai="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
        reference=config["reference"],
        regions=config["regions"],
    output:
        "results/vcf/lofreq/{sample_name}_lofreq.vcf.gz",
    log:
        "results/logs/lofreq_call/{sample_name}_lofreq.log"
    params:
        ref=lambda wildcards, input: f"{input.reference}",
        extra=lambda wildcards, input: f"-l {input.regions}",
    threads: config["threads"]
    wrapper:
        "v5.9.0/bio/lofreq/call"

rule gatk_get_pileup_summaries:
    input:
        bam="results/final_bams/{sample_name}_sorted_bqsr.bam",
        idx="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
        intervals=config["regions"],
        variants=config["germline"],
    output:
        "results/pileups/{sample_name}_pileups.table",
    threads: 2
    resources:
        mem_mb=1024,
    params:
        extra="",
    log:
        "results/logs/gatk/get_pileups/{sample_name}_summary.log",
    wrapper:
        "v6.2.0/bio/gatk/getpileupsummaries"

rule gatk_calculate_contamination:
    input:
        tumor="results/pileups/{sample_name}_pileups.table",
    output:
        "results/contamination/{sample_name}_contamination.table",
    threads: 2
    resources:
        mem_mb=1024,
    log:
        "results/logs/gatk/calculate_contamination/{sample_name}_contamination.log",
    params:
        extra="",
    wrapper:
        "v6.2.0/bio/gatk/calculatecontamination"

rule mutect2:
    input:
        fasta=config["reference"],
        map="results/final_bams/{sample_name}_sorted_bqsr.bam",
        idx="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
        interval=config["regions"],
        pon=config["pon"],
    output:
        vcf="results/vcf/mutect2/{sample_name}_mutect2.vcf.gz",
        idx="results/vcf/mutect2/{sample_name}_mutect2.vcf.gz.tbi",
    threads: config["threads"]
    resources:
        mem_mb=8192,
    params:
        extra=lambda wildcards, input: f"-L {input.interval} --max-reads-per-alignment-start 0 --create-output-variant-index",
    log:
        "results/logs/mutect2/{sample_name}_mutect2.log",
    wrapper:
        "v5.9.0/bio/gatk/mutect"