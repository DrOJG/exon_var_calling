configfile: "config/config.yaml"

rule haplotype_caller:
    input:
        bam="results/final_bams/{sample_name}_sorted_bqsr.bam",
        ref=config["reference"],
        known=config["known"]  # optional
    output:
        vcf="results/vcf/haplotypecaller/{sample}.vcf",
    log:
        "results/logs/gatk/haplotypecaller/{sample}.log",
    params:
        extra="-L chr17 --max-reads-per-alignment-start 0",
    threads: config["threads"]
    resources:
        mem_mb=1024,
    wrapper:
        "v5.8.3/bio/gatk/haplotypecaller"