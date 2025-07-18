configfile: "config/config.yaml"

rule gatk_baserecalibrator:
    input:
        bam="results/aligned_bams/{sample_name}_sorted.bam",
        ref=config["reference"],
        dict=config["refdict"],
        known=config["known"],  # optional known sites - single or a list
    output:
        recal_table="results/recal/{sample_name}.grp",
    log:
        "results/logs/gatk/baserecalibrator/{sample_name}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v5.8.3/bio/gatk/baserecalibrator"


rule gatk_applybqsr:
    input:
        bam="results/aligned_bams/{sample_name}_sorted.bam",
        ref=config["reference"],
        dict=config["refdict"],
        recal_table="results/recal/{sample_name}.grp",
    output:
        bam="results/final_bams/{sample_name}_sorted_bqsr.bam",
    log:
        "results/logs/gatk/gatk_applybqsr/{sample_name}.log",
    params:
        extra="--add-output-sam-program-record",  # optional
        java_opts="",  # optional
        embed_ref=True,  # embed the reference in cram output
    resources:
        mem_mb=1024,
    wrapper:
        "v5.8.3/bio/gatk/applybqsr"


rule samtools_index_final_bam:
    input:
        "results/final_bams/{sample_name}_sorted_bqsr.bam",
    output:
        "results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
    log:
        "results/logs/samtools_index_final_bam/{sample_name}.log",
    params:
        extra="",  # optional params string
    threads: config["threads"]  # This value - 1 will be sent to -@
    wrapper:
        "v5.9.0/bio/samtools/index"