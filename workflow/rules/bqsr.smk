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