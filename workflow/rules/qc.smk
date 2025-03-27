import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule samtools_flagstat:
    input:
        "results/final_bams/{sample_name}_sorted_bqsr.bam",
    output:
        "results/flagstat/{sample_name}_sorted_bqsr.bam.flagstat",
    log:
        "results/logs/flagstat/{sample_name}.log",
    wrapper:
        "v5.9.0/bio/samtools/flagstat"

rule multiqc_all:
    input:
        expand("results/logs/fastp/{sample}_{exon}_fastp.json",
                sample=samples["sample"],
                exon=samples["exon"]),
        expand("results/flagstat/{sample}_{exon}_sorted_bqsr.bam.flagstat",
                sample=samples["sample"],
                exon=samples["exon"]),
        expand("results/snpeff_stats/{sample}_{caller}_snpeff.csv",
                sample=samples["sample"],
                caller=["hapcaller", "freebayes", "lofreq", "mutect2"]),
    output:
        "results/multiqc/multiqc_final.html",
        directory("results/multiqc/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "results/logs/multiqc.log",
    wrapper:
        "v5.9.0/bio/multiqc"