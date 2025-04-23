import pandas as pd

include: "workflow/rules/preprocess.smk"
include: "workflow/rules/bqsr.smk"
include: "workflow/rules/varcall.smk"
include: "workflow/rules/vcf_process.smk"
include: "workflow/rules/qc.smk"
configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule all:
    input:
        expand("results/final_bams/{sample}_{exon}_sorted_bqsr.bam",
                sample=samples["sample"],
                exon=samples["exon"]),
        expand("results/vcf/final/{sample}_{caller}_merged_filtered_snpeff.vcf.gz",
                sample=samples["sample"],
                caller=["hapcaller", "freebayes", "lofreq", "mutect2"]),
        expand("results/vcf/final/{sample}_{caller}_merged_filtered_snpeff.vcf.gz.tbi",
                sample=samples["sample"],
                caller=["hapcaller", "freebayes", "lofreq", "mutect2"]),
        "results/multiqc/multiqc_final.html",
        expand("results/vcf/final/summaries/{sample}_variant_summary.csv",
                sample=samples["sample"]),