import pandas as pd

include: "workflow/rules/preprocess.smk"
include: "workflow/rules/bqsr.smk"
include: "workflow/rules/varcall.smk"
include: "workflow/rules/vcf_process.smk"

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule all:
    input:
        expand("results/final_bams/{sample}_{exon}_sorted_bqsr.bam",
                sample=samples["sample"],
                exon=samples["exon"]),
        expand("results/vcf/final/{sample}_{caller}_merged_filtered_snpeff.vcf.gz",
                sample=samples["sample"],
                caller=["hapcaller", "freebayes"])