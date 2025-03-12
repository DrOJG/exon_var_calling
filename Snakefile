import pandas as pd

include: "workflow/rules/preprocess.smk"
include: "workflow/rules/bqsr.smk"

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule all:
    input:
        expand("results/final_bams/{id}_sorted_bqsr.bam", id=samples["id"])