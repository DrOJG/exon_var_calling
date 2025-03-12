import pandas as pd

include: "workflow/rules/preprocess.smk"

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule all:
    input:
        expand("results/aligned_bams/{id}_sorted.bam", id=samples["id"])
