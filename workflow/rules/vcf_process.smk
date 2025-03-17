import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule bcftools_merge_haplotypecaller:
    input:
        calls=expand("results/vcf/haplotypecaller/{sample}_{item.exon}_hapcaller.vcf.gz",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
        idx=expand("results/vcf/haplotypecaller/{sample}_{item.exon}_hapcaller.vcf.gz.tbi",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
    output:
        "results/vcf/merged/{sample}_hapcaller_merged.vcf.gz",
    log:
        "{sample}_all.log",
    params:
        uncompressed_bcf=False,
        extra="--merge both --write-index",
    wrapper:
        "v5.8.3/bio/bcftools/merge"