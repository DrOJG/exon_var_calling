import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule bcftools_index_freebayes:
    input:
        "results/vcf/freebayes/{sample_name}_freebayes.vcf.gz",
    output:
        "results/vcf/freebayes/{sample_name}_freebayes.vcf.gz",
    log:
        "results/bcftools_index/{sample_name}.log",
    wrapper:
        "v5.9.0/bio/bcftools/index"

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
        "results/logs/merge_hapcaller/{sample}_hapcaller_merged.log",
    params:
        uncompressed_bcf=False,
        extra="--merge both --write-index",
    wrapper:
        "v5.9.0/bio/bcftools/merge"

rule bcftools_merge_freebayes:
    input:
        calls=expand("results/vcf/freebayes/{sample}_{item.exon}_freebayes.vcf.gz",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
        idx=expand("results/vcf/freebayes/{sample}_{item.exon}_freebayes.vcf.gz.tbi",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
    output:
        "results/vcf/merged/{sample}_freebayes_merged.vcf.gz",
    log:
        "results/logs/merge_freebayes/{sample}_freebayes_merged.log",
    params:
        uncompressed_bcf=False,
        extra="--merge both --write-index",
    wrapper:
        "v5.9.0/bio/bcftools/merge"

rule filter_merged_vcfs:
    input:
        "results/vcf/merged/{sample}_{caller}_merged.vcf.gz",
    output:
        "results/vcf/filtered/{sample}_{caller}_merged_filtered.vcf.gz",
    log:
        "results/log/bcftools_filter/{sample}_{caller}_merged_filtered.vcf.gz.log",
    params:
        filter="-i 'FORMAT/DP>=10 && QUAL>=20'",
        extra=f"-R {config["regions"]}, --write-index",
    wrapper:
        "v5.9.0/bio/bcftools/filter"

rule snpeff_download:
    output:
        directory("snpeff/hg38")
    log:
        "results/logs/snpeff/download/hg38.log"
    params:
        reference="hg38"
    resources:
        mem_mb=1024
    wrapper:
        "v5.9.0/bio/snpeff/download"

rule snpeff:
    input:
        calls="results/vcf/filtered/{sample}_{caller}_merged_filtered.vcf.gz", # (vcf, bcf, or vcf.gz)
        db="snpeff/hg38" # path to reference db downloaded with the snpeff download wrapper
    output:
        calls="results/vcf/final/{sample}_{caller}_merged_filtered_snpeff.vcf.gz",   # annotated calls (vcf, bcf, or vcf.gz)
        stats="results/snpeff_stats/{sample}_{caller}_snpeff.html",  # summary statistics (in HTML), optional
        csvstats="results/snpeff_stats/{sample}_{caller}_snpeff.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/{sample}_{caller}.log"
    resources:
        java_opts=f"-XX:ParallelGCThreads={config["threads"]}",
        mem_mb=4096
    wrapper:
        "v5.9.0/bio/snpeff/annotate"