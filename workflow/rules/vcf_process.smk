import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule bcftools_index_freebayes:
    input:
        "results/vcf/freebayes/{sample_name}_freebayes.vcf.gz",
    output:
        "results/vcf/freebayes/{sample_name}_freebayes.vcf.gz.tbi",
    log:
        "results/logs/bcftools_index/{sample_name}_freebayes.log",
    wrapper:
        "v5.9.0/bio/bcftools/index"

rule bcftools_index_lofreq:
    input:
        "results/vcf/lofreq/{sample_name}_lofreq.vcf.gz",
    output:
        "results/vcf/lofreq/{sample_name}_lofreq.vcf.gz.tbi",
    log:
        "results/logs/bcftools_index/{sample_name}_lofreq.log",
    wrapper:
        "v5.9.0/bio/bcftools/index"

rule gatk_filtermutectcalls:
    input:
        vcf="results/vcf/mutect2/{sample_name}_mutect2.vcf.gz",
        ref=config["reference"],
    output:
        vcf=temp("results/vcf/mutect2_filter/{sample_name}_mutect2_prefilter.vcf.gz"),
        idx=temp("results/vcf/mutect2_filter/{sample_name}_mutect2_prefilter.vcf.gz.tbi"),
    log:
        "results/logs/gatk/filtermutect/{sample_name}_mutect_filter.log",
    params:
        extra="--max-alt-allele-count 4 --max-events-in-region 10 --create-output-variant-index",
    resources:
        mem_mb=4096,
    wrapper:
        "v5.9.0/bio/gatk/filtermutectcalls"

rule bcftools_merge_haplotypecaller:
    input:
        calls=expand("results/vcf/haplotypecaller/{sample}_{item.exon}_hapcaller.vcf.gz",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
        idx=expand("results/vcf/haplotypecaller/{sample}_{item.exon}_hapcaller.vcf.gz.tbi",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
    output:
        temp("results/vcf/merged/{sample}_hapcaller_merged.vcf.gz"),
        temp("results/vcf/merged/{sample}_hapcaller_merged.vcf.gz.tbi"),
    log:
        "results/logs/merge_hapcaller/{sample}_hapcaller_merged.log",
    params:
        uncompressed_bcf=False,
        extra="--merge both --write-index=tbi",
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
        temp("results/vcf/merged/{sample}_freebayes_merged.vcf.gz"),
        temp("results/vcf/merged/{sample}_freebayes_merged.vcf.gz.tbi"),
    log:
        "results/logs/merge_freebayes/{sample}_freebayes_merged.log",
    params:
        uncompressed_bcf=False,
        extra="--merge both --write-index=tbi",
    wrapper:
        "v5.9.0/bio/bcftools/merge"

rule bcftools_merge_lofreq:
    input:
        calls=expand("results/vcf/lofreq/{sample}_{item.exon}_lofreq.vcf.gz",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
        idx=expand("results/vcf/lofreq/{sample}_{item.exon}_lofreq.vcf.gz.tbi",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
    output:
        temp("results/vcf/merged/{sample}_lofreq_merged.vcf.gz"),
        temp("results/vcf/merged/{sample}_lofreq_merged.vcf.gz.tbi"),
    log:
        "results/logs/merge_lofreq/{sample}_lofreq_merged.log",
    params:
        uncompressed_bcf=False,
        extra="--merge both --write-index=tbi",
    wrapper:
        "v5.9.0/bio/bcftools/merge"

rule bcftools_merge_mutect:
    input:
        calls=expand("results/vcf/mutect2_filter/{sample}_{item.exon}_mutect2_prefilter.vcf.gz",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
        idx=expand("results/vcf/mutect2_filter/{sample}_{item.exon}_mutect2_prefilter.vcf.gz.tbi",
                        item=lookup(query="sample == '{sample}'", within=samples),
                        allow_missing=True),
    output:
        # Merge mutect2 vcfs directly to filtered as already filtered
        temp("results/vcf/filtered/{sample}_mutect2_merged_filtered.vcf.gz"),
        temp("results/vcf/filtered/{sample}_mutect2_merged_filtered.vcf.gz.tbi"),
    log:
        "results/logs/merge_mutect/{sample}_mutect2_merged.log",
    params:
        uncompressed_bcf=False,
        extra="--merge both --write-index=tbi",
    wrapper:
        "v5.9.0/bio/bcftools/merge"

rule filter_merged_vcfs:
    input:
        vcf="results/vcf/merged/{sample}_{caller}_merged.vcf.gz",
        index="results/vcf/merged/{sample}_{caller}_merged.vcf.gz.tbi",
        regions=config["regions"],
    output:
        vcf=temp("results/vcf/filtered/{sample}_{caller}_merged_filtered.vcf.gz"),
    log:
        "results/logs/bcftools_filter/{sample}_{caller}_merged_filtered.vcf.gz.log",
    params:
        filter="-i 'INFO/DP>=10 && QUAL>=20'",
    wrapper:
        "v5.9.0/bio/bcftools/filter"

rule bcftools_index_filtered:
    input:
        "results/vcf/filtered/{sample}_{caller}_merged_filtered.vcf.gz",
    output:
        temp("results/vcf/filtered/{sample}_{caller}_merged_filtered.vcf.gz.tbi"),
    log:
        "results/logs/bcftools_index/{sample}_{caller}_merged_filtered.log",
    wrapper:
        "v5.9.0/bio/bcftools/index"

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
        "results/logs/snpeff/{sample}_{caller}.log"
    resources:
        java_opts=f"-XX:ParallelGCThreads={config["threads"]}",
        mem_mb=4096
    wrapper:
        "v5.9.0/bio/snpeff/annotate"

rule bcftools_index_final_vcf:
    input:
        "results/vcf/final/{sample}_{caller}_merged_filtered_snpeff.vcf.gz",
    output:
        "results/vcf/final/{sample}_{caller}_merged_filtered_snpeff.vcf.gz.tbi",
    log:
        "results/logs/bcftools_index/{sample}_{caller}_merged_filtered_snpeff.log",
    wrapper:
        "v5.9.0/bio/bcftools/index"

rule summarise_vcf_files:
    input:
        expand("results/vcf/final/{sample}_{caller}_merged_filtered_snpeff.vcf.gz",
                caller=["hapcaller", "freebayes", "lofreq", "mutect2"],
                allow_missing=True),
    output:
        "results/vcf/final/summaries/{sample}_variant_summary.csv"
    threads: 2
    conda:
        "ruleenvs/pysam.yml"
    script:
        "../scripts/summarise_variants.py"