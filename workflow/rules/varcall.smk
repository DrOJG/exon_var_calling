configfile: "config/config.yaml"

rule haplotype_caller:
    input:
        bam="results/final_bams/{sample_name}_sorted_bqsr.bam",
        idx="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
        ref=config["reference"],
        known=config["known"] # optional
    output:
        vcf="results/vcf/haplotypecaller/{sample_name}_hapcaller.vcf.gz",
        idx="results/vcf/haplotypecaller/{sample_name}_hapcaller.vcf.gz.tbi"
    log:
        "results/logs/gatk/haplotypecaller/{sample_name}.log",
    params:
        extra="-L chr17 --max-reads-per-alignment-start 0 --create-output-variant-index",
    threads: config["threads"]
    resources:
        mem_mb=8192,
    wrapper:
        "v5.8.3/bio/gatk/haplotypecaller"

rule freebayes_bed:
    input:
        alns="results/final_bams/{sample_name}_sorted_bqsr.bam",
        idxs="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
        ref=config["reference"],
        regions=config["regions"],
    output:
        vcf = "results/vcf/freebayes/{sample_name}_freebayes.vcf.gz",
    log:
        "results/logs/freebayes/{sample_name}.log",
    params:
    
    threads: config["threads"]
    resources:
        mem_mb=8192,
    wrapper:
        "v5.9.0/bio/freebayes"




# Waiting to see if deepvariant gets fixed
# rule deepvariant:
#     input:
#         bam="results/final_bams/{sample_name}_sorted_bqsr.bam",
#         idx="results/final_bams/{sample_name}_sorted_bqsr.bam.bai",
#         ref=config["reference"],
#     output:
#         vcf="results/vcf/deepvariant/{sample_name}_deepvar.vcf.gz",
#     threads: config["threads"]
#     log:
#         directory("results/logs/deepvariant/{sample_name}/")
#     container:
#         "docker://google/deepvariant:1.8.0"
#     shell:
#         """
#         /opt/deepvariant/bin/run_deepvariant \
# 	    --model_type=WGS \
# 	    --ref={input.ref} \
# 	    --reads={input.bam} \
# 	    --regions "chr17" \
# 	    --output_vcf={output.vcf} \
# 	    --num_shards={threads} \
#         --logging_dir={log}
#         """
