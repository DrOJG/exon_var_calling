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
        mem_mb=1024,
    wrapper:
        "v5.8.3/bio/gatk/haplotypecaller"

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
