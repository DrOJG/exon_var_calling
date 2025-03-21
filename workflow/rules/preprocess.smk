import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table("config/samples.tsv").set_index("id", drop=False)

def build_fastq_dicts(samples):
    """
    Takes in samples df, builds a sample id in "{sample}_{exon}" format,
   and returns two dicts (fq1_dict, fq2_dict) in form {sample_id: fastq_file}
    """
    samples["sample_name"] = samples[["sample","exon"]].apply("_".join, axis=1)

    fq1_dict = dict(zip(samples.sample_name, samples.fq1))
    fq2_dict = dict(zip(samples.sample_name, samples.fq2))

    return(fq1_dict, fq2_dict)

fq1_dict, fq2_dict = build_fastq_dicts(samples)

rule fastp_trim_and_filter:
    input:
        fq1=lambda wildcards: "data/" + fq1_dict[wildcards.sample_name],
        fq2=lambda wildcards: "data/" + fq2_dict[wildcards.sample_name],
    
    output:
        fq1_trim=temp("results/trimmed_fastq/{sample_name}_R1_trimmed.fastq.gz"),
        fq2_trim=temp("results/trimmed_fastq/{sample_name}_R2_trimmed.fastq.gz"),
        json="results/logs/fastp/{sample_name}_fastp.json",
        html="results/logs/fastp/{sample_name}_fastp.html"
    
    threads: config["threads"]
    
    conda:
        "ruleenvs/fastp.yml"
    log:
        "results/logs/fastp_logs/{sample_name}_fastp.log"
    shell:
        """
        fastp --thread {threads} -i {input.fq1} \
        -I {input.fq2} \
        -o {output.fq1_trim} \
        -O {output.fq2_trim} \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --json {output.json} \
        --html {output.html} \
        --trim_poly_g > {log} 2>&1
        """

rule bwa_align_and_sort:
    input:
        fq1_trim="results/trimmed_fastq/{sample_name}_R1_trimmed.fastq.gz",
        fq2_trim="results/trimmed_fastq/{sample_name}_R2_trimmed.fastq.gz",
        ref=config["reference"]
    output:
        temp("results/aligned_bams/{sample_name}_sorted.bam")

    params:
        rg=r"@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA"

    threads: config["threads"]
    
    conda:
        "ruleenvs/bwa.yml"
    shell:
        """
        bwa mem -t {threads} -R "{params.rg}" "{input.ref}" "{input.fq1_trim}" "{input.fq2_trim}" |
        samtools sort -@ {threads} -m 1G -o {output} - 
        """