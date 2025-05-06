import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

def build_fastq_dicts(samples):
    """
    Takes in samples df, builds a sample id in "{sample}_{exon}" format,
   and returns two dicts (fq1_dict, fq2_dict) in form {sample_id: fastq_file}
    """
    samples["sample_name"] = samples[["sample","exon"]].apply("_".join, axis=1)

    fq1_dict = dict(zip(samples.sample_name, samples.fq1))
    fq2_dict = dict(zip(samples.sample_name, samples.fq2))

    return(fq1_dict, fq2_dict)

def get_cutadapt_flags(tsv_file):
    """
    Reads in the primer adapters file and prepares them into adapter strings for cutadapt
    """
    r1_flags = []
    r2_flags = []
    with open(tsv_file) as f:
        next(f)  # skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                r1_flags.append(f"-a {fields[1]}")
                r2_flags.append(f"-A {fields[2]}")
    return ' '.join(r1_flags), ' '.join(r2_flags)

fq1_dict, fq2_dict = build_fastq_dicts(samples)

r1_flags_str, r2_flags_str = get_cutadapt_flags(config["primers"])


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

rule cutadapt_remove_primers:
    input: 
        fq1_trim="results/trimmed_fastq/{sample_name}_R1_trimmed.fastq.gz",
        fq2_trim="results/trimmed_fastq/{sample_name}_R2_trimmed.fastq.gz",
    
    output:
        fq1_trim_ca=temp("results/trimmed_fastq/{sample_name}_R1_trimmed_cutadapt.fastq.gz"),
        fq2_trim_ca=temp("results/trimmed_fastq/{sample_name}_R2_trimmed_cutadapt.fastq.gz"),

    params:
        r1_primers=r1_flags_str,
        r2_primers=r2_flags_str,
    
    log:
        "results/logs/cutadapt/{sample_name}_cutadapt.log"
    
    conda:
        "ruleenvs/cutadapt.yml"
    
    threads: config["threads"]
    
    shell:
        """
        cutadapt {params.r1_primers} \
        {params.r2_primers} \
        --discard-untrimmed \
        --times 1 \
        --action=trim \
        -m 150 \
        -o {output.fq1_trim_ca} \
        -p {output.fq2_trim_ca} \
        {input.fq1_trim} {input.fq1_trim} > {log} 2>&1
        """


rule bwa_align_and_sort:
    input:
        fq1_trim_ca="results/trimmed_fastq/{sample_name}_R1_trimmed_cutadapt.fastq.gz",
        fq2_trim_ca="results/trimmed_fastq/{sample_name}_R2_trimmed_cutadapt.fastq.gz",
        ref=config["reference"]
    
    output:
        temp("results/aligned_bams/{sample_name}_sorted.bam")

    params:
        rg=r"@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA"

    threads: config["threads"]
    
    log:
        "results/logs/bwa/{sample_name}_bwa.log"
    conda:
        "ruleenvs/bwa.yml"
    shell:
        """
        bwa mem -t {threads} -R "{params.rg}" "{input.ref}" "{input.fq1_trim_ca}" "{input.fq2_trim_ca}" |
        samtools sort -@ {threads} -m 1G -o {output} - > {log} 2>&1
        """