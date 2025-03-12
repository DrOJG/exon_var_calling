import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("id", drop=False)

rule fastp_trim_and_filter:
    input:
        fq1="data/{id}_R1_001.fastq.gz",
        fq2="data/{id}_R2_001.fastq.gz"
    
    output:
        fq1_trim="results/trimmed_fastq/{id}_R1_001_trimmed.fastq.gz",
        fq2_trim="results/trimmed_fastq/{id}_R2_001_trimmed.fastq.gz",
        json="results/logs/fastp/{id}_fastp.json",
        html="results/logs/fastp/{id}_fastp.html"
    
    threads: config["threads"]
    
    conda:
        "/home/oliver/Python_projects/snakemake_workflows/exon_var_calling/workflow/envs/fastp.yml"
    
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
        --trim_poly_g
        """

rule bwa_align_and_sort:
    input:
        fq1_trim="results/trimmed_fastq/{id}_R1_001_trimmed.fastq.gz",
        fq2_trim="results/trimmed_fastq/{id}_R2_001_trimmed.fastq.gz",
        ref=config["reference"]
    output:
        "results/aligned_bams/{id}_sorted.sam"

    params:
        rg=r"@RG\tID:{id}\tSM:{id}"

    threads: config["threads"]
    
    conda:
        "/home/oliver/Python_projects/snakemake_workflows/exon_var_calling/workflow/envs/bwa.yml"
    shell:
        "bwa mem -t {threads} -R '{params.rg}' '{input.ref}' '{input.fq1_trim}' '{input.fq2_trim}' > {output}"
       # """
       # bwa mem -t {threads} -R "{params.rg}" "{input.ref}" "{input.fq1_trim}" "{input.fq2_trim}" |
       #  samtools sort -@ {threads} -m 1G -o {output} -
      #  """