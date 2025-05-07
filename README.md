# **Snakemake workflow for single gene exon variant calling**


### **Inputs**
**This workflow will take the following inputs:** 
* Folder of raw, paired fastq files (R1 and R2 for each sample).  
* Sample table with columns: [id, sample, exon, fq1, fq2]:
  - **id**: a unique identifier for indexing, usual just use numbers 1, 2, 3... etc.
  - **sample**: A sample name, that will be shared between different exons of the same sample. Will be the final name of merged samples and vcfs. No spaces or underscores.
  - **exon**: Exon identifier, e.g. ex1, Exon2. No spaces or underscores.
  - **fq1**: full filename of read 1 FASTQ file, including extension. do not include path.
  - **fq2**: full filename of read 2 FASTQ file, including extension. do not include path.
* Reference genome. Should be indexed with BWA, samtools faidx and Picard CreateSequenceDictionary. Currently script does not do this automatically.
* Common sites VCF compatible with your reference genome.
* A BED file with the target regions of your sequencing. Can either be your specific amplicons or a region encompassing your whole gene of interest
* A tab-delimited table with columns [amplicon, r1_adapter, r2_adapter]:
  - **amplicon**: An amplicon identifier for your exon PCRs. Name not used, just for reference.
  - **r1_adapter**: a linked adapter for cutadapt to remove PCR primers from read 1, in form: ^FWDPRIMER...RCREVPRIMER.
  - **r2_adapter**: a linked adapter for cutadapt to remove PCR primers from read 2, in form: ^REVPRIMER...RCFWDPRIMER.
* config.yaml file with general options plus paths to reference, sample table, primer adapters table etc.  


### **Outputs**
**Final outputs will be:**
* Individual aligned BAM files for every sequenced sample. BAM files are sorted, and base score recalibrated.
* Merged and filtered VCF files for every samples for all exons. Varient calling with: Haplotypecaller, Mutect2, Lofreq and Freebayes. Annotated with SNPEff.
* MultiQC report summarising outputs.
* A summary CSV file for each sample with all variants called by the 4 tools.


### **Pipeline**
**Will follow same schema we have been using for exon analysis Bash scripts:**  
1. Trim and quality filter with fastp.  
2. Remove PCR primers from reads with cutadapt
3. Align each file to reference genome (Single chromosome hg38 ref) with BWA-MEM (Maybe BWA-MEM2 as it is faster).  
4. Carry out base quality score recalibration with GATK baserecalibrator and ApplyBQSR.  
5. Index output BAM files with samtools.
6. Variant calling with Haplotypecaller, Mutect2, Lofreq and Freebayes.
7. Pre-filter Mutect2 VCFs with FilterMutectCalls.
8. Merge all exons from each sample with BCFtools merge.
9. Filter merged VCFs for depth and quality with BCFtools filter.
10. Annotate variants with SNPeff.
11. Summarise VCF files into single CSV for readability (custom script).
11. Get mapping stats for BAM files with samtools flagstat.
12. Summarise QC and metrics with MultiQC.

### **Useage**

* Make sure to update the the config.yaml file with paths to inputs, and to create the sample and primer tables.
* Create a conda env from the provided `snakemake.yml` file with:
  ```
  conda env create -f ./workflow/envs/snakemake.yml
  ```

* For local use, the workflow can be run with:
  ```
  snakemake --sdm conda
  ```
* A runscript is provided which can be used for job submission on SGE computer clusters. Note that the workflow currently only runs on a single cluster node.



