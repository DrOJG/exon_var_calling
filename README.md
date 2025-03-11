# **Snakemake workflow for single gene exon variant calling**
## WIP
This workflow will take the following inputs:  
* Folder of raw, paired fastq files (R1 and R2 for each sample).  
* Sample table with columns: [id, sample, exon].  
* Reference genome (We can make script check for and download this if needed).  
* config.yml file with general options plus links to ref and samp table.  

Final output will be:  
* Individual aligned BAM files for every sequenced sample.  
* Merged and filtered VCF files for every samples for all exons.  
* Plot of variants on gene map? Check if possible.  

Will follow same schema we have been using for exon analysis Bash scripts:  
1. Trim and quality filter with fastp.  
2. multiqc of fastp outputs.  
3. Align each file to reference genome (Single chromosome hg38 ref) with BWA-MEM (Maybe BWA-MEM2 as it is faster).  
4. Carry out base quality score recalibration with GATK baserecalibrator and ApplyBQSR.  
5. Index output BAM files with samtools.  
6. Samtools flagstat and multiqc on flagstat outputs.  
7. Variant calling with ? (need to decide best tool).  
8. Post-filter VCFs with bcftools (or variant-caller-specific tool).  
9. Merge per-sample exon VCFs with bcftools merge.  
10. (?) Plot metrics and variants.  


