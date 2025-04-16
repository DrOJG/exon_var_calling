#!/usr/bin/env python3
import re
from pathlib import Path
from pysam import VariantFile 
import pandas as pd

class VCFSummaryBuilder:
    def __init__(self, VariantRecord, inputPath):
        self.rcrd = VariantRecord
        self.filename = inputPath.stem
        self.varCaller = self.assign_var_caller()
        self.ADTuple, self.exonIdx = self.get_allele_depths()
        self.AF = self.get_allele_frequencies()
        self.altAlleles = self.format_alt_alleles()
        self.readsRef, self.readsAlt, self.pcentAlt = self.calculate_percent_frequency()

    def assign_var_caller(self):
        if re.search('mutect', self.filename) is not None:
            return "Mutect2"
        elif re.search('hapcaller', self.filename) is not None:
            return "HaplotypeCaller"
        elif re.search('lofreq', self.filename) is not None:
            return "Lofreq"
        elif re.search('freebayes', self.filename) is not None:
            return "Freebayes"
        else:
            raise ValueError("Unrecognised variant caller format in input")
    
    def get_allele_depths(self):
        # Lofreq has DP4 list in info which has reads as (ref_fwd, ref_rev, alt_fwd, alt_rev)
        if self.varCaller == "Lofreq":
            dp4 = self.rcrd.info.get("DP4")
            ADTuple = (dp4[0] + dp4[1], dp4[2] + dp4[3])
            exonIdx = None
            return(ADTuple, exonIdx)
        # All other var callers, get allele depths from per-sample format
        else:   
            adList = []
            for sample in self.rcrd.samples.items():
                adRecs = [x[1] for x in sample[1].items() if x[0] == 'AD' and x[1] != (None,)]
                if adRecs:
                    adList.extend(adRecs)
            # In cases where >1 exons have reads, find the main exon by its higher read depth
            if len(adList) > 1:
                filterAdList = []
                for tuple in adList:
                    filteredTuple = [x for x in tuple if x is not None]
                    filterAdList.append(filteredTuple)
                sumAdList = list(map(sum, filterAdList))
                exonIdx = sumAdList.index(max(sumAdList))
            else:
                exonIdx = 0
            
            ADTuple = adList[exonIdx]
            return(ADTuple, exonIdx)
    
    def get_allele_frequencies(self):
        if self.varCaller == "Mutect2":
            afList = []
            for sample in self.rcrd.samples.items():
                afRecs = [x[1] for x in sample[1].items() if x[0] == 'AF' and x[1] != (None,)]
                if afRecs:
                    # If multiple alleles, keep af as tuple, else convert to float
                    if len(self.rcrd.alts) == 1:
                        afRecs = map(lambda x: float(x[0]), afRecs)
                    afList.extend(afRecs)
            AF = afList[self.exonIdx]
            return AF
        else:
            AF = self.rcrd.info.get("AF")
            if isinstance(AF, tuple):
                if len(AF) == 1:
                    AF = AF[0]
            return AF

    def format_alt_alleles(self):
        if len(self.rcrd.alts) == 1:
            altAlleles = self.rcrd.alts[0]
        else:
            altAlleles = self.rcrd.alts
        return altAlleles
    
    def calculate_percent_frequency(self):
        readsRef = self.ADTuple[0]
        if len(self.ADTuple) > 2:
            readsAlt = self.ADTuple[1:(len(self.ADTuple)-1)]
            pcentAlt = tuple(map(lambda x: x * 100, self.AF))
        else:
            readsAlt = self.ADTuple[1]
            pcentAlt = self.AF * 100
        return(readsRef, readsAlt, pcentAlt)
    
    def build_output_dict(self):
        
        joiner = lambda t: "|".join(str(x) for x in t)
        
        if isinstance(self.altAlleles, tuple):
            altAllelesFormatted = joiner(self.altAlleles)
            readsAltFormatted = joiner(self.readsAlt)
            AFFormatted = joiner(self.AF)
            pcentAltFormatted = joiner(self.pcentAlt)
        else:
            altAllelesFormatted = self.altAlleles
            readsAltFormatted = self.readsAlt
            AFFormatted = self.AF
            pcentAltFormatted = self.pcentAlt

        rowDict = {"chromosome": self.rcrd.chrom,
                   "position": self.rcrd.pos,
                   "ref": self.rcrd.ref,
                   "alt": altAllelesFormatted,
                   "var_caller": self.varCaller,
                   "reads_ref": self.readsRef,
                   "reads_alt": readsAltFormatted,
                   "freq_alt": AFFormatted,
                   "percent_alt": pcentAltFormatted}
        return rowDict
    
if __name__ == "__main__":

    testVCF1 = Path("/home/oliver/Documents/NGS_analysis/CD7_exon_results/results_TvT10-BM_20250328/vcf/final/TvT10-BM_mutect2_merged_filtered_snpeff.vcf.gz")
    testVCF2 = Path("/home/oliver/Documents/NGS_analysis/CD7_exon_results/G944_freebayes_merged_filtered_snpeff.vcf.gz")
    testVCF3 = Path("/home/oliver/Documents/NGS_analysis/CD7_exon_results/results_TvT10-BM_20250328/vcf/final/TvT10-BM_lofreq_merged_filtered_snpeff.vcf.gz")
    testVCF4 = Path("/home/oliver/Documents/NGS_analysis/CD7_exon_results/results_TvT10-BM_20250328/vcf/final/TvT10-BM_hapcaller_merged_filtered_snpeff.vcf.gz")
    vcfList = [testVCF1, testVCF2, testVCF3, testVCF4]
    
    outDF = pd.DataFrame(columns=["chromosome",
                                "position",
                                "ref",
                                "alt",
                                "var_caller",
                                "reads_ref",
                                "reads_alt",
                                "freq_alt",
                                "percent_alt"])
    
    for input in vcfList:
        with VariantFile(input, "rb") as vcf:
            for rcrd in vcf.fetch():
                varSummary = VCFSummaryBuilder(rcrd, input)
                outDF.loc[len(outDF)] = varSummary.build_output_dict()

    outDF.sort_values(["chromosome", "position"],
                      ascending=[True, True],
                      inplace=True)
    print(outDF)