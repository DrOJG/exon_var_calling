#!/usr/bin/env python3

from pathlib import Path
from pysam import VariantFile 
import pandas as pd

testVCF = Path("/home/oliver/Documents/NGS_analysis/CD7_exon_results/results_TvT10-BM_20250328/vcf/final/TvT10-BM_mutect2_merged_filtered_snpeff.vcf.gz")
#testVCF = Path("/home/oliver/Documents/NGS_analysis/CD7_exon_results/G944_freebayes_merged_filtered_snpeff.vcf.gz")
outDF = pd.DataFrame(columns=["chromosome",
                              "position",
                              "ref",
                              "alt",
                              "var_caller",
                              "reads_ref",
                              "reads_alt",
                              "freq_alt",
                              "percent_alt"])

with VariantFile(testVCF, "rb") as vcf:
    for rcrd in vcf.fetch():
        adList = []
        afList = []
        for sample in rcrd.samples.items():
            adRecs = [x[1] for x in sample[1].items() if x[0] == 'AD' and x[1] != (None,)]
            afRecs = [x[1] for x in sample[1].items() if x[0] == 'AF' and x[1] != (None,)]
            if adRecs:
                adList.extend(adRecs)
            if afRecs:
                # If multiple alleles, keep af as tuple, else convert to float
                if len(rcrd.alts) == 1:
                    afRecs = map(lambda x: float(x[0]), afRecs)
                afList.extend(afRecs)

        if len(adList) > 1:
            sumAdList = list(map(sum, adList))
            exonIdx = sumAdList.index(max(sumAdList))
        else:
            exonIdx = 0
        
        ADTuple = adList[exonIdx]
        AF = afList[exonIdx]
        
        readsRef = ADTuple[0]
        if len(ADTuple) > 2:
            readsAlt = ADTuple[1:(len(ADTuple)-1)]
            pcentAlt = tuple(map(lambda x: x *100, AF))
        else:
            readsAlt = ADTuple[1]
            pcentAlt = AF *100

        rowDict = {"chromosome": rcrd.chrom,
                   "position": rcrd.pos,
                   "ref": rcrd.ref,
                   "alt": rcrd.alts,
                   "var_caller": "mutect2",
                   "reads_ref": readsRef,
                   "reads_alt": readsAlt,
                   "freq_alt": AF,
                   "percent_alt": AF * 100}
        outDF.loc[len(outDF)] = rowDict

print(outDF)