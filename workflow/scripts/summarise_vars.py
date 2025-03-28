#!/usr/bin/env python3

from pathlib import Path
from pysam import VariantFile 
import pandas as pd

testVCF = Path("/home/oliver/Documents/NGS_analysis/CD7_exon_results/results_TvT10-BM_20250328/vcf/final/TvT10-BM_mutect2_merged_filtered_snpeff.vcf.gz")

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
                afList.extend(afRecs)
        print(afList)
        adSum = [sum(x) for x in zip(*adList)]
        calcAF = adSum[1] / sum(adSum)
        print(calcAF)
        #print(adSum)
        # rowDict = {"chromosome": rcrd.chrom,
        #            "position": rcrd.pos,
        #            "ref": rcrd.ref,
        #            "alt": list(rcrd.alts),
        #            "var_caller": "mutect2",
        #            "reads_ref": rcrd.format["AD"][0],
        #            "reads_alt": rcrd.format["AD"][1],
        #            "freq_alt": rcrd.format["AF"],
        #            "percent_alt": rcrd.format["AF"] * 100}
        # outDF.loc[len(outDF)] = rowDict

#print(outDF)