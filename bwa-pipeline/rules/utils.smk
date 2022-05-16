import sys
import os
from pathlib import Path
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
import os

min_version("5.9.1")

units = pd.read_csv("input.tsv", sep = "\t", dtype=str).set_index(["patient", "samp", "readgroup"], drop=False)
units = units.sort_index()
units['key'] = units[['patient', 'samp', 'readgroup']].agg('.'.join, axis=1)
keys = units.key.unique().tolist()
smps = units.samp.unique().tolist()
rgs = units.readgroup.unique().tolist()
pats = units.patient.unique().tolist()

##  Make patient/sample/readgroup path
paths = [key.replace(".", "/") for key in keys]
for path in paths:
    os.makedirs(path, exist_ok=True)

idx = "/data/personalized_medicine/ref/GRCh38.d1.vd1.fa"
dbsnp = "/data/databases/dbSNP/UCSC_dbSNP153Common_hg38_combined.chr.vcf.gz"

def get_fastq(wildcards):
    return {'fastq1' : units.loc[(wildcards.patient, wildcards.sample, wildcards.readgroup), 'f1'], 
            'fastq2' : units.loc[(wildcards.patient, wildcards.sample, wildcards.readgroup), 'f2']}

""" def get_bams(wildcards):
    samps = units.loc[(wildcards.patient), 'samp'].unique().tolist()
    rgs = units.loc[(wildcards.patient), 'readgroup'].unique().tolist()
    bams = [f"bams/{wildcards.patient}.{samp}.{rg}.merged.bam" for samp,rg in zip(samps,rg)]
    return bams
"""