import sys
import os
import calendar
import time
import pathlib
from pathlib import Path
import pandas as pd
from snakemake.utils import validate, min_version

min_version("5.9.1")
validate(config, "../schemas/config.schema.yml")

units = pd.read_csv(config["input"], sep = "\t", dtype=str).set_index(["patient", "samp", "readgroup"], drop=False)
units = units.sort_index()
units['key'] = units[['patient', 'samp', 'readgroup']].agg('.'.join, axis=1)
units['f1'] = units['f1'].apply(lambda x: pathlib.Path(x).absolute())
units['f2'] = units['f2'].apply(lambda x: pathlib.Path(x).absolute())

out = pathlib.Path(config["outdir"]).absolute()
# tmp = pathlib.Path(config["tempdir"]).absolute()

## change wd
# ts = str(calendar.timegm(time.gmtime()))
# tmpdir = f"{tmp}/{ts}"
# pathlib.Path(tmpdir).mkdir(parents=True, exist_ok=True)
# pathlib.Path(f"{tmpdir}/.logs").mkdir(parents=True, exist_ok=True)
#out.mkdir(parents=True, exist_ok=True)

"""
The cleanup step in the end of the pipeline will remove any intermidiate file. If old files are in the outdir, they're considered 
as intermidiate, thus deleted. If considering removing the FileExistsError, please consider to adjust the code to avoid accidental file removals!
"""
if out.exists():
    raise FileExistsError(f"The output directory {out} exists, please use a different name or remove it") 

pathlib.Path(f"{out}/.logs").mkdir(parents=True)
workdir: out

## Get unique lists
keys = units.key.unique().tolist()
smps = units.samp.unique().tolist()
rgs = units.readgroup.unique().tolist()
pats = units.patient.unique().tolist()

##  Make patient/sample/readgroup path
paths = [key.replace(".", "/") for key in keys]
for path in paths:
    os.makedirs(path, exist_ok=True)

idx = config["idx"]
dbsnp = config["dbsnp"] 

def get_fastq(wildcards):
    return {'fastq1' : units.loc[(wildcards.patient, wildcards.sample, wildcards.readgroup), 'f1'], 
            'fastq2' : units.loc[(wildcards.patient, wildcards.sample, wildcards.readgroup), 'f2']}

