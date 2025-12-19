#! /usr/bin/env python

import pandas as pd
from modules.utils import getTargetInfo

configfile: "config.yaml"
samplesheet = pd.read_table(config["samples"],dtype=str).set_index("sample_name", drop=False)

def get_fastq(wildcards):
  fastqs = samplesheet.loc[(wildcards.sample), ['fq1', 'fq2']].dropna()
  return {"r1": fastqs.fq1, "r2": fastqs.fq2}


