#!/usr/bin/env python

import pandas as pd

samples = pd.read_table("samplesheet.tsv").set_index("sample_name", drop=False)
SAMPLES = samples.index.values

#=============================================#

def getTargetInfo(config):
    targetFiles = []
    targetFiles.extend([_getSTARaligns(config)])
    return targetFiles

#=============================================#

def _getSTARaligns(config):
    """ensure that bam and indexes are built"""
    #STAR alignment sorted.bam file, its index, and transcript count
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/STAR/" + sample + "/" + sample + ".sorted.bam")
        ls.append("analysis/" + config["token"] + "/STAR/" + sample + "/" + sample + ".sorted.bam.bai")
    ls.append("analysis/" + config["token"] + "/STAR/STAR_Align_Report.png")   
    return ls
