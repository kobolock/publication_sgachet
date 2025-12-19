#!/usr/bin/env python

import pandas as pd

samples = pd.read_table("samplesheet.tsv").set_index("sample_name", drop=False)
SAMPLES = samples.index.values

#=============================================#

def getTargetInfo(config):
    targetFiles = []
#    targetFiles.extend([_fastqc(config)]),
    targetFiles.extend([
#                        _bbduk(config), # PAS DE TRIMMING POUR LE RNASEQ
#                        _fastqc(config),
                        _getSTARaligns(config),
#                        _rRNAmetrics(config),
#                        _readQC(config),
#                        _get_cuff(config),
#                        _get_cluster(config),
#                        _getDE(config),
#                        _getpathway(config),
#                        _report_html(config),
#                        _getvarscan2mpileup(config),
#                        _preparation_input_files_for_rnaseqcnv(config),
#                        _run_rnaseqcnv(config)
                        ])
    return targetFiles

#=============================================#

## Attention, le trimming n'est pas recommandé pour le RNAseq
#def _bbduk(config):
#    ls = []
#    for sample in SAMPLES:
#        ls.append("analysis/" + config["token"] + "/bbduk/" + sample + "_R1_cleaned.fastq.gz")
#        ls.append("analysis/" + config["token"] + "/bbduk/" + sample + "_R2_cleaned.fastq.gz")
#    return ls


def _fastqc(config):
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/fastqc/" + sample + "_R1_fastqc.html")
        ls.append("analysis/" + config["token"] + "/fastqc/" + sample + "_R2_fastqc.html")
        #ls.append("analysis/" + config["token"] + "/fastqc/" + sample + "_R1_cleaned_fastqc.html")
        #ls.append("analysis/" + config["token"] + "/fastqc/" + sample + "_R2_cleaned_fastqc.html") 
    return ls


def _getSTARaligns(config):
    """ensure that bam and indexes are built"""
    #STAR alignment sorted.bam file, its index, and transcript count
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/STAR/" + sample + "/" + sample + ".sorted.bam")
        ls.append("analysis/" + config["token"] + "/STAR/" + sample + "/" + sample + ".sorted.bam.bai")
    ls.append("analysis/" + config["token"] + "/STAR/STAR_Align_Report.png")   
    return ls


def _rRNAmetrics(config):
    ls = []
    ls.append("analysis/" + config["token"] + "/STAR_rRNA/STAR_rRNA_Align_Report.csv")
    return ls


def _readQC(config):
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/RSeQC/read_distrib/" + sample + ".txt")
        ls.append("analysis/" + config["token"] + "/RSeQC/gene_body_cvg/" + sample + "/" + sample + ".geneBodyCoverage.curves.png")
        ls.append("analysis/" + config["token"] + "/RSeQC/junction_saturation/" + sample + "/" + sample + ".junctionSaturation_plot.pdf")
       # ls.append("analysis/RSeQC/insert_size/"+sample+"/"+sample+".histogram.pdf")
    ls.append("analysis/" + config["token"] + "/RSeQC/read_distrib/read_distrib.matrix.tab")
    ls.append("analysis/" + config["token"] + "/RSeQC/read_distrib/read_distrib.png")
    ls.append("analysis/" + config["token"] + "/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png") 
    return ls


def _get_cuff(config):
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/cufflinks/" + sample + "/" + sample + ".genes.fpkm_tracking")
        ls.append("analysis/" + config["token"] + "/cufflinks/" + sample + "/" + sample + ".isoforms.fpkm_tracking")
    ls.append("analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.csv")
    ls.append("analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.gct")
    ls.append("analysis/" + config["token"] + "/cufflinks/Cuff_Isoform_Counts.csv")
    return ls


def _get_cluster(config):
    ls = []
    #ls.append("analysis/" + config["token"] + "/plots/pca_plot.pdf")
    ls.append("analysis/" + config["token"] + "/plots/heatmapSS_plot.pdf")
    ls.append("analysis/" + config["token"] + "/plots/heatmapSF_plot.pdf")
    #for metacol in config["metacols"]:
    #    ls.append("analysis/" + config["token"] + "/plots/images/pca_plot_" + metacol + ".png")
    return ls


def _getDE(config):
    ls=[]
    for comparison in config['comps']:
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/"+comparison+".limma.csv")
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/"+comparison+".deseq.csv")
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/"+comparison+".deseq.sum.csv")
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/"+comparison+".limma.annot.csv")
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/"+comparison+".deseq.annot.csv")
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/"+comparison+".deseq.rlog.csv")
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/deseq_limma_fc_corr.csv")
        ls.append("analysis/" + config["token"] + "/diffexp/"+comparison+"/deseq_limma_fc_corr.png")
        ls.append("analysis/" + config["token"] + "/diffexp/" + comparison + "/" + comparison +"_volcano.pdf")
    ls.append("analysis/" + config["token"] + "/diffexp/de_summary.csv")
    ls.append("analysis/" + config["token"] + "/diffexp/de_summary.png")
    return ls


def _getpathway(config):
    ls=[]
    for comparison in config['comps']:
        ls.append("analysis/" + config["token"] + "/diffexp/" + comparison + "/" + comparison + ".goterm.done")
        #ls.append("analysis/" + config["token"] + "/diffexp/" + comparison + "/" + comparison + ".kegg.done")
    return ls


def _getvarscan2mpileup(config):
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/varscan2mpileup/"+sample+".snp-indel.vcf")
    return ls


def _preparation_input_files_for_rnaseqcnv(config):
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/rnaseqcnv/input/"+sample+".counts")
        ls.append("analysis/" + config["token"] + "/rnaseqcnv/input/"+sample+".vcf")
    ls.append("analysis/" + config["token"] + "/rnaseqcnv/input/config.txt")
    ls.append("analysis/" + config["token"] + "/rnaseqcnv/input/metadata.txt")
    return ls

def _run_rnaseqcnv(config):
    ls = []
    for sample in SAMPLES:
        ls.append("analysis/" + config["token"] + "/rnaseqcnv/{sample}/{sample}_CNV_main_fig.png")
    return ls


def _report_html(config):
    ls=[]
    for comparison in config['comps']:
        ls.append("analysis/" + config["token"] + "/report/" + comparison + ".report.done")
    return ls


