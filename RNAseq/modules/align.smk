strand_command=""
rRNA_strand_command=""

strand_command="--outFilterIntronMotifs RemoveNoncanonical"
rRNA_strand_command="--outFilterIntronMotifs RemoveNoncanonical"
gz_command="--readFilesCommand zcat"
_mates = ['mate1', 'mate2'] 
_keepPairs = "KeepPairs" if len(_mates) == 2 else ""


rule run_STAR:
    input:
        unpack(get_fastq)
    output:
        bam = protected("analysis/{token}/STAR/{sample}/{sample}.sorted.bam"),
        counts = "analysis/{token}/STAR/{sample}/{sample}.counts.tab",
        log_file = "analysis/{token}/STAR/{sample}/{sample}.Log.final.out",
        sjtab = "analysis/{token}/STAR/{sample}/{sample}.SJ.out.tab"
    threads: 10
    resources:
        mem_mb=8000
    params:
        stranded=strand_command,
        gz_support=gz_command,
        prefix=lambda wildcards: "analysis/{token}/STAR/{sample}/{sample}".format(sample=wildcards.sample,token=config["token"]),
        readgroup=lambda wildcards: "ID:{sample} PL:ILLUMINA LB:{sample} SM:{sample}".format(sample=wildcards.sample),
        keepPairs=_keepPairs,
        star_index=config['star_index'],
        gtf=config['gtf_file'],
        threads = 10
    benchmark:
        "analysis/{token}/benchmarks/{sample}/{sample}.run_STAR.txt"
    log: 
        stderr = "analysis/{token}/logs/{sample}/star/{sample}_STAR.stderr.log",
        stdout = "analysis/{token}/logs/{sample}/star/{sample}_STAR.stdout.log"
    shell:
        "STAR" 
        " --runMode alignReads"
        " --runThreadN {params.threads}"
        " --genomeDir {params.star_index}"
        " --readFilesIn {input} {params.gz_support}" 
        " --outFileNamePrefix {params.prefix}."
        " --outSAMstrandField intronMotif"
        " --outSAMmode Full --outSAMattributes All {params.stranded}"
        " --outSAMattrRGline {params.readgroup}"
        " --outSAMtype BAM SortedByCoordinate"
        " --limitBAMsortRAM 60000000000" 
        " --quantMode GeneCounts"
        " --twopassMode Basic"
        " --sjdbGTFfile {params.gtf}"
        " 1> {log.stdout} 2> {log.stderr}"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"


rule index_bam:
    """INDEX the {sample}.sorted.bam file"""
    input:
        "analysis/{token}/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/{token}/STAR/{sample}/{sample}.sorted.bam.bai"
    message: "Indexing {wildcards.sample}.sorted.bam"
    benchmark:
        "analysis/{token}/benchmarks/{sample}/{sample}.index_bam.txt"
    shell:
        "samtools index {input}"


rule generate_STAR_report:
    input:
        star_log_files=expand( "analysis/{token}/STAR/{sample}/{sample}.Log.final.out", sample=SAMPLES, token=config["token"]),
        star_gene_count_files=expand( "analysis/{token}/STAR/{sample}/{sample}.counts.tab", sample=SAMPLES, token=config["token"])
    output:
        csv="analysis/{token}/STAR/STAR_Align_Report.csv",
        gene_counts="analysis/{token}/STAR/STAR_Gene_Counts.csv"
    message: "Generating STAR report"
    benchmark:
        "analysis/{token}/benchmarks/generate_STAR_report.txt"
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell("perl modules/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}")
        shell( "perl modules/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )

