import logging
import re
import shutil
import os
EXTERNAL_BAM_VERSION = "1.0.0"
GATK_ALIGNER_VER = "1.0.1"

def collate_output_fastqs(wildcards):
    out_path = checkpoints.SPLIT_BAM.get(**wildcards).output[0]
    newpath = re.sub("split_bams", "splitfq", out_path)
    return expand(newpath + "/{i}_1.fq",
                  i=glob_wildcards(os.path.join(out_path, "{i}.bam")).i)

def run_rname_to_rid(wildcards):
    sample_id = genomicsapi.translate.stringtoid(wildcards.sample)
    selection = genomicsapi.select.multi_select(db = db, table = "runs", filters = {"sample_id":sample_id, "rname":wildcards.run}, bench = config["bench"])
    return(selection[0][0])


rule SPLIT_BAM:
    input:
        get_external_bam_path ## Need to write a custom function to return the path to a bam based on wildcards
    output:
        temp(directory("data/studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/split_bams/"))
    singularity: "docker://lculibrk/bwa-mem2:latest"
    log:
        "data/studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/split_bams/split.log"
    threads: 8
    resources:
        threads = 8,
        cpus = 8,
        mem_mb = 8000,
        runtime = 240,
        slurm_partition = config["partitions"]["short"],
        iotasks = 5,
    shell:
        "samtools split {input} -f {output}/{wildcards.sample}-%#.%. -@ {resources.threads}"

def get_split_bam(wildcards):
    outs = checkpoints.SPLIT_BAM.get(**wildcards).output[0]
    ## Get sample ID
    rid = genomicsapi.translate.stringtoid(wildcards.run)
    rname = genomicsapi.select.multi_select(db = db, table = "runs", filters = {"id":rid}, bench = config["bench"])[0][1]
    ## Use sample ID to get run ID & rname
    bamfile = os.path.join(outs, f"{rname}.bam")
    #print(bamfile)
    return(bamfile)

rule EXTERNAL_BAM_TO_FASTQ:
    input:
        "data/studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/split_bams/"
        #get_split_bam
        #"data/studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/{analysis}/split_bams/{rname}.bam"
    output:
        reads1 = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/reads1.fq"),
        reads2 = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/reads2.fq"),
        readsU = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/readsU.fq",
        readsS = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/readsS.fq",
    singularity: "docker://lculibrk/bwa-mem2:latest"
    log:
        "data/studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/split.log",
    resources:
        iotasks = 2,
        mem_mb = 20000,
        runtime = 600,
        slurm_partition = config["partitions"]["short"],
    params:
        rname = lambda w: genomicsapi.select.multi_select(db = db, table = "runs", filters = {"id":genomicsapi.translate.stringtoid(w.run)}, bench = config["bench"])[0][1]
    priority: 1
    shell:
        "samtools collate -u -O {input}/{params.rname}.bam | samtools fastq -1 {output.reads1} -2 {output.reads2} -0 {output.readsU} -s {output.readsS} - &> {log};"

