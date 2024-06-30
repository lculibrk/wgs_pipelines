import os
import numpy as np


rule UBAM:
    ## Input: pair of R1/R2 fastq files, output: unaligned BAM
    input:
        lambda wildcards: gateway("FASTQ", wildcards.run, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])
    output:
        ubam = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    log:
        "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam.log"
    resources:
        runtime = 60,
        mem_mb = 3000,
        slurm_partition = config["partitions"]["short"],
    params:
        TEMP_DIR = config["temp_dir"]
        TEMP_BAM_PATH = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis +  "/markdups/",
    priority: 9999
    shell:
        """
        gatk FastqToSam \
            --FASTQ {input[0]} \
            --FASTQ2 {input[1]} \
            --OUTPUT {output.ubam} \
            --READ_GROUP_NAME {wildcards.run} \
            --SAMPLE_NAME {wildcards.sample} \
            --LIBRARY_NAME {wildcards.run} \
            --PLATFORM ILLUMINA > {log} 2>&1
        """

rule MARK_ADAPTERS:
    ## Marks reads with adapter sequences
    input:
        ubam = rules.UBAM.output
    output:
        obam = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.bam"),
        metrics = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.metrics"
    log:
        "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.bam.log"
    params:
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/markadapts/",
    benchmark: 
        "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.benchmark"
    resources:
        runtime = 240,
        mem_mb = 2000,
        iotasks = 2,
        slurm_partition = config["partitions"]["short"],
    singularity: "docker://lculibrk/gatk_alignment:latest"
    priority: 1000
    shell:
        """
        mkdir -p {params.tmp}
        gatk MarkIlluminaAdapters \
            -I {input.ubam} \
            -O {output.obam} \
            -M {output.metrics} \
            -TMP_DIR {params.tmp} 2> {log}
        rm -rf {params.tmp}
        """

rule MARKED_SAM_TO_FASTQ:
    input:
        bam = rules.MARK_ADAPTERS.output.obam,
    output:
        ifastq = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/interleaved.fastq")
    log:
        samtofastq = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.samtofastq.log",
    params:
        align_threads = ALIGN_THREADS - 2,
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/marktofq/",
    singularity: "docker://lculibrk/gatk_alignment:latest"
    benchmark: "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.samtofastq.benchmark"
    resources:
        runtime = 480,
        iotasks = 2,
        slurm_partition = config["partitions"]["short"],
    priority: 1001
    shell:
        """
        mkdir -p {params.tmp}
        gatk SamToFastq \
            -I {input.bam} \
            -FASTQ {output.ifastq} \
            -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
            -TMP_DIR {params.tmp} 2> {log.samtofastq}
        rm -rf {params.tmp}
        """

rule BWA:
    input:
        ifastq = rules.MARKED_SAM_TO_FASTQ.output.ifastq,
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
    output:
        rawbam = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.raw.sam")
    threads: ALIGN_THREADS
    resources:
        threads = ALIGN_THREADS,
        cpus = ALIGN_THREADS,
        mem_mb = 30000,
        iotasks = 4,
        slurm_partition = config["partitions"]["medium"],
        tmpdisk = lambda wc, input: int(np.round(input.size_mb))
    log:
        bwa = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.log",
    singularity: "docker://lculibrk/gatk_alignment:latest"
    benchmark: "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.benchmark"
    params:
        tmpfastqpath = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/",
        tmpfastq = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/input.fastq"
    priority: 1002
    shell:
        """
            set +eo pipefail;
            mkdir -p {params.tmpfastqpath};
            cp {input.ifastq} {params.tmpfastq};
            bwa mem -M -t {threads} -p {input.bwa_idxbase} {params.tmpfastq} > {output} 2> {log.bwa};
            if [ $? -ne 0 ]; then
                rm {params.tmpfastq};
                set -eo pipefail;
                exit 69;
            fi;
            rm {params.tmpfastq}
        """

rule MERGEBAMALIGNMENT:
    input:
        bwa = rules.BWA.output.rawbam,
        fa_base = lambda wildcards: FA_PATHS[wildcards.reference],
        ubam = rules.UBAM.output,
    output:
        bam = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bam")
    log:
        mergebamalignment = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.mergebamalignment.log"
    params:
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/mergebamalignment/",
    resources:
        threads = 1,
        cpus = 1,
        mem_mb = 3000,
        iotasks = 4,
        slurm_partition = config["partitions"]["short"],
    singularity: "docker://lculibrk/gatk_alignment:latest"
    benchmark: "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.mergebamalignment.benchmark"
    priority: 9999
    shell:
        """
        mkdir -p {params.tmp}
        gatk MergeBamAlignment \
            -ALIGNED_BAM {input.bwa} \
            -UNMAPPED_BAM {input.ubam} \
            -OUTPUT {output.bam} \
            -R {input.fa_base} -CREATE_INDEX true -ADD_MATE_CIGAR true \
            -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true \
            -INCLUDE_SECONDARY_ALIGNMENTS true -MAX_INSERTIONS_OR_DELETIONS -1 \
            -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS \
            -TMP_DIR {params.tmp} 2> {log.mergebamalignment}
        rm -rf {params.tmp}
        """

rule GATK_MARKDUPS:
    ## Mark duplicate reads/alignments
    input:
        rules.MERGEBAMALIGNMENT.output,
    output:
        bam = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam"),
#        metrics = "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.metrics"
    threads: 1    
    params:
        TEMP_DIR = TEMP_DIR,
        TEMP_BAM_PATH = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis +  "/markdups/",
        TEMP_BAM_FILE = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/markdups/input.bam"
    resources:
        threads = 8,
        cpus = ALIGN_THREADS,
        mem_mb = 35000,
        iotasks = 4,
        slurm_partition = config["partitions"]["medium"],
        runtime = 60*24*3,
        tmpdisk = lambda wc, input: int(np.round(input.size_mb))
    log:
        "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.log"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    benchmark: "data/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.benchmark"
    priority: 1003
    shell:
        """
        set +eo pipefail;
        rm -rf {output.bam}.parts/;
        mkdir -p {params.TEMP_BAM_PATH};
        cp {input} {params.TEMP_BAM_FILE};
        gatk MarkDuplicatesSpark \
            -I {params.TEMP_BAM_FILE} \
            -O {output.bam} \
            --conf 'spark.executor.cores={resources.threads}' \
            --conf 'spark.local.dir={params.TEMP_DIR}' &> {log};
        if [ $? -ne 0 ]; then
            rm {params.TEMP_BAM_FILE};
            set -eo pipefail;
            exit 69;
        fi
        rm -rf {params.TEMP_BAM_FILE}
        """

rule GATK_MERGE:
    ## Merge multi-run samples into single bams
    ## Currently merges ALL bam runs from a sample into a single bam. If we want subsets, will need to add that functionality later
    input:
        aggregate_runs
    output:
        merge = temp("data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam")
    log:
        "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.log"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    threads: 2
    resources:
        threads = 2,
        cpus = 2,
        runtime = 480,
        iotasks = 4,
        slurm_partition = config["partitions"]["short"],
    params:
        inputlist = lambda wildcards, input: f"-I {input}" if isinstance(input, str) else "-I " + " -I ".join(input)
    benchmark: "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.benchmark"
    priority: 1004
    shell:
        """
        gatk MergeSamFiles \
            {params.inputlist} -O {output.merge} \
            --USE_THREADING true \
             &>> {log}
        """

rule BAM2CRAM:
    input:
        merge = rules.GATK_MERGE.output,
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
    output:
        cram = "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram",
        crai = "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.crai",
        ref_md5 = "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.ref.md5",
        readme = "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/CRAM_README.txt",
    log:
        "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.log"
    singularity: "docker://lculibrk/bwa-mem2:latest"
    threads: 1
    resources:
        runtime = 480,
        slurm_partition = config["partitions"]["short"],
        mem_mb = 4000,
    benchmark: "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.bench"
    priority: 1005
    shell:
        """
        md5sum {input.bwa_idxbase} > {output.ref_md5};
        echo    'This CRAM was generated using the fasta file located at {input.bwa_idxbase}. The md5 hash of the reference is at {output.ref_md5}. \
                The matching reference fasta at the specified directory is REQUIRED for proper decompression of the file' > {output.readme};
        samtools view -@ {threads} -C -T {input.bwa_idxbase} {input.merge} > {output.cram} 2> {log}
        chmod 750 {output};
        samtools index {output.cram} &>> {log}
        """

rule GATK_BAM_DONE:
    input:
        rules.BAM2CRAM.output,
    output:
        "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.done"
    log:
        "data/studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.log"
    params:
        db = config["db"]
    resources:
        runtime = 10,
        slurm_partition = config["partitions"]["short"],
    priority: 1007
    shell:
        "python scripts/mark_complete.py --id {wildcards.analysis} --db {params.db} {input} >> {log}; touch {output}"