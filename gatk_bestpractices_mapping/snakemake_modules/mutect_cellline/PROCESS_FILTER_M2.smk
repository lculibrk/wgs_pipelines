rule PILEUP_SUMMARIES:
    input:
        BAM = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])[0],
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        SUMMARY = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/pileup.txt",
    singularity: "docker://lculibrk/gatk_alignment:latest"
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/pileup.log",
    benchmark:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/pileup.benchmark",
    resources:
        mem_mb = 3000,
        slurm_partition = config["partitions"]["short"],
        runtime = 60 * 4,
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.BAM} \
            -R {input.fa} \
            -V small_exac_common_3.vcf \
            -L small_exac_common_3.vcf \
            -O {output.SUMMARY} &> {log}
        """

rule COMBINE_MUTECT2_CELLLINE_VCFS:
    input:
        FILES = lambda wildcards: expand("data/studies/{{study}}/samples/{{sample}}/analyses/MUTECT_CELLLINE/{{analysis}}/mutect2/{{reference}}/{{type}}/raw_{chrom}.vcf", chrom = chromosomes[wildcards.reference])
    output:
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/merged.vcf"
    params:
        inputlist = lambda wildcards, input: " -I ".join([input]) if isinstance(input, str) else " -I ".join(input.FILES)
    singularity: "docker://lculibrk/gatk_alignment:latest"
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/merge_vcf.log"
    resources:
        slurm_partition = config["partitions"]["short"],
        runtime = 240,
        mem_mb = 3000
    shell:
        "gatk MergeVcfs -I {params.inputlist} -O {output.vcf} &> {log}"

rule MERGE_MUTECT2_CELLLINE_STATS:
    input:
        FILES = lambda wildcards: expand("data/studies/{{study}}/samples/{{sample}}/analyses/MUTECT_CELLLINE/{{analysis}}/mutect2/{{reference}}/{{type}}/raw_{chrom}.vcf.stats", chrom = chromosomes[wildcards.reference])
    output:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/merged.vcf.stats",
    params:
        inputlist = lambda wildcards, input: " --stats ".join([input]) if isinstance(input, str) else " --stats ".join(input.FILES)
    singularity: "docker://lculibrk/gatk_alignment:latest"
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/stats.log",
    resources:
        slurm_partition = config["partitions"]["short"],
        runtime = 240,
        mem_mb = 3000
    shell:
        """
        gatk MergeMutectStats \
            -stats {params.inputlist} -O {output}
        """


rule MUTECT2_CELLLINE_READ_ORIENTATION:
    input:
        files = lambda wildcards: expand("data/studies/{{study}}/samples/{{sample}}/analyses/MUTECT_CELLLINE/{{analysis}}/mutect2/{{reference}}/{{type}}/f1r2_{chrom}.tar.gz", chrom = chromosomes[wildcards.reference])
    output:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/read_orientation_model.tar.gz",
    singularity: "docker://lculibrk/gatk_alignment:latest"
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/read_orientation_model.log",
    resources:
        mem_mb = 3000,
        slurm_partition = config["partitions"]["short"],
        runtime = 60 * 4,
    params:
        inputlist = lambda wildcards, input: " -I ".join([input]) if isinstance(input, str) else " -I ".join(input.files)
    shell:
        "gatk LearnReadOrientationModel -I {params.inputlist} -O {output} &> {log}"


rule MUTECT2_CELLLINE_CALCULATE_CONTAMINATION:
    input:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/pileup.txt"
    output:
        contam = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/contamination.table",
        segments = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/segments.table",
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/contamination.log"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    resources:
        mem_mb = 3000,
        slurm_partition = config["partitions"]["short"],
        runtime = 60 * 4,
    shell:
        """
        gatk CalculateContamination -I {input} \
            -O {output.contam} --tumor-segmentation {output.segments} &> {log}
        """    


rule MUTECT2_CELLLINE_FILTERMUTECTCALLS:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        raw = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/merged.vcf",
        contam = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/contamination.table",
        segments = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/segments.table",
        mergedstats = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/merged.vcf.stats",
        orientation = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/read_orientation_model.tar.gz",
    output:
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered.vcf",
        stats = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filter.stats",
    singularity: "docker://lculibrk/gatk_alignment:latest"
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered.log"
    resources:
        mem_mb = 5000,
        slurm_partition = config["partitions"]["short"],
        runtime = 240,
    shell:
        """
        gatk FilterMutectCalls -V {input.raw} \
            -R {input.fa} \
            -O {output.vcf} \
            --contamination-table {input.contam} \
            --tumor-segmentation {input.segments} \
            --ob-priors {input.orientation} \
            --filtering-stats {output.stats} &> {log}
        """