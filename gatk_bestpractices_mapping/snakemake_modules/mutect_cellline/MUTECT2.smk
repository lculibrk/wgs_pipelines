##
#Standard Mutect2 run on a parent-daughter cell line
rule MUTECT2_DAUGHTER_STANDARD:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])[0],
        parent_merge = lambda wildcards: parent_cell(wildcards)
    output:
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/raw_{chrom}.vcf",
        stats= "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/raw_{chrom}.vcf.stats",
        tgz = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/f1r2_{chrom}.tar.gz"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    resources:
        threads = 1,
        mem_mb = lambda wildcards, attempt: 4000 * (1 + ((attempt-1)/2)),
        iotasks = 2,
        runtime = 24*60,
        slurm_partition = config["partitions"]["medium"],
        att = lambda wildcards, attempt: attempt
    benchmark:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/variants_{chrom}.resources",
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/variants_{chrom}.log"
    shell:
        """
            NORMAL_NAME=$(gatk GetSampleName -I {input.parent_merge} -R {input.fa} -O /dev/stdout 2> {log}); \
            gatk --java-options '-Xmx{resources.mem_mb}M' \
                Mutect2 -R {input.fa} \
                -L {wildcards.chrom} \
                -I {input.cell_merge} \
                -I {input.parent_merge} \
                --normal $NORMAL_NAME \
                --panel-of-normals pon.vcf \
                --germline-resource gnomad.vcf \
                --native-pair-hmm-threads 1 \
                --f1r2-tar-gz {output.tgz} \
                --output {output.vcf} &>> {log}.{resources.att}
        """

##
#Comprehensive mutect2 run incl. germline sites, for assessing sharedness etc between daughters and for genotyping parents
rule MUTECT2_DAUGHTER_GERMLINE_SITES:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])[0],
        parent_merge = lambda wildcards: parent_cell(wildcards)
    output:
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/germ/raw_{chrom}.vcf",
        stats= "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/germ/raw_{chrom}.vcf.stats",
        tgz = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/germ/f1r2_{chrom}.tar.gz"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    resources:
        threads = 1,
        mem_mb = lambda wildcards, attempt: 4000 * (1 + ((attempt-1)/2)),
        iotasks = 2,
        runtime = 24*60,
        slurm_partition = config["partitions"]["medium"],
        att = lambda wildcards, attempt: attempt
    benchmark:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/germ/variants_{chrom}.resources",
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/germ/variants_{chrom}.log"
    shell:
        """
            NORMAL_NAME=$(gatk GetSampleName -I {input.parent_merge} -R {input.fa} -O /dev/stdout 2> {log}) ; \
            gatk --java-options '-Xmx{resources.mem_mb}M' \
                Mutect2  -R {input.fa} \
                -L {wildcards.chrom} \
                -I {input.cell_merge} \
                -I {input.parent_merge} \
                --normal $NORMAL_NAME \
                --genotype-germline-sites \
                --panel-of-normals pon.vcf  \
                --germline-resource gnomad.vcf \
                --native-pair-hmm-threads 1 \
                --f1r2-tar-gz {output.tgz} \
                --output {output.vcf} &>> {log}.{resources.att}
        """

##
# Super comprehensive mutect2 run to genotype all parental sites where a mutation was found in any daughter, for assessing sharedness between parental lineages and daughters
rule MUTECT2_PARENTAL_FORCE_SITES:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])[0],
        daughter_vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/daughters_merged_vcf/{chrom}.vcf",
    output:
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parental/raw_{chrom}.vcf",
        stats= "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parental/raw_{chrom}.vcf.stats",
        tgz = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parental/f1r2_{chrom}.tar.gz"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    resources:
        threads = 1,
        mem_mb = lambda wildcards, attempt: 4000 * (1 + ((attempt-1)/2)),
        iotasks = 2,
        runtime = 24*60,
        slurm_partition = config["partitions"]["medium"],
        att = lambda wildcards, attempt: attempt
    benchmark:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parental/variants_{chrom}.resources",
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parental/variants_{chrom}.log"
    shell:
        """
            gatk --java-options '-Xmx{resources.mem_mb}M' \
                Mutect2 -R {input.fa} \
                -L {wildcards.chrom} \
                -I {input.merge} \
                --alleles {input.daughter_vcf} \
                --force-call-filtered-alleles \
                --genotype-germline-sites \
                --panel-of-normals pon.vcf  \
                --germline-resource gnomad.vcf \
                --native-pair-hmm-threads 1 \
                --f1r2-tar-gz {output.tgz} \
                --output {output.vcf} &> {log}.{resources.att}
        """

