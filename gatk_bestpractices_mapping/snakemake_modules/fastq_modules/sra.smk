SRA_PIPE_VERSION = "1.0.0"


rule PREFETCH:
    output:
        sra = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/SRA/{analysis}/sra/{sra}/{sra}.sra"),
    retries: 3
    params:
        SCRATCH_DIR = SCRATCH_DIR,
        SRA_PIPE_VERSION = SRA_PIPE_VERSION
    resources:
        iotasks = 2,
        runtime = 700,
        slurm_partition = config["partitions"]["short"],
    log:
        "data/studies/{study}/samples/{sample}/runs/{run}/analyses/SRA/{analysis}/sra/{sra}/{sra}.log"
    shell:
        "prefetch {wildcards.sra} -O data/studies/{wildcards.study}/samples/{wildcards.sample}/runs/{wildcards.run}/analyses/SRA/{wildcards.analysis}/sra/ -f ALL -r no &> {log};"

rule DUMP:
    input:
        sra = rules.PREFETCH.output
    output:
        temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/SRA/{analysis}/sra/{sra}_1.fastq"),
        temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/SRA/{analysis}/sra/{sra}_2.fastq")
    params:
        SCRATCH_DIR = SCRATCH_DIR,
        SRA_PIPE_VERSION = SRA_PIPE_VERSION,
        TEMP_PATH = lambda wildcards: f"{TEMP_DIR}{wildcards.sra}/"
    resources:
        iotasks = 2,
        slurm_partition = config["partitions"]["short"],
        runtime = 700,
    log:
        "data/studies/{study}/samples/{sample}/runs/{run}/analyses/SRA/{analysis}/sra/{sra}.log"
    shell:
        """
        fasterq-dump {input.sra} --split-files -t {params.TEMP_PATH} -O data/studies/{wildcards.study}/samples/{wildcards.sample}/runs/{wildcards.run}/analyses/SRA/{wildcards.analysis}/sra/ &>> {log}; rm -rf ./{wildcards.sra}/
        """

def sra_function(wildcards):
    result = genomicsapi.select.simple_select(db = db, table = "runs", filter_column = "id", filter_value = genomicsapi.translate.stringtoid(wildcards.run), bench = config["bench"])
    sra_id = result[0][6]
    if not sra_id or sra_id == "NULL":
        raise ValueError(f"SRA pipeline is being used, but run {wildcards.run} has no associated SRA ID in the database")
    fq_path = f"data/studies/{wildcards.study}/samples/{wildcards.sample}/runs/{wildcards.run}/analyses/SRA/{wildcards.analysis}/sra/{sra_id}_"
    fq_paths = [fq_path + num + ".fastq" for num in ["1", "2"]]
    return(fq_paths)

rule COPY_SRA:
    input:
        sra_function
    output:
        reads1 = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/SRA/{analysis}/fq/reads1.fq"),
        reads2 = temp("data/studies/{study}/samples/{sample}/runs/{run}/analyses/SRA/{analysis}/fq/reads2.fq")
    params:
        DB = db,
    resources:
        iotasks = 2,
        slurm_partition = config["partitions"]["short"],
        runtime = 700,
    shell:
        """
        mv {input[0]} {output.reads1}; mv {input[1]} {output.reads2}; \
        python scripts/mark_complete.py -i {wildcards.analysis} -d {params.DB} {output.reads1}
        """