import os
import numpy as np
import glob

MUTECT_CELLLINE_PIPELINE_VERSION = "1.0.0"
MUTECT_CELLLINE_PATH = os.path.join(basedir, "snakemake_modules", "mutect_cellline")


include: os.path.join(MUTECT_CELLLINE_PATH, "MUTECT2.smk")
include: os.path.join(MUTECT_CELLLINE_PATH, "PROCESS_FILTER_M2.smk")
include: os.path.join(MUTECT_CELLLINE_PATH, "MUTECT2_PARENTAL.smk")


def parent_cell(wildcards):
    sample_id = genomicsapi.translate.stringtoid(wildcards.sample)
    db_line = genomicsapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample_id, bench = config["bench"])
    #print(db_line)
    parent_id = genomicsapi.translate.idtostring(db_line[0][5], "MPS")
    parent_merge = gateway("WGS_MERGE_BAM", parent_id, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])[0]
    return(parent_merge)

rule MUTECT2_CELLLINE_VCFTOTABLE:
    input:
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered_renamed.vcf",
    output:
        tab = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/table_raw.txt"
    resources:
        mem_mb = 3000,
        slurm_partition = config["partitions"]["short"],
        runtime = 60 * 4,
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/table.log"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    shell:
        """
        gatk VariantsToTable -V {input.vcf} \
            -O {output.tab} -raw --split-multi-allelic &> {log}
        """

def get_parents_for_sample(wildcards):
    ## Get parent "comprehensive" VCF paths
    ## 1. API call to get the daughter's parent ID
    ## 2. For that parent ID, get all samples from that study with the same cell line ID and parent ID
    #
    ## Get daughter sample table line
    study_id = genomicsapi.translate.stringtoid(wildcards.study)
    daughter_id = genomicsapi.translate.stringtoid(wildcards.sample)
    daughter_line = genomicsapi.select.multi_select(db = db, table = "samples", filters = {"id":daughter_id}, headers = True, bench = config["bench"])
    ## Get the line for that daughter's annotated parent and get the cell ID from that
    parent_id = np.array(daughter_line[1])[np.where(np.array(daughter_line[0]) == "sample_parent_id")][0]
    cell_id = np.array(daughter_line[1])[np.where(np.array(daughter_line[0]) == "cell_id")][0]
    ## Now get the parent of the parent, so we can get all the possible parents
    parental_sample = genomicsapi.select.multi_select(db = db, table = "samples", filters = {"id":parent_id}, headers = True, bench = config["bench"])
    parental_parent = np.array(parental_sample[1])[np.where(np.array(daughter_line[0]) == "sample_parent_id")][0]
    ## Now get all the appropriate parental lines
    all_parental_lines = genomicsapi.select.multi_select(db = db, table = "samples", filters = {"study_id":study_id, "sample_parent_id":parental_parent, "cell_id":cell_id}, bench = config["bench"])
    ## Loop through to get the analysis paths and thereby the path to the parents' variant table
    parental_paths = []
    for parent in all_parental_lines:
        ## First gateway to ensure the parent is made
        gateway("MUTECT", genomicsapi.translate.idtostring(parent[0], "MPS"), SCRATCH_DIR, config["PROD_DIR"], db = config["db"])
        ## Get the analysis dir
        analysis = genomicsapi.select.multi_select(db = db, table = "analyses", filters = {"samples_id":parent[0], "pipeline_name":"MUTECT_CELLLINE"}, bench = config["bench"])
        analysis = analysis[0]
        vcf_path = os.path.join(analysis[10], genomicsapi.translate.idtostring(analysis[0], "MPA"), "mutect2", wildcards.reference, "parental/table_raw.txt")
        parental_paths.append(vcf_path)
    return(parental_paths)



rule MUTECT2_CELLLINE_PARENT_TABLE:
    input:
        get_parents_for_sample
    output:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parents.table"
    resources:
        mem_mb = 2000,
        slurm_partition = config["partitions"]["short"],
        runtime = 5
    shell:
        """
        for i in {input}; do echo $i >> {output}; done;
        """

def parent_sample_id(wildcards):
    sample_id = genomicsapi.translate.stringtoid(wildcards.sample)
    db_line = genomicsapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample_id, bench = config["bench"])
    #print(db_line)
    parent_id = genomicsapi.translate.idtostring(db_line[0][5], "MPS")
    return(parent_id)

def full_cell_line_cohort(wildcards):
    study_id = genomicsapi.translate.stringtoid(wildcards.study)
    samples_ids = genomicsapi.select.multi_select(db, "samples", {"study_id":study_id}, bench = config["bench"])
    o_list = []
    for sample in samples_ids:
        #o_list.extend(gateway("MUTECT", genomicsapi.translate.idtostring(sample[0], "MPS"), scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = db))
        ## get the analysis ID for this sample
        analysis = genomicsapi.select.multi_select(db, "analyses", {"samples_id":sample[0], "pipeline_name":"MUTECT_CELLLINE"})
        if len(analysis) < 1:
            print(sample)
            continue
        ## TODO: make this hg-agnostic
        analysis = analysis[0]
        p = analysis[10]
        aid = genomicsapi.translate.idtostring(analysis[0], "MPA")
        o_list.extend(glob.glob(os.path.join(p, aid, "mutect2/hg19/*/table_raw.txt")))
    return(o_list)
    
rule M2_SBS_TABLES:
    output:
        daughter = "data/studies/{study}/analyses/MUTECT_CELLLINE/daughters_table.txt",
        parental = "data/studies/{study}/analyses/MUTECT_CELLLINE/parents_table.txt",
    threads: 1
    resources:
        slurm_partition = config["partitions"]["short"],
        cpus = 1,
        threads = 1,
        mem_mb = 2000,
        runtime = 15,
    run:
        if not os.path.exists(f"data/studies/{wildcards.study}/analyses/MUTECT_CELLLINE/"):
            os.makedirs(f"data/studies/{wildcards.study}/analyses/MUTECT_CELLLINE/")
        study_id = genomicsapi.translate.stringtoid(wildcards.study)
        samples_ids = genomicsapi.select.multi_select(db, "samples", {"study_id":study_id})
        for sample in samples_ids:
            cell_name = genomicsapi.select.multi_select(db, "cells", {"id":sample[6]})[0][1]
            ## Get analysis row
            analysis = genomicsapi.select.multi_select(db = db, table = "analyses", filters = {"samples_id":sample[0], "pipeline_name":"MUTECT_CELLLINE"}, bench = config["bench"])
            if len(analysis) == 0:
                continue
            analysis = analysis[0]
            aid = genomicsapi.translate.idtostring(analysis[0], "MPA")
            ap = analysis[10]
            ## Test if this sample is used as a reference for any other samples
            daughters = genomicsapi.select.simple_select(db = db, table = "samples", filter_column = "sample_parent_id", filter_value = sample[0], bench = config["bench"])
            parent_id = sample[5]
            if parent_id is not None and parent_id != "NULL": ## if this sample has a parent ie. it's a daughter
                files = ["mutect2/hg19/std/table_raw.txt", "mutect2/hg19/germ/table_raw.txt"]
                files = [os.path.join(ap, aid, f) for f in files]
                parent_name = genomicsapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = parent_id)[0][1]
                with open(output.daughter, "a") as f:
                    f.write("%s\n" % '\t'.join([files[0], files[1], sample[1], parent_name, cell_name]))
            if len(daughters) > 0: ## ie. if this is a parent
                files = ["mutect2/hg19/parental/table_raw.txt"]
                with open(output.parental, "a") as f:
                    f.write("%s\n" % '\t'.join([files[0], sample[1], cell_name]))



rule MUTECT2_CELLLINE_SBS_GROUP_FILTER:
    input:
        daughter = "data/studies/{study}/analyses/MUTECT_CELLLINE/daughters_table.txt",
        parental = "data/studies/{study}/analyses/MUTECT_CELLLINE/parents_table.txt",
        data = full_cell_line_cohort
    output:
        "data/studies/{study}/analyses/MUTECT_CELLLINE/filtering_done.txt"
    threads: 4
    resources:
        mem_mb = 150000,
        slurm_partition = config["partitions"]["highmem"],
        threads = 4,
        cpus = 4,
        runtime = 720,
    params:
        script_path = MUTECT_CELLLINE_PATH + "/scripts/filter_mutations.R"
    log:
        "data/studies/{study}/analyses/MUTECT_CELLLINE/filtering.log"
    shell:
        """
        module load r/4.1.2;
        Rscript {params.script_path} -d {input.daughter} -p {input.parental} -t {resources.threads} 2> {log};
        echo 'Rscript {params.script_path} -d {input.daughter} -p {input.parental} -t {resources.threads}' > {output}
        """

rule PROC_FILE:
    input:
        filt = "data/studies/{study}/analyses/MUTECT_CELLLINE/filtering_done.txt",
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/filtered_renamed.vcf"
    output:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/variants_final.vcf"
    threads: 1
    resources:
        mem_mb = 3000,
        slurm_partition = config["partitions"]["short"],
        threads = 1,
        cpus = 1,
        runtime = 10,
    params:
        script_path = os.path.join(MUTECT_CELLLINE_PATH, "scripts/cellline_vcf.py"),
        proc_vcf = lambda w: f"data/studies/{w.study}/samples/{w.sample}/analyses/MUTECT_CELLLINE/{w.analysis}/mutect2/{w.reference}/proc/variants.vcf",
        cont = "docker://lculibrk/gatk_alignment:latest"
    log:
         "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/variants_final.log"
    shell:
        """
        cat <(python {params.script_path} -v {input.vcf} -c {input.filt}) {params.proc_vcf} > {output} &&
        singularity run -B / {params.cont} gatk ValidateVariants -V {output} &> {log}
        """
        

rule CELLLINE_OF_ORIGIN:
    input:
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])[0],
        FA = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        vcf = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/qc_genotyping.vcf",
        tbl1 ="data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/qc_genotyping.txt",
        tbl2 ="data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/line_of_origin.txt",
    singularity:
        "docker://lculibrk/bcftools/bcftools_latest.sif"
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/qc_genotyping.log",
    resources:
        slurm_partition = config["partitions"]["short"],
        runtime = 240,
        mem_mb = 3000
    shell:
        """
        bcftools mpileup \
         -R resources/genotype_positions.txt \
         --fasta-ref {input.FA} \
         -A {input.cell_merge} 2> {log} | bcftools call -c 2>> {log} > {output.vcf};
        gatk VariantsToTable -V {output.vcf} -O {output.tbl1} 2>> {log};
        module load r/4.1.2;
        cat {output.tbl1} | Rscript scripts/qc_identity.R > {output.tbl2} 2>> {log}
        """

rule MUTECT2_CELLLINE_PAIRED_FILTER:
    input:
        germ = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/germ/table_raw.txt",
        std = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/table_raw.txt",
        tbl = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parents.table"
    output:
        tab = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/variants.txt"
    resources:
        mem_mb = 100000,
        slurm_partition = config["partitions"]["highmem"],
        runtime = 30
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/table.log"
    params:
        parent = lambda w: parent_sample_id(w)
    shell:
        "module load r/4.1.2; Rscript scripts/filter_paired.R --germline {input.germ} --input {input.std} --table {input.tbl} --parent {params.parent} --outtable {output.tab} 2> {log}"

## At the end of all this, the invocation of this module should be able to figure out if a sample is purely parental or not. If it's purely parental, then it shouldn't trigger the pipeline execution. 