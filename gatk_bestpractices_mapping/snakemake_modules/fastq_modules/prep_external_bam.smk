import logging
import re
import shutil
import os

localrules: COUNT_RG_ASSIGN_RUNS

rule COUNT_RG_ASSIGN_RUNS:
    input:
        get_external_bam_path ## Need to write a custom function to return all the bams to be loaded into the runs table
    output:
        "data/studies/{study}/samples/{sample}/analyses/LOAD_EXT_BAM/{analysis}/db_load_runs/runs.loaded"
    params:
        workdir = "data/studies/{study}/samples/{sample}/analyses/LOAD_EXT_BAM/{analysis}/db_load_runs",
        db = config["db"]
    shell:
        """
            samtools view -H {input} > {params.workdir}/header.txt
            NRUNS=$(grep '^@RG' {params.workdir}/header.txt | wc -l)
            for ((i = 0; i < NRUNS; i++)); do
                python scripts/add_element.py -t run -r {wildcards.sample}-$i --study {wildcards.study} --sample {wildcards.sample} --db {params.db} --source external_bam
                echo "Iteration $((i + 1))"
            done
            touch {output}
        """