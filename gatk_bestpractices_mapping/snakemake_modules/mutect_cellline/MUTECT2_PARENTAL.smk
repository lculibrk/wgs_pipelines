rule COPY_RENAME_DAUGHTER_CALLS:
    input:
        "data/studies/{study}/samples/{daughter_sample}/analyses/MUTECT_CELLLINE/{daughter_analysis}/mutect2/{reference}/germ/raw_{chrom}.vcf"
    output:
        "data/studies/{study}/samples/{daughter_sample}/analyses/MUTECT_CELLLINE/{daughter_analysis}/mutect2/{reference}/germ/renamed_{chrom}.vcf"
    log:
        "data/studies/{study}/samples/{daughter_sample}/analyses/MUTECT_CELLLINE/{daughter_analysis}/mutect2/{reference}/germ/renamed_{chrom}.log"
    resources:
        threads = 1,
        mem_mb = 5000,
        runtime = 15,
        slurm_partition = config["partitions"]["short"],
    shell:
    ## First copy over the header
    ## Next, grep out the header, cut out the parental column, then use awk to rename the daughter column to {parental-sample-name}_d, then remove variants with ref/alt length longer than the reads
    ## Otherwise, it does not play nice with mutect2
        """
        grep '##' {input} > {output} 2> {log};
        TUMOR_NAME=$(grep '##tumor_sample' {input} | sed 's/##tumor_sample=//g');
        column_number=$(awk -v search="$TUMOR_NAME" -F'\\t' '{{ for (i=1; i<=NF; i++) {{ if ($i == search) {{ print i; exit }} }} }}' {input})
        grep -v '##' {input} | cut -f 1,2,3,4,5,6,7,8,9,$column_number | awk -v OFS="\\t" 'NR==1{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, "daughter"}}NR>1 && length($4) == 1 && length($5) == 1{{print $0}}' >> {output} 2>> {log}
        """    

rule RENAME_SAMPLES_VCF:
    input:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered.vcf",
    output:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered_renamed.vcf",
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered_renamed.log",
    resources:
        slurm_partition = config["partitions"]["short"],
        runtime = 10,
        mem_mb = 1000
    shell:
        """
        set +o pipefail;
        TUMOR_NAME=$( grep '##tumor_sample' {input} | sed 's/##tumor_sample=//g');
        NORMAL_NAME=$(grep '##normal_sample' {input} | sed 's/##normal_sample=//g');
        grep '##' {input} > {output}
        if [ ! -z $NORMAL_NAME ]; then
            grep '#CHROM' {input} | sed "s/$TUMOR_NAME/tumor/g" | sed "s/$NORMAL_NAME/normal/g" >> {output};
        elif [ -z $NORMAL_NAME ]; then
            grep '#CHROM' {input} | sed "s/$TUMOR_NAME/tumor/g" >> {output};
        fi
        grep -v '#' {input} >> {output}
        """

def all_daughter_calls_function(wildcards):
    ## Function to take in wildcards (ie. the sample) and output a list of all the daughter calls as renamed above
    ##
    daughters = genomicsapi.select.multi_select(db = db, table = "samples", filters = {"sample_parent_id":genomicsapi.translate.stringtoid(wildcards.sample)}, bench = config["bench"])
    paths = []
    for daughter in daughters:
        ## id
        daught_id = daughter[0]
        ## Run gateway to ensure the analysis dir is created
        gateway("MUTECT", genomicsapi.translate.idtostring(daught_id, "MPS"), SCRATCH_DIR, config["PROD_DIR"], db = config["db"])
        ## Now we need the relevant analysis ID for the daughter sample
        analysis = genomicsapi.select.multi_select(db = db, table = "analyses", filters = {"pipeline_name":"MUTECT_CELLLINE", "samples_id":daught_id}, bench = config["bench"])
        p = analysis[0][10]
        chromstring = "renamed_" + str(wildcards.chrom) + ".vcf"
        vcf_path = os.path.join(p, genomicsapi.translate.idtostring(analysis[0][0], "MPA"), "mutect2", wildcards.reference, "germ", chromstring)
        paths.append(vcf_path)
    return(paths)


rule COMBINE_MUTECT2_CELLLINE_DAUGHTER_CALLS:
    input:
        all_daughter_calls_function
    output:
        f = "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/daughters_merged_vcf/{chrom}.vcf"
    singularity: "docker://lculibrk/gatk_alignment:latest"
    resources:
        threads = 1,
        mem_mb = 5000,
        runtime = 15,
        slurm_partition = config["partitions"]["short"],
    log:
        "data/studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/daughters_merged_vcf/{chrom}.log"
    params:
        inputlist = lambda wildcards, input: " -I ".join([input]) if isinstance(input, str) else " -I ".join(input)
    shell:
        """
        gatk MergeVcfs -I {params.inputlist} -O {output.f} &> {log}
        """