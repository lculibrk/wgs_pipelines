# wgs_pipelines

Database, API, and pipelines for WGS analysis, written by [Luka Culibrk](https://github.com/lculibrk)

## genomicsDB

the schema is found under `db/schema/schema.sql`. To use the pipelines, you'll need to add `studies` entities to track projects that describe groups of samples, `samples` entities to track individual biological samples (ie. NCBI BioSamples), and `runs` to track sequencing runs (ie. multiple sequencing runs of a library)

To use the API, `pip install db/api/`

## Pipelines
```cd gatk_bestpractices_mapping```

### Quick Execute
```python ./runner.py [--idfile ID_FILE | --id SAMPLE_ID] --pipeline PIPELINE_NAME [snakemake options]```

### Overview

The pipelines are a collection of [Snakemake](https://snakemake.readthedocs.io/en/stable/) modules, each of which does a specific task in genomics, e.g. WGS alignment, or variant calling in a cell line parent-daughter context. 

Granted, many different workflows make, for example, a VCF. Which workflow is run, is determined from metadata that is stored on the DB under `db/`

The runner function, `runner.py` is given a sample (or list of samples) alongside an endpoint. For example, `WGS_MERGE_BAM` is the endpoint for making a WGS cram. 
The function determines based on the sample metadata on the genomicsdb which module should be used to generate the cram file - for example, it will determine whether raw data is on EGA or SRA, based on its corresponding entries on the `runs` table, and use the corresponding pipeline.
Similarly, for a parent-daughter cell line experiment, the pipeline determines the sample's matched parental sample as well as its sibling clones based on the database metadata. 
Refer to the database/API repo for information on how the data are stored there.

Current endpoints:
- FASTQ generation
- Alignment
- SNV calling with Mutect2

Each endpoint has at least one module that can generate that endpoint. 

FASTQ:
- SRA-hosted data (`SRA`)
- EGA-hosted data (`EGA`)
- Local BAM files to be remapped (`EXTERNAL_BAM`)

Alignment:
- GATK WGS best practices CRAM (`GATK_BAM`)

SNV calling:
- Parent-daughter cell line (`MUTECT_CELLLINE`)


## Developer documentation

### How the pipelines work

These pipelines are intended to be as simple to execute as possible. `gateway()`,  a convenience function, exists to enable the quick accession of the endpoint of an analysis for a given entity (study/sample/run). For example, a call to `gateway()` for the FASTQ endpoint determines the appropriate source for the FASTQ files based on their entries in the `genomicsdb`. Similarly, a `gateway()` call for the `MUTECT` endpoint should determine whether the samples are cell lines or patients, and thereby determine which module to use to satisfy the `MUTECT` endpoint.

Overall, the `genomicsdb` should store all the metadata required to handle how the pipelines are run. The pipeline framework should, therefore, utilize this information to reproducibly and consistently execute analyses on data entities in a generalizable fashion. 

### To add a new pipeline:
1. Create a new Snakefile module under `modules/`
2. Modify `modules/db_deps.py` to include this new pipeline:
   
   a) modify the `db_deps` dictionary. It expects a key:list pair, where the list elements are the DB tables needed for this analysis. For most analyses downstream of CRAM, this should be studies and samples.

   b) add the pipeline to the `module_outputs` dictionary. This maps the name of the pipeline to the type of analysis. For example, we would run a different Mutect pipeline for biopsies vs cell lines, so this dict tracks that both analyses would create a `MUTECT` endpoint. Similarly, we have multiple ways to make a FASTQ, depending on the data source, and this is documented accordingly in `module_outputs`.
   
   c) Modify `module_inputs`. This maps each module with the **endpoint** that is required for input. For example, the mapping module needs FASTQ. It doesn't particularly care where the FASTQ came from, just that it gets made.

3. Add it to `gateway()` in `lib/input_functions.py`. If it's a new endpoint, you'll need to add code to handle that endpoint. Otherwise, if it's a new module for an existing endpoint, you need to add the appropriate code to handle this module and decide when it should be executed (as opposed to another module to satisfy the endpoint). 




