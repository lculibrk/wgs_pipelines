## Probably this should all go in the database eventually but for now this is easier for prototyping

db_deps = {
    "SOMATIC": [
        "studies",
        "samples"
    ],
    "FASTQ": [
        "studies",
        "samples",
        "runs"
    ],
    "WGS_MERGE_BAM": [
        "studies",
        "samples"
    ],
    "MUTECT": [
        "studies",
        "samples"
    ],
    "EXTERNAL_BAM":[
        "studies",
        "samples"
    ],
    "LOAD_EXTERNAL_BAM":[
        "studies",
        "samples"
    ],
    "MUTECT_CELLLINE":[
        "studies",
        "samples"
    ],
    "INDEL":[
        "studies",
        "samples"
    ]

}

## Describes the endpoint of a pipeline
## Controlled vocabulary: [FASTQ, MERGE_BAM]
module_outputs = {
    "SOMATIC":"SOMATIC",
    "TESTS":"FASTQ",
    "SRA":"FASTQ",
    "EGA":"FASTQ",
    "GATK_BAM":"WGS_MERGE_BAM",
    "MUTECT_CELLLINE":"MUTECT",
    "LOAD_EXTERNAL_BAM":"LOAD_EXTERNAL_BAM",
    "EXTERNAL_BAM":"EXTERNAL_BAM",
    "INDEL":"INDEL",
}

## Invert the dict

filetypes = {}
for key, value in module_outputs.items():
    if value not in filetypes:
        filetypes[value] = [key]
    else:
        filetypes[value].append(key)

## Empty list if the module requires no input from a module
## Controlled vocabulary: [FASTQ, WGS_MERGE_BAM]
module_inputs = {
    "TESTS":[],
    "SRA":[],
    "EGA":[],
    "GATK_BAM":["FASTQ"],
    "MUTECT_CELLLINE":["WGS_MERGE_BAM"],
    "EXTERNAL_BAM":[],
    "LOAD_EXTERNAL_BAM":[],
    "INDEL":["WGS_MERGE_BAM"]
}

## e.g. we can see that GATK_BAM requires a FASTQ module as input and outputs a WGS_MERGE_BAM