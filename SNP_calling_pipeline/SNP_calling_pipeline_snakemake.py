##########################################################
# 2024/1/9
# Heng Li
# m17715629787@163.com
# SNP calling pipeline with snakemake
# snakemake -s SNP_calling_pipeline_snakemake.py --cores 1
##########################################################

example_project="example_data"

#Please note that this must be replaced with the absolute path of your "gatk" execution software. try(linux):which -a gatk
GATK_PATH="/home/wtzhang/soft/gatk-4.1.0.0/gatk"

rule all:
    input:
        expand("out/{example_project}/out.txt", example_project=example_project)

rule donwload:
    input:
        "00_data/{example_project}/sra.txt"
    output:
        "out/{example_project}/out.txt"
    params:
        project="{example_project}",
        species="tilapia",
        source="EBI",
        WGETF="kong",
        PREPATH="kong",
        DATA_TYPE="WGS",
        CONTROL_FLOW="parent",
        CHILDREN_INDEX="children",
        CHILDREN_LIST="kong",
        PARENTS_INDEX="parent",
        SPECIESINDEX="species/"
    #conda:
    #    "./FishSNP_env.yaml"
    shell:
        #"./SNP_calling_pipeline/NCBIdataProcessingGroup.sh {params.project} \
        #{input} tilapia {input} parent kong children parent WGS kong EBI kong> {log} 2>&1"
        "./NCBIdataProcessingGroup.sh {params.project} \
        {input} {params.species} {input} {params.PARENTS_INDEX} {params.CHILDREN_LIST} \
        {params.CHILDREN_INDEX} {params.CONTROL_FLOW} {params.DATA_TYPE} {params.PREPATH} \
        {params.source} {params.WGETF} {params.SPECIESINDEX} {GATK_PATH} >{output} 2>&1"
