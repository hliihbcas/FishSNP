###################################
# 2024/1/9
# Lei Zhang & Heng Li
# 2211696720@qq.com/m17715629787@163.com
# Hardy_Weinberg_test
# snakemake -s hardy_pipeline_snakemake.py --cores 1
###################################

example_project="example_data"


rule all:
    input:
        expand("{example_project}/out.txt", example_project=example_project)

rule donwload:
    input:
        "{example_project}/multiple_snp.pass.vcf"
    output:
        "{example_project}/out.txt"
    params:
        project="{example_project}"
    #conda:
    #    "../FishSNP_01.yaml"
    shell:
        "bash plink_pipline.sh {params.project} >{output} 2>&1"
