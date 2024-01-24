############################################
# 2024/1/9
# Lei Zhang & Heng Li
# 2211696720@qq.com/m17715629787@163.com
# Mendelian_test
# snakemake -s mendelian_test_snakemake.py --cores 1
#############################################

example_project="example_data"


rule all:
    input:
        expand("{example_project}/out.txt", example_project=example_project)

rule donwload:
    input:
        "{example_project}/multiple_sort_vcf.vcf",
        "{example_project}/sample.txt"
    output:
        "{example_project}/out.txt"
    params:
        project="{example_project}",
        pvalue="0.05"
    #conda:
    #   "../FishSNP_01.yaml"
    shell:
        "bash mendelian_test.sh {input[0]} {input[1]} {params.project} {params.pvalue} >{output} 2>&1"
