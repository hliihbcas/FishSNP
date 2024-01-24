# FishSNP
The code repository of the FishSNP database.  
The script process, including the test data, can also be downloaded at the following address ().  
The script should be run in the Linux environment, and make sure that the conda in the environment can be used normally.  

## Import the conda environment
1.Before the process begins, you need to import the conda environment and run the following example in the active conda environment  
$ cd FishSNP-main  
$ conda env create -f FishSNP_01.yaml  
$ source activate FishSNP_02  

## SNP calling pipeline
The process includes downloading, quality control, SNP calling, etc.  

**Before you start the snakemake process, here are a few things you need to do**  
1.Install GATK 4.0 in advance.  
2.Install libcrypto.so.1.0.0, libbz2.so.1.0 in {Your lib path}/FishSNP_01/lib/  
3.Change the "GATK_PATH" variable in the script (./SNP_calling_pipeline/SNP_calling_pipeline_snakemake.py) to the absolute path of the GATK in your environment.  
4.Download the corresponding genome index file (tilapia) from the following link (), and extract it in the ./SNP_calling_pipeline/species/tilapia/.  

**Getting started (Linux)**  
$ cd ./SNP_calling_pipeline  
$ snakemake -s SNP_calling_pipeline_snakemake.py --cores 2

## Mendelian test
The process includes the Mendelian test.  

**Before you start the snakemake process, here are a few things you need to do**  
1.gunzip ./Mendelian_test/example_data/multiple_sort_vcf.vcf.gz

**Getting started (Linux)**   
$ cd Mendelian_test  
$ snakemake -s mendelian_test_snakemake.py --cores 1

## Hardy Weinberg test
The process includes the Hardy Weinberg test.  
**Getting started (Linux)**  
$ cd Hardy_Weinberg_test  
$ snakemake -s Hardy_Weinberg_test/hardy_pipeline_snakemake.py --cores 1

