#!/usr/bin/env bash
#####   mapping and GATK4
#### USAGE: sh gatk4_DNA_merge_filter.sh [SAMPLE_LIST][PARENTS_CHILDREN][GENOME]

### sample.list    
SAMPLE_LIST=$1   ## parents.list   / children.list
PARENTS_CHILDREN=$2 ## parents or children  or multiple
GENOME=$3
PROJECT=$4
gatk4=$5
## outfile
wkdir=`pwd`
OUT_=${wkdir}/01_SNPcalling/$PROJECT

#source activate py2   ### py3 to py2

### merge gvcf
mkdir ${wkdir}/01_SNPcalling/$PROJECT/multiple
MUL=${wkdir}/01_SNPcalling/$PROJECT/multiple
awk -v A=${OUT_} '{print A"/"$0}' ${SAMPLE_LIST} | sed 's/$/\/raw\.g\.vcf/' |sed 's/^/-V /' |echo `cat -` > ${MUL}/${PARENTS_CHILDREN}_raw_gvcf.list
${gatk4} --java-options "-Xmx100G" CombineGVCFs -R ${GENOME} `cat ${MUL}/${PARENTS_CHILDREN}_raw_gvcf.list` -O ${MUL}/${PARENTS_CHILDREN}_raw.g.vcf  #-XL ${MUL}/mnp.intervals   # -Xms20G -Xmn5G -Xss1m -XX:PermSize=256m -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit
${gatk4} --java-options "-Xmx100G" GenotypeGVCFs -R ${GENOME} -V ${MUL}/${PARENTS_CHILDREN}_raw.g.vcf -O ${MUL}/${PARENTS_CHILDREN}_raw.vcf  #[--dbsnp ${DBSNP}]

##SelectVariants
${gatk4} --java-options "-Xmx10G" SelectVariants -select-type SNP -R ${GENOME} -V ${MUL}/${PARENTS_CHILDREN}_raw.vcf  -O ${MUL}/${PARENTS_CHILDREN}_snp.vcf
#${gatk4} --java-options "-Xmx10G" SelectVariants -select-type INDEL -R ${GENOME} -V ${MUL}/${PARENTS_CHILDREN}_raw.vcf  -O ${MUL}/${PARENTS_CHILDREN}_indel.vcf



##Variant filter:hard filter(VariantFiltration) and VQSR
${gatk4} --java-options "-Xmx10G" VariantFiltration  -R ${GENOME} -V ${MUL}/${PARENTS_CHILDREN}_snp.vcf -O ${MUL}/${PARENTS_CHILDREN}_snp.filter.vcf --filter-expression "QD<2.0 || MQ<35.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0" --filter-name LowQualFilter --missing-values-evaluate-as-failing TRUE --verbosity ERROR
grep "^#" ${MUL}/${PARENTS_CHILDREN}_snp.filter.vcf > ${MUL}/vcf.head
grep -v "^#" ${MUL}/${PARENTS_CHILDREN}_snp.filter.vcf |awk '{if($7~/PASS/)print}' > ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf.
cat ${MUL}/vcf.head ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf. > ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf
rm ${MUL}/vcf.head ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf.

##VariantEval
#${GATK} -R ${GENOME} -T VariantEval --eval ${SAMPLE}_06.snp.filter.vcf

### merge parents and children. vcf-merge
#bgzip ${PARENTS_VCF} > ${PARENTS_VCF}.gz
#bgzip ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf > ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf.gz
#tabix -p vcf ${PARENTS_VCF}.gz
#tabix -p vcf ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf.gz
#vcf-merge ${PARENTS_VCF}.gz ${MUL}/${PARENTS_CHILDREN}_snp.pass.vcf.gz > ${MUL}/parents_children.vcf

