#!/usr/bin/env bash
#####   mapping and GATK4
#### USAGE: sh gatk4_DNA_groups.sh [SAMPLE_LIST][PARENTS_CHILDREN][GENOME]
#### eg. sh gatk4_DNA_groups.sh parents.list parents children.list children
### sample
#SAMPLE=$1
### data
#FQ1=$2
#FQ2=$3
#PARENTS_VCF=$4

### sample.list
PARENTS_LIST=$1   ## parents.list   / children.list
PARENTS_INDEX=$2 ## eg. parents
CHILDREN_LIST=$3
CHILDREN_INDEX=$4  ## eg.children
PROJECT=$5
GENOME=$6
CONTROL=$7
DATA_TYPE=$8
gatk4=$9

PERL5LIB=""
echo $PERL5LIB
echo "i am here"
#GENOME=$5
#GENOME=/home/Xiaxq/users/wtzhang/201904_reseq/02_SNPcalling_87/pilon2.fasta
#GENOME=/home/Xiaxq/users/wtzhang/201904_finless/finless_eel/Denovo_genome/genome.fasta
#GENOME=${GENOME3}
#DBSNP=
wkdir=`pwd`
### snp calling for parents, *.g.vcf
if [ $CONTROL = "parent" ]
then
for i in `cat ${PARENTS_LIST}`
do
  echo "i am here inside"
  #nomal pipeline
  #sh gatk4_DNA_parents.sh ${PROJECT} ${i} ${GENOME} ${DATA_TYPE} &

  #pipeline without saving bam file
  sh gatk4_DNA_parents_savebam.sh ${PROJECT} ${i} ${GENOME} ${DATA_TYPE} ${gatk4} &
done

wait
echo "parent gatk over"

#wait
### merge parents snp, *.vcf
### filter, *.pass.vcf
sh gatk4_DNA_merge_filter.sh ${PARENTS_LIST} ${PARENTS_INDEX} ${GENOME} ${PROJECT} ${gatk4}

### snp calling for children, *.g.vcf
else
for j in `cat ${CHILDREN_LIST}` # `head -87 ${CHILDREN_LIST} |tail -24`
do
   sh gatk4_DNA_children.sh ${PROJECT} ${j} ${GENOME}  ${wkdir}/01_SNPcalling/${PROJECT}/multiple/${PARENTS_INDEX}_snp.pass.vcf ${DATA_TYPE} ${gatk4} &
done
wait
echo "child gatk over"
### 
sh gatk4_DNA_merge_filter.sh ${CHILDREN_LIST} ${CHILDREN_INDEX} ${GENOME} ${PROJECT} ${gatk4}
fi
