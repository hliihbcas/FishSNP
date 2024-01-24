#!/usr/bin/env bash
#####   mapping and GATK4
#### USAGE: sh gatk4_DNA_children.sh [sample][fq][genome]

### sample
PROJECT=$1
SAMPLE=$2
### data
GENOME=$3
PARENTS_VCF=$5

DATA_TYPE=$4
gatk4=$6
#GENOME=$5

#N=${N_GC3}
#DBSNP=

N_threads=10

## outfile
wkdir=`pwd`
mkdir ${wkdir}/01_SNPcalling
mkdir ${wkdir}/01_SNPcalling/$PROJECT
mkdir ${wkdir}/01_SNPcalling/$PROJECT/$SAMPLE
OUT=${wkdir}/01_SNPcalling/$PROJECT/$SAMPLE

mkdir ${wkdir}/temp
mkdir  `pwd`/temp/$PROJECT/
mkdir `pwd`/temp/$PROJECT/tmp_${SAMPLE}
TEMP=`pwd`/temp/$PROJECT/tmp_${SAMPLE}


#source activate py2   ### py3 to py2

######################################################################################################################################################
####################  Data pre-processing for variant discovery
##################################################################################################
### mapping bwa
##bwa mem -t 32 -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:WES\tPL:Illumina" ${INDEX} ${FQ}/${SAMPLE}.fq1.gz ${FQ}/${SAMPLE}.fq2.gz >${OUT}/${SAMPLE}.sam
##################################################################################################

### mapping bowtie2
#bowtie2-build ${GENOME} ${GENOME%.*}  --threads 60 #
N=`ls ${wkdir}/00_data/PRJNA561042_2/${SAMPLE}*fastx.gz |wc -l`

if [ $N -eq 1 ]
then bowtie2  -x ${GENOME%.*} -U `ls ${wkdir}/00_data/PRJNA561042_2/${SAMPLE}_fastx.gz |sed ":a;N;s/\n/,/g;ta"` -S ${OUT}/bowtie2.sam -p ${N_threads}
else bowtie2  -x ${GENOME%.*} -1 `ls ${wkdir}/00_data/PRJNA561042_2/${SAMPLE}_1_fastx.gz |sed ":a;N;s/\n/,/g;ta"` -2 `ls ${wkdir}/00_data/PRJNA561042_2/${SAMPLE}_2_fastx.gz |sed ":a;N;s/\n/,/g;ta"` -S ${OUT}/bowtie2.sam -p ${N_threads}
fi


### sort
${gatk4} --java-options "-Xmx20G" SortSam -SO coordinate -I ${OUT}/bowtie2.sam -O ${OUT}/bowtie2.bam --TMP_DIR=$TEMP --VALIDATION_STRINGENCY SILENT
#-Djava.io.tmpdir=./    #

FILE_SIZE=`du ${OUT}/bowtie2.bam | awk '{print $1}'`
if [ $FILE_SIZE -ge 10240 ]
then
  rm ${OUT}/bowtie2.sam &
  rm -r ${wkdir}/00_data/$PROJECT/${SAMPLE}_*
fi

### samtools index
samtools index ${OUT}/bowtie2.bam   #

${gatk4} --java-options "-Xmx20G" AddOrReplaceReadGroups -I ${OUT}/bowtie2.bam -O ${OUT}/addGroup.bam -ID ${SAMPLE} -LB reseq -PL illumine -PU hiseq -SM ${SAMPLE} -SO coordinate #VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true RGPL=illumina RGID=NoThornID RGSM=NoThorn RGLB=NoThornLB RGPU=NoThornPU 
##${picard} AddOrReplaceReadGroups I=${OUT}/${SAMPLE}.bam O=${OUT}/${SAMPLE}_addGroup.bam ID=${SAMPLE} LB=reseq PL=illumine PU=hiseq SM=${SAMPLE} SORT_ORDER=coordinate
if [ ${DATA_TYPE} = "rad-seq" ]
then
echo "i am rad-seq and i no MarkDuplicates steps"
${gatk4} --java-options "-Xmx20G"  FixMateInformation -I ${OUT}/addGroup.bam -O ${OUT}/marked_fixed.bam -SO coordinate  --TMP_DIR=$TEMP#
else
    echo "my name is de nove data"
    ${gatk4} --java-options "-Xmx20G" MarkDuplicates -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 -I ${OUT}/addGroup.bam  -O ${OUT}/marked.bam -M ${OUT}/${SAMPLE}.metrics
    ${gatk4} --java-options "-Xmx20G"  FixMateInformation -I ${OUT}/marked.bam -O ${OUT}/marked_fixed.bam -SO coordinate --TMP_DIR=$TEMP
fi

### MarkDuplicates
#${gatk4} --java-options "-Xmx10G" MarkDuplicates -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 -I ${OUT}/addGroup.bam  -O ${OUT}/marked.bam -M ${OUT}/${SAMPLE}.metrics    #
### FixMateInformation
#${gatk4} --java-options "-Xmx10G"  FixMateInformation -I ${OUT}/marked.bam -O ${OUT}/marked_fixed.bam -SO coordinate  --TMP_DIR=$TEMP #

FILE_SIZE=`du ${OUT}/marked_fixed.bam | awk '{print $1}'`
if [ $FILE_SIZE -ge 10240 ]
then
 rm ${OUT}/bowtie2.bam ${OUT}/addGroup.bam ${OUT}/marked.bam ${OUT}/bowtie2.bam.bai &
fi

### samtools index 2
samtools index ${OUT}/marked_fixed.bam   #

### bedtools genomecov
#bedtools genomecov -ibam ${OUT}/marked_fixed.bam -bga >${OUT}/bedtools.genomecov &

### HaplotypeCaller
${gatk4} --java-options "-Xmx20G"  HaplotypeCaller -R ${GENOME} -I ${OUT}/marked_fixed.bam -O ${OUT}/raw.g.vcf -ERC GVCF -L ${PARENTS_VCF} -stand-call-conf 30 #--bam-output ${OUT}/HaplotypeCaller.bam #--sample-ploidy 2 --max-mnp-distance 0 #-L ${PARENTS_VCF}

rm -r $TEMP &

### coverage
#awk '{if($4==0)print}' ${OUT}/bedtools.genomecov |awk -v T=$N  -v S=$i '{a=$3-$2;sum=sum+a}END{print S"\t"(T-sum)/T}' >> ${OUT_}/all.coverage &
