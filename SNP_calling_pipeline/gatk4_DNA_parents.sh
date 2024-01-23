#!/usr/bin/env bash
#####   mapping and GATK4
#### USAGE: sh gatk4_DNA_parents.sh [sample][fq][genome]
### sample
PROJECT=$1
SAMPLE=$2
### data
#GENOME=$4
GENOME=$3
DATA_TYPE=$4
#N=${N_GC3}
#DBSNP=
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
N_threads=10
source activate py2   ### py3 to py2
### index
## faidx
#samtools faidx ${GENOME}  ##for bowtie. bwa index  #
## dict     only a time
#${gatk4} --java-options "-Xmx10G" CreateSequenceDictionary -R ${GENOME} -O ${GENOME%.*}.dict

######################################################################################################################################################
####################  Data pre-processing for variant discovery
### mapping bwa
#bwa index ${GENOME}
#bwa mem -t 32 -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:WES\tPL:Illumina" ${INDEX} ${FQ}/${SAMPLE}.fq1.gz ${FQ}/${SAMPLE}.fq2.gz >${OUT}/${SAMPLE}.sam

### mapping bowtie2
#bowtie2-build ${GENOME} ${GENOME%.*} --threads=50  #
## bowtie2 -a -f bowtie2 ${FASTA} -p 5 > ${OUT}/${SAMPLE}.sam
##bowtie2  -x ${GENOME%/*}/bowtie2 -1 ${FQ1} -2 ${FQ2} -S ${OUT}/bowtie2.sam -p 6   #
N=`ls ${wkdir}/00_data/$PROJECT/${SAMPLE}*fastx.gz |wc -l`

if [ $N -eq 1 ]
then bowtie2  -x ${GENOME%.*} -U `ls ${wkdir}/00_data/$PROJECT/${SAMPLE}_fastx.gz |sed ":a;N;s/\n/,/g;ta"` -S ${OUT}/bowtie2.sam -p ${N_threads}
else bowtie2  -x ${GENOME%.*} -1 `ls ${wkdir}/00_data/$PROJECT/${SAMPLE}_1_fastx.gz |sed ":a;N;s/\n/,/g;ta"` -2 `ls ${wkdir}/00_data/$PROJECT/${SAMPLE}_2_fastx.gz |sed ":a;N;s/\n/,/g;ta"` -S ${OUT}/bowtie2.sam -p ${N_threads}
fi
### sort
${gatk4} --java-options "-Xmx20G" SortSam -SO coordinate -I ${OUT}/bowtie2.sam -O ${OUT}/bowtie2.bam  --TMP_DIR=$TEMP #-Djava.io.tmpdir=./    #
#--TMP_DIR
FILE_SIZE=`du ${OUT}/bowtie2.bam | awk '{print $1}'`
if [ $FILE_SIZE -ge 102400 ]
then
  rm ${OUT}/bowtie2.sam &
  rm -r ${wkdir}/00_data/$PROJECT/${SAMPLE}_*
fi
### samtools index
samtools index ${OUT}/bowtie2.bam   #

### alignment
#samtools flagstat ${OUT}/${SAMPLE}.bam > ${OUT}/${SAMPLE}.alignment.flagstat
#samtools stats ${OUT}/${SAMPLE}.bam > ${OUT}/${SAMPLE}.alignment.stat
#echo plot-bamstats -p ${SAMPLE}_QC ${OUT}/${SAMPLE}.alignment.stat

##add groups
##GATK2.0 and later will no longer support the mutation detection of headless files. This step can be performed during BWA alignment, and can be accomplished by selecting the -r parameter. If the -r parameter is not selected during BWA alignment, you can add this step.
${gatk4} --java-options "-Xmx20G" AddOrReplaceReadGroups -I ${OUT}/bowtie2.bam -O ${OUT}/addGroup.bam -ID ${SAMPLE} -LB reseq -PL illumine -PU hiseq -SM ${SAMPLE} -SO coordinate #VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true RGPL=illumina RGID=NoThornID RGSM=NoThorn RGLB=NoThornLB RGPU=NoThornPU 
##${picard} AddOrReplaceReadGroups I=${OUT}/${SAMPLE}.bam O=${OUT}/${SAMPLE}_addGroup.bam ID=${SAMPLE} LB=reseq PL=illumine PU=hiseq SM=${SAMPLE} SORT_ORDER=coordinate
##ID str: Enters the reads set ID. LB: read set library name; PL: Sequencing platform (illunima or solid); PU: name of subordinate unit of sequencing platform (name of run); SM: Sample name.

if [ ${DATA_TYPE} = "rad-seq" ]
then
### MarkDuplicates
#${gatk4} --java-options "-Xmx20G" MarkDuplicates -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 -I ${OUT}/addGroup.bam  -O ${OUT}/marked.bam -M ${OUT}/${SAMPLE}.metrics    #
### FixMateInformation
echo "i am rad-seq and i no MarkDuplicates steps" 
${gatk4} --java-options "-Xmx20G"  FixMateInformation -I ${OUT}/addGroup.bam -O ${OUT}/marked_fixed.bam -SO coordinate  --TMP_DIR=$TEMP#
else
    echo "my name is de nove data"
    ${gatk4} --java-options "-Xmx20G" MarkDuplicates -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 -I ${OUT}/addGroup.bam  -O ${OUT}/marked.bam -M ${OUT}/${SAMPLE}.metrics
    ${gatk4} --java-options "-Xmx20G"  FixMateInformation -I ${OUT}/marked.bam -O ${OUT}/marked_fixed.bam -SO coordinate --TMP_DIR=$TEMP
fi

FILE_SIZE=`du ${OUT}/marked_fixed.bam | awk '{print $1}'`
if [ $FILE_SIZE -ge 102400 ]
then
 rm ${OUT}/bowtie2.bam ${OUT}/addGroup.bam ${OUT}/marked.bam ${OUT}/bowtie2.bam.bai &
fi

### samtools index 2
samtools index ${OUT}/marked_fixed.bam   #

### bedtools genomecov
#bedtools genomecov -ibam ${OUT}/marked_fixed.bam -bga >${OUT}/bedtools.genomecov &

### HaplotypeCaller
${gatk4} --java-options "-Xmx20G"  HaplotypeCaller -R ${GENOME} -I ${OUT}/marked_fixed.bam -O ${OUT}/raw.g.vcf -ERC GVCF  -stand-call-conf 30 #--bam-output ${OUT}/HaplotypeCaller.bam # --sample-ploidy 30 --max-mnp-distance 0 #-L ${PARENTS_VCF}
## ${gatk4} --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller -R ${GENOME} -I ${OUT}/${SAMPLE}_marked_fixed.bam --dbsnp ${DBSNP} -O ${OUT}/${SAMPLE}_raw.db.vcf

rm -r $TEMP &

### coverage
#awk '{if($4==0)print}' ${OUT}/bedtools.genomecov |awk -v T=$N  -v S=$i '{a=$3-$2;sum=sum+a}END{print S"\t"(T-sum)/T}' >> ${OUT_}/all.coverage&
