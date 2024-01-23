# !/usr/bin/env bash
PROJECT=$1
SRRLIST=$2
SPECIES=$3

PARENTS_LIST=$4
PARENTS_INDEX=$5
CHILDREN_LIST=$6
CHILDREN_INDEX=$7

CONTROL=$8
DATA_TYPE=$9
PREPATH=${10}
DATABASE=${11}
WGETPATH=${12}
SPECIESINDEX=${13}
gatk4=${14}
INDEX=${SPECIESINDEX}${SPECIES}/
SPECIESFA=${INDEX}${SPECIES}.fa
wkdir=`pwd`
PERL5LIB=""
echo $PERL5LIB
echo $PREPATH

if [ ! -f ${INDEX}*.fai ];then
	samtools faidx $SPECIESFA
fi

if [ ! -f ${INDEX}*.dict ];then
	${gatk4} --java-options "-Xmx10G" CreateSequenceDictionary -R $SPECIESFA -O ${SPECIESFA%.*}.dict
fi

if [ ! -f ${SPECIESFA%.*}.1.bt2 ];then
	bowtie2-build $SPECIESFA ${SPECIESFA%.*} --threads=50
fi


#GENOME_INDEX=$3 ##eg.GENOME_cf, GENOME_sm, GENOME_gc, GENOME_nt
#GFF_INDEX=$4 ##eg.GFF_cf,...

#mkdir ${wkdir}/00_data
mkdir ${wkdir}/00_data/${PROJECT}

### download data
for i in `cat ${SRRLIST}`
do
#WGETPATH
    if [ $WGETPATH = "nochange" ];then	
    	echo ${PREPATH}/${i}/${i}.1
    	wget -c ${PREPATH}/${i}/${i}.1 -O ${wkdir}/00_data/${PROJECT}/${i}.sra &
    else
    	WGETURL=`srapath ${i}`
    	wget -c ${WGETURL} -O ${wkdir}/00_data/${PROJECT}/${i}.sra &
    fi
done

wait
sleep 10s

for i in `cat ${SRRLIST}`
do
vdb-validate ${wkdir}/00_data/${PROJECT}/${i}.sra 1>${wkdir}/00_data/${PROJECT}/${i}_check.log 2>${wkdir}/00_data/${PROJECT}/${i}_check.err &
done

wait
sleep 10s
for i in `cat ${SRRLIST}`
do
    fastq-dump --split-3 --gzip ${wkdir}/00_data/${PROJECT}/${i}.sra -O ${wkdir}/00_data/${PROJECT}/ &
done


wait
sleep 10s
for i in `cat ${SRRLIST}`
do
    rm ${wkdir}/00_data/${PROJECT}/${i}.sra &
done



wait
for i in `cat ${SRRLIST}`
do
  #sh ./download.sh ${PROJECT} ${i}
  if [ $DATABASE = "NCBI" ];then
	echo "NCBI"
  	sh ./filter_fastx.sh ${PROJECT} ${i} &
  else
     if [ $DATABASE = "EBI" ];then
	echo "EBI"
        sh ./filter_fastx_ebi.sh ${PROJECT} ${i} &
     else
        echo "DDR"
        sh ./filter_fastx_ddr.sh ${PROJECT} ${i} &
     fi
  fi
#  parent SNP calling
# sh ./gatk4_DNA_parents.sh ${PROJECT} ${i} ${SPECIESFA} &
done

wait
#  offspring SNP calling
./gatk4_DNA_groups.sh ${PARENTS_LIST} ${PARENTS_INDEX} ${CHILDREN_LIST} ${CHILDREN_INDEX} ${PROJECT} ${SPECIESFA} ${CONTROL} ${DATA_TYPE} ${gatk4}

wait

for i in `cat ${SRRLIST}`
do
#remenber to active
#rm -r ${wkdir}/00_data/$PROJECT/${i}_* &
done

wait


echo "ALL is down"

#Single end sequencing

#wait
#sh gatk4_DNA_merge_filter.sh ${PROJECT} ${SRRLIST} ${SPECIESFA} "allmerge" 1>stand.txt 2>err.txt
#sh gatk4_DNA_merge_filter.sh ${PROJECT} ${TEMP} ${SPECIESFA} "allmerge" 1>stand.txt 2>err.txt

### check fail 
#ls ${wkdir}/00_data/${PROJECT}/*.sra |sed 's/\.sra//g' |awk -F'/' '{print $NF}' > ${wkdir}/00_data/${PROJECT}/next_SRR.list
#NN=`less -S ${wkdir}/00_data/${PROJECT}/next_SRR.list |wc -l `
#while [[ $NN -gt 0 ]]
#do
#  for j in `cat ${wkdir}/00_data/${PROJECT}/next_SRR.list`
#  do
#    sh ./download.sh ${PROJECT} ${j}
#    sh ./filter_fastx.sh ${PROJECT} ${j} 
 #   sh ./salmon.sh ${PROJECT} ${i} $3 $4  &
#  done
#  ls ${wkdir}/00_data/${PROJECT}/*.sra |sed 's/\.sra//g' |awk -F'/' '{print $NF}' > ${wkdir}/00_data/${PROJECT}/next_SRR.list
#  NN=`less -S ${wkdir}/00_data/${PROJECT}/next_SRR.list |wc -l `
#done 
