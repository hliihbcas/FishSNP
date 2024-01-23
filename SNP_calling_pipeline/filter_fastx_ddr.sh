#!/usr/bin/env bash
### input file, a list file of SRA accessions.
PROJECT=$1 ##eg.cf_rnaseq, gc_bsseq
SRR=$2
wkdir=`pwd`

#mkdir ${wkdir}/00_data
#mkdir ${wkdir}/00_data/${PROJECT}

### Filter data
cd  ${wkdir}/00_data/${PROJECT}
N=`ls ${wkdir}/00_data/${PROJECT}/${SRR}*fastq.gz |wc -l ` 
if [[ $N -eq 2 ]];then
  zcat ${wkdir}/00_data/${PROJECT}/${SRR}_?.fastq.gz   |fastq_quality_filter -i - -o ${wkdir}/00_data/${PROJECT}/${SRR}_fastx.gz -q 20 -p 70 -z -Q 33
  zcat ${wkdir}/00_data/${PROJECT}/${SRR}_fastx.gz |grep  '^@DRR' |cut -d' ' -f1 |sort |uniq -c |awk  '{if($1==2){print $2 >> "'${wkdir}'/00_data/'${PROJECT}'/'${SRR}'_highQpairedreads.list"}else{print $2 >> "'${wkdir}'/00_data/'${PROJECT}'/'${SRR}'_highQsinglereads.list"}}'
 ## sed -i 's/^@//g'  ${wkdir}/00_data/${PROJECT}/${SRR}_highQpairedreads.list
 ## sed -i 's/^@//g'  ${wkdir}/00_data/${PROJECT}/${SRR}_highQsinglereads.list
 ## zcat ${wkdir}/00_data/${PROJECT}/${SRR}_1.fastq.gz   |perl ~/users/wtzhang/script/List_Fasta.pl ${wkdir}/00_data/${PROJECT}/${SRR}_highQpairedreads.list - QF |gzip - > ${wkdir}/00_data/${PROJECT}/${SRR}_1_fastx.gz  &
 ## zcat ${wkdir}/00_data/${PROJECT}/${SRR}_2.fastq.gz   |perl ~/users/wtzhang/script/List_Fasta.pl ${wkdir}/00_data/${PROJECT}/${SRR}_highQpairedreads.list - QF |gzip - > ${wkdir}/00_data/${PROJECT}/${SRR}_2_fastx.gz  &
 ## zcat ${wkdir}/00_data/${PROJECT}/${SRR}_1.fastq.gz ${wkdir}/00_data/${PROJECT}/${SRR}_2.fastq.gz |perl ~/users/wtzhang/script/List_Fasta.pl ${wkdir}/00_data/${PROJECT}/${SRR}_highQsinglereads.list - QF |gzip - > ${wkdir}/00_data/${PROJECT}/${SRR}_se_fastx.gz  &
  zcat ${wkdir}/00_data/${PROJECT}/${SRR}_1.fastq.gz   |grep -A 3 -wFf ${wkdir}/00_data/${PROJECT}/${SRR}_highQpairedreads.list - |sed '/^--$/d' |gzip - > ${wkdir}/00_data/${PROJECT}/${SRR}_1_fastx.gz  &
  zcat ${wkdir}/00_data/${PROJECT}/${SRR}_2.fastq.gz   |grep -A 3 -wFf ${wkdir}/00_data/${PROJECT}/${SRR}_highQpairedreads.list - |sed '/^--$/d' |gzip - > ${wkdir}/00_data/${PROJECT}/${SRR}_2_fastx.gz  &
  #zcat ${wkdir}/00_data/${PROJECT}/${SRR}_1.fastq.gz ${wkdir}/00_data/${PROJECT}/${SRR}_2.fastq.gz |grep -A 3 -wFf ${wkdir}/00_data/${PROJECT}/${SRR}_highQsinglereads.list - |sed '/^--$/d' |gzip - > ${wkdir}/00_data/${PROJECT}/${SRR}_se_fastx.gz  &
  wait
  echo "good" 
  #rm ${wkdir}/00_data/${PROJECT}/${SRR}_fastx.gz  # ${wkdir}/00_data/${PROJECT}/${SRR}_highQpairedreads.list ${wkdir}/00_data/${PROJECT}/${SRR}_highQsinglereads.list &
elif [[ $N -eq 1 ]];then
  zcat ${wkdir}/00_data/${PROJECT}/${SRR}*.gz    | fastq_quality_filter -i - -o ${wkdir}/00_data/${PROJECT}/${SRR}_fastx.gz -q 20 -p 70 -z -Q 33
fi 
cd ${wkdir}
## Remove raw fastq data
FILTERED_SIZE=`du ${wkdir}/00_data/${PROJECT}/${SRR}*_fastx.gz | awk '{print $1}' |head -1`
if [[ ${FILTERED_SIZE} -ge 10720 ]] ;then
  #rm ${wkdir}/00_data/${PROJECT}/${SRR}.sra  &
  rm ${wkdir}/00_data/${PROJECT}/${SRR}*fastq.gz  &
  rm ${wkdir}/00_data/${PROJECT}/${SRR}_highQpairedreads.list &
  rm ${wkdir}/00_data/${PROJECT}/${SRR}_highQsinglereads.list &
fi
