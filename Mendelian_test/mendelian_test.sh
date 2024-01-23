vcf=$1
samplelist=$2
PROJECT=$3
PVALUE=$4
##1. Sample sort and filter
python3 ./vcf_sort_filter.py ${vcf} ${vcf%%.*}_sort.vcf -f ${samplelist} -a 
echo "a"
##2.
#./vcf2genotype.py -i ${vcf%%.*}_sort.vcf -o ${vcf%%.*}_sort.geno
python3 ./vcf2genotype_cdy.py ${vcf%%.*}_sort.vcf ${vcf%%.*}_sort.geno
echo "b"
##3.
python3 ./hgmap_treat.py ${vcf%%.*}_sort.geno -o ${vcf%%.*}_sort
echo "c"
##4.
python3 ./hgmap_merge.py ${vcf%%.*}_sort.raw -o ${vcf%%.*}_sort
echo "b"
##5.
python3 ./hgmap_chiqtest.py -c ${PVALUE} ${vcf%%.*}_sort_subBin.wraw -o ${vcf%%.*}_sort_subBin 
echo "e"
