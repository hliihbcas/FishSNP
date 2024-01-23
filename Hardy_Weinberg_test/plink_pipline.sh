Filepath=$1
Input_file=${Filepath}/multiple_snp.pass.vcf
plink  --vcf ${Input_file} --make-bed --out ${Filepath}/multiple_snp.pass --allow-extra-chr

#hardy test
plink --bfile ${Filepath}/multiple_snp.pass --hardy --allow-extra-chr

mkdir ./${Filepath}/plink_test

#filter
plink --bfile ${Filepath}/multiple_snp.pass --hwe 0.000001 --make-bed --out ./${Filepath}/plink_test/multiple_snp.pass --allow-extra-chr

grep  "##" ${Filepath}/multiple_snp.pass.vcf > ${Filepath}/multiple_snp.pass.vcf1
#body vcf
grep -v "##" ${Filepath}/multiple_snp.pass.vcf > ${Filepath}/multiple_snp.pass.vcf2
#paste ready file
#sed '1d' plink.hwe > plink.hwe2
awk '{print$7"\t"$8"\t"$9}' plink.hwe > ${Filepath}/plink.hwe2
#paste file : multiple_hardy_pass_snp.vcf
paste -d "\t" ${Filepath}/plink.hwe2 ${Filepath}/multiple_snp.pass.vcf2 > ${Filepath}/multiple_hardy_pass_snp.vcf
# : multiple_hardy_pass_snp.vcf2
grep -v "#" ${Filepath}/multiple_hardy_pass_snp.vcf > ${Filepath}/multiple_hardy_pass_snp.vcf2
grep  "#"  ${Filepath}/multiple_hardy_pass_snp.vcf > ${Filepath}/head_chrom.txt

#according to Pvalue,extract the pass snp : multiple_hardy_pass_snp.vcf3
awk '{if($3 >= 0.000001) print $0}' ${Filepath}/multiple_hardy_pass_snp.vcf2 > ${Filepath}/multiple_hardy_pass_snp.vcf3
#full file:multiple_hardy_pass_snp_final.vcf
cat ${Filepath}/multiple_snp.pass.vcf1 ${Filepath}/head_chrom.txt ${Filepath}/multiple_hardy_pass_snp.vcf3 > ${Filepath}/multiple_hardy_pass_snp_final.vcf
