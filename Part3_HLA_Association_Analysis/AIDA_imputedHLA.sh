imputed="/mnt/icbs_shared_storage_general_2022/Benjamaporn.sri/AIDA/HLA_2026/HLA/chr6.dose"

plink2 \
 --vcf ${imputed}.vcf.gz dosage=DS \
 --make-pgen --out ${imputed}

# this will create ${imputed}.{pgen,psam,pvar}
# when you only test for HLA alleles and amino acids (modify as necessry)
cat ${imputed}.pvar | grep -v "#" | grep -E 'HLA_|AA_' | cut -f3 >  test_markers.txt
cat ${imputed}.pvar | grep -v "#" | grep -E 'HLA_|AA_' | awk '{print $3,"T"}' >  test_markers_alleles.txt # this is to make sure to output the dosages of the "presence" of the allele coded as T, but not the "absence" coded as A.

plink2 --vcf ${imputed}.vcf.gz dosage=DS --export A --extract test_markers.txt --export-allele test_markers_alleles.txt --out ${imputed}
