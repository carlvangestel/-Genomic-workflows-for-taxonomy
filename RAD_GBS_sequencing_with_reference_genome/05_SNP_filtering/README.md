In the final step, we generate a highly accurate VCF that retains only biallelic SNPs with high-quality genotypes and no missing data across all individuals. This kind of stringent filtering is useful when performing analyses that are sensitive to genotyping error or missingness. The filtering involves two stages: First, we use bcftools view to keep only SNPs that are biallelic, using the flags -m2 and -M2 to specify a minimum and maximum of two alleles per site. Next, we pipe the result into bcftools filter, applying several quality thresholds:

    A minimum genotype quality (GQ) of 30, which helps ensure confidence in the called genotypes.
    A minor allele frequency (MAF) filter to retain variants with frequencies between 0.01 and 0.99, removing both extremely rare and fixed variants.
    A requirement that there are no missing genotypes across samples (F_MISSING=0).

Since this filtering is applied only to SNPs, we can use the VCF file project.snps.vcf.gz generated in the previous step as input.

#!/bin/bash

cd ~/project

# Load modules
module load BCFtools

# Make VCF with stringent filtering
VCF_SNPS="./vcf/project.snps.vcf.gz"
VCF_STRINGENT="./vcf/project.stringent.vcf.gz"
bcftools view -m2 -M2 -v snps "$VCF_SNPS" | bcftools filter -i 'MIN(GQ)>=30 && MAF>=0.01 && MAF<=0.99 && F_MISSING=0' -Oz -o "$VCF_STRINGENT"
tabix -p vcf "$VCF_STRINGENT"
