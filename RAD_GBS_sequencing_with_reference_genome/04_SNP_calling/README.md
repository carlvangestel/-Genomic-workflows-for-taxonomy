## SNP calling

Using the BAM files generated in the previous step, we can now identify polymorphic genomic positions—known as single nucleotide polymorphisms (SNPs)—and infer each individual's genotype at those sites. Common tools for SNP calling include the **Genome Analysis Toolkit (GATK)** and **bcftools**. **GATK** is widely used for high-quality variant calling in human genomes, while **bcftools** is well-suited for non-model organisms with moderate-sized genomes (<3 Gb), offering comparable accuracy for these organisms (1).

SNP calling with bcftools is a two-step process:
1.	Counting reads supporting each allele at variable positions.
2.	Inferring genotypes using likelihood-based methods.

The output is a standardized **Variant Call Format (VCF)** file, which stores the detected variants and genotypes as well as information on the read depth and quality of the genotypes. 
Given a set of BAM files, we will use the following bcftools command to perform SNP calling and to generate a compressed VCF (=BCF) file:

```bash
#!/bin/bash

cd ~/project

# Load modules
module load BCFtools

# Define output VCF
VCF_RAW="./vcf/project.raw.vcf.gz"

