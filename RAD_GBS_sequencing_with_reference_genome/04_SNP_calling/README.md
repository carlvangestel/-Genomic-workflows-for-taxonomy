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

# SNP calling
bcftools mpileup --min-MQ 30 -a AD,DP,SP -Ou -f "$GENOME" sample01.rmd.bam sample02.rmd.bam sample03.rmd.bam ... | bcftools call -f GQ,GP -m -Oz -o "$VCF_RAW"
```

The first part of the command, `bcftools mpileup`, generates the allelic count data for each individual. The `-f` flag specifies the reference genome, while `-a AD,DP,SP` adds annotations for allele depths (`AD`), total read depth (`DP`), and strand bias (`SP`) to the output. We also want to avoid SNPs to be called from reads that map ambiguously and filtering them out by `–min-MQ 30` (minimum read mapping quality > 30). The `-Ou` flag outputs an uncompressed VCF to standard output. The BAM files are listed at the end.
This output is then piped into `bcftools call`, which performs the actual variant calling. The `-f GQ,GP` flag includes genotype quality (`GQ`) and genotype probabilities (`GP`) in the output VCF. The `-m` flag invokes the multiallelic caller, and `-Oz` specifies compressed VCF (=BCF) output (`vcf.gz`), which is recommended as these VCF files can be very large as they include information of every genomic position. Finally, `-o` defines name of the output VCF. To avoid confusion, it is highly recommended that the output file name extension is `.vcf.gz` if the output vcf is compressed (`-Oz`). If you want the output file to be uncompressed (`-Ou`), use `.vcf` as output file extension. 
We finally index the compressed VCF for fast querying in downstream applications with the **tabix** tool that is included in **bcftools**. 

```bash
tabix -p vcf "$VCF_RAW"
```

