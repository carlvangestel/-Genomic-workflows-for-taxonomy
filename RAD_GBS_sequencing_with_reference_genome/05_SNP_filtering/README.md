## SNP filtering

### 1. VCF of SNPs shared in 80% of individuals
### 1. VCF of SNPs shared in all indivudals
### 1. VCF of thinned SNPs

The VCF file generated in the previous step contains genotype information for all genomic positions, including both variable (polymorphic) and invariant sites, as well as genotype calls made with low quality. While this raw VCF (hence the proposed filename extension `.raw.vcf.gz` ) is useful as a starting point, this unfiltered VCF is generally not suitable for most downstream analyses, which typically require only high-confidence, polymorphic sites.
Filtering a VCF file is a critical step that can greatly influence analytical results, especially in population genomics, phylogenetics, or variant association studies. Care must be taken when applying filters to avoid excluding important variants or retaining low-quality data. In the following sections, we will generate three different VCF files from the raw input using **bcftools**, each filtered for a specific goal or application. 

### 1. VCF without indels

For applications such as generating aligned genomic consensus sequences (e.g., for phylogenetic inference), we need a VCF that retains all positions, including invariant sites, but excludes variants that disrupt sequence alignment — specifically insertions and deletions (indels) and multi-nucleotide polymorphisms (MNPs).
Indels are particularly problematic in this context because they alter the coordinate space of the alignment, potentially creating frameshifts or misalignments across samples. To avoid this, we remove indels and MNPs using the `-V` flag with **bcftools view**, which allows us to exclude specified variant types. Importantly, this filtering step retains invariant (non-polymorphic) sites, which are required to reconstruct the full consensus sequence for each indivdiual. We also set called genotypes that are based on a low sequencing depth (e.g. a sequencing depth of at least 5 = `DP<5`) as missing (`./.`) to avoid that read errors are called as SNPs. Setting the individual genotypes is done in the second part of the command through the `+setGT` plugin. 

```bash
#!/bin/bash

cd ~/project

# Load modules
module load BCFtools

# Make VCF without indels
VCF_RAW="./vcf/project.raw.vcf.gz"
VCF_NOINDEL="./vcf/project.noindel.minDP5.vcf.gz"
bcftools view -V indels,mnps "$VCF_RAW" -Ou | bcftools +setGT -Oz -o "$VCF_NOINDEL" -- -t q -n . -i 'FMT/DP<5'
tabix -p vcf "$VCF_NOINDEL"
```

### 2. VCF with SNPs only

The raw VCF is typically very large because it includes all genomic positions, most of which are invariant. For many analyses, such as population structure, diversity metrics, or genome-wide association studies, we are only interested in variable positions, specifically single nucleotide polymorphisms (SNPs).
To extract only SNPs, we use the `-v` snps option with **bcftools view**, which retains polymorphic sites of the specified type. Unlike the previous command where we used `-V` (uppercase) to exclude variant types, the lowercase `-v` explicitly includes only the types we want, in this case, SNPs.

```bash
#!/bin/bash

cd ~/project

# Load modules
module load BCFtools

# Make VCF with only variable sites (SNPs)
VCF_RAW="./vcf/project.raw.vcf.gz"
VCF_SNPS="./vcf/project.snps.vcf.gz"
bcftools view -v snps -Oz -o "$VCF_SNPS" "$VCF_RAW"
tabix -p vcf "$VCF_SNPS"
```

### 3. VCF with high quality genotypes in all individuals

In the final step, we generate a highly accurate VCF that retains only biallelic SNPs with high-quality genotypes and no missing data across all individuals. This kind of stringent filtering is useful when performing analyses that are sensitive to genotyping error or missingness. The filtering involves two stages:
First, we use **bcftools view** to keep only SNPs that are biallelic, using the flags `-m2` and `-M2` to specify a minimum and maximum of two alleles per site. Next, we pipe the result into **bcftools filter**, applying several quality thresholds:
- A minimum genotype quality (`GQ`) of 30, which helps ensure confidence in the called genotypes.
- A minor allele frequency (`MAF`) filter to retain variants with frequencies between 0.01 and 0.99, removing both extremely rare and fixed variants.
- A requirement that there are no missing genotypes across samples (`F_MISSING=0`).

Since this filtering is applied only to SNPs, we can use the VCF file `project.snps.vcf.gz` generated in the previous step as input.

```bash
#!/bin/bash

cd ~/project

# Load modules
module load BCFtools

# Make VCF with stringent filtering
VCF_SNPS="./vcf/project.snps.vcf.gz"
VCF_STRINGENT="./vcf/project.stringent.vcf.gz"
bcftools view -m2 -M2 -v snps "$VCF_SNPS" | bcftools filter -i 'MIN(GQ)>=30 && MAF>=0.01 && MAF<=0.99 && F_MISSING=0' -Oz -o "$VCF_STRINGENT"
tabix -p vcf "$VCF_STRINGENT"
``` 
