## SNP filtering

The VCF file generated in the previous step contains genotype information for all genomic positions, including both variable (polymorphic) and invariant sites, as well as genotype calls made with low quality. While this raw VCF (hence the proposed filename extension `.raw.vcf.gz`) is useful as a starting point, this unfiltered VCF is generally not suitable for most downstream analyses, which typically require only high-confidence, polymorphic sites.
Filtering a VCF file is a critical step that can greatly influence analytical results, especially in population genomics, phylogenetics, or variant association studies. Care must be taken when applying filters to avoid excluding important variants or retaining low-quality data. In the following sections, we will use **BCFtools and VCFtools** to generate a VCF containing only high quality SNPs. Next, we will set a lower limit to the number of individuals that need to share a SNP to avoid any bias due to missing data. Finally, as many genomic analyses may require a set of independent SNPs as input, we will further thin the vcf to reduce the amount of linkage disequilibrium (i.e. we will select (approximately) a single SNP per RADtag).

### 1. A VCF with high quality SNPs
In this first filter step we remove multiallelic sites, indels and multi-nucleotide polymorphisms and retain only biallelic SNPs with sufficient depth (to exclude possible sequencing errors are being interpreted as SNPs).

```bash
#!/bin/bash
module load BCFtools
VCF_RAW="./vcf/project.raw.vcf.gz"
VCF_HQ="./vcf/project.HQ.minDP15.vcf.gz"
bcftools view -V indels,mnps -v snps -m2 -M2 "$VCF_RAW" -Ou | bcftools +setGT -Ou -- -t q -n . -i 'FMT/DP<15 || FMT/GQ<30' | bcftools filter -i 'MAF>=0.01 && MAF<=0.99 -oz -o "$VCF_HQ" 
tabix -p vcf "$VCF_HQ"
```

The raw VCF file was filtered with bcftools to retain only high-quality biallelic SNPs and to mask unreliable genotype calls. First, bcftools view was used with the options -V indels,mnps and -v snps to exclude insertions/deletions (indels) and multi-nucleotide polymorphisms (MNPs) and to retain only single-nucleotide polymorphisms (SNPs). The options -m2 -M2 restricted the dataset to biallelic variants, meaning that only sites with exactly one reference allele and one alternative allele were kept. The filtered output was written in uncompressed BCF format (-Ou) and piped directly to the next processing step.

Next, the bcftools +setGT plugin was used to modify low-confidence genotype calls. The option -t q specified that genotypes should be filtered according to a quality expression, while -n . replaced genotypes failing the filter with missing values (./.). The expression 'FMT/DP<15 || FMT/GQ<30' identified genotype calls with sequencing depth below 15 reads or a genotype quality below 30, and these low quality genotypes were masked as missing. The resulting VCF was compressed in BGZF format using -Oz.

Finally, bcftools filter was used to constrain the minor allele frequency (MAF) to values between 0.01 and 0.99, thereby removing extremely rare variants and nearly fixed sites. The filtered variants were written to the output file specified by $VCF_HQ.

The final VCF file was indexed with tabix -p vcf, enabling efficient random access to genomic regions within the compressed VCF.

### 2. A VCF shared in 80% of individuals

RAD-seq datasets typically contain missing genotypes because restriction sites may be absent in some individuals, sequencing depth can vary among samples, and loci may not be recovered consistently across all samples. As a result, many SNPs are genotyped in only a subset of the samples. Therefore, to reduce the impact of missing data on downstream analyses, SNPs are filtered to retain only those present in, for example, at least 80% of the samples. This is equivalent to requiring a maximum missingness of 20% per locus. 
It is important to note that the 80% threshold is not a universal standard but rather an arbitrary chosen compromise between data completeness and data retention. The optimal threshold may differ between projects (depending on sequencing depth, number of samples, the study species, ...). It is often advisable to evaluate the effects of different missing data thresholds (e.g., SNPs shared in 70%, 80%, 90%, or even 100% of all samples) on the number of retained SNPs. The choice of threshold therefore reflects a balance between maximizing the number of SNPs available for analysis and minimizing potential biases introduced by excessive missing data.

```
#!/bin/bash
module load BCFtools
VCF_HQ="./vcf/project.HQ.minDP15.vcf.gz"
VCF_80shared="./vcf/project.HQ.minDP15.80shared.vcf.gz"
bcftools view -i 'F_MISSING<=0.20' "$VCF_HQ" -Oz -o "$VCF_80shared"
tabix -p vcf "$VCF_80shared"
```
Here, we used previous vcf containing high quality SNPs as input and added the filter option 'F_MISSING' to allow you to set the criteria of how much missing data you allow per site (maximim 20% missing data per site was allowed in the code above).

### 3. A thinned VCF

RAD loci often contain multiple SNPs that are physically close and therefore highly linked. Many downstream analyses (e.g., PCA, STRUCTURE, ADMIXTURE, DAPC) assume SNPs are approximately independent. Thinning helps reduce the overrepresentation of genomic regions containing many SNPs and decreases the effect of linkage disequilibrium (LD). We will use the option 'thinning' implemented in VCFTools to retain only SNPs that are at least 1200 bp apart from one another (approximate length of a single RAD tag). Note that we use a distance-based thinning as an approximation and do not estimate the actual LD between SNPs which is often sufficient for most analyses.       

```
#!/bin/bash
module load VCFtools
VCF_80shared="./vcf/project.HQ.minDP15.80shared.vcf.gz"
"$VCF_thin"="./vcf/project.HQ.minDP15.80shared.thin.vcf.gz"
vcftools --vcf "$VCF_80shared" --thinning 1200 --recode
mv output.recode.vcf "$VCF_thin"
tabix -p vcf "$VCF_thin"
```




