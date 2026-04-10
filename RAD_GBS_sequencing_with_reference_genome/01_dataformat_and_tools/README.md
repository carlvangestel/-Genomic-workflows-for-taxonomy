## Set-up the directory structure

Before we start, it is most convenient to set up an organized directory structure to store the different datatypes used throughout this workflow. Using the bash command `mkdir <foldername>`, you can create the following directories in your main working directory, which we refer to as `~/project`:


| Directory      | Description |
|----------------|-------------|
| `bam`          | Stores the alignment files in BAM format. |
| `clean`        | Stores cleaned reads (after adapter trimming and quality filtering), if preprocessing is performed. |
| `genome`       | Stores the reference genome sequence in FASTA format. |
| `reads`        | Stores the raw sequencing reads in `fastq.gz` format. |
| `samples`      | Contains a file `samples.txt` listing all sample names (individual IDs). These names should match the basenames of the files in both the `reads` (or `clean`) directories. |
| `scripts`      | Stores your custom scripts used in the workflow. |
| `vcf`          | Stores Variant Call Format (VCF) files containing SNP and other variant information for all individuals. |

## Required programs and tools

| Tool      | Description |
|-----------|-------------|
| BCFTOOLS  | SNP calling and filtering vcf files (https://www.htslib.org/) |
| BWA       | Mapping low-divergent short reads to a reference genome (https://github.com/lh3/bwa) |
| SAMTOOLS  | Suite of programs for interacting with high-throughput sequencing data (https://www.htslib.org) |
|FastQC| Versatile tool to perform quality control checks on raw sequence data (https://www.bioinformatics.babraham.ac.uk/projects/fastqc)|
|MultiQC| Summarizes multiple FastQC outputs into a single interatcive report (https://github.com/MultiQC/MultiQC)|
|Stacks| Tools specifically designed for processing restriction enzyme-based data (https://catchenlab.life.illinois.edu/stacks)|
|Fastp| Processing tool to filter, trim or merge sequence reads (https://github.com/opengene/fastp)|
|PAUP*| Phylogenetic software program to estimate a species trees from SNPs using the SVDQuartets method (https://paup.phylosolutions.com)| 
| R packages| *adegenet* (Exploratory Analysis of Genetic and Genomic Data): contains a set of explorative tools such as multivariate methods to provide low-dimensional visualization of SNP variaton (https://cran.r-project.org/web/packages/adegenet/index.html)|

## Input data formats

### Sequencing reads

The standard data format produced by Illumina or DNB sequencing platforms consists of a **pair of FASTQ files** for each sample, typically distinguished by a “1” and “2” in their filenames. These two files contain the **forward** and **reverse** sequencing reads of the same DNA fragments.  

Importantly, the order of reads is preserved between the two files: the *n*-th record in the R1 file corresponds to the *n*-th record in the R2 file, representing the two ends of the same DNA fragment.

Each sequencing read in FASTQ format is represented by **four lines**:

1. **Read identifier** (starts with `@`)
2. **Nucleotide sequence**
3. **A plus sign (`+`)**
4. **Per-base quality scores** (ASCII-encoded Phred scores)

Phred scores reflect the probability that a base call is incorrect (see [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score)).

Here are some example FASTQ entries

**Forward read** (`sampleID_1.fq.gz`)

```text
@E250063195L1C001R00300000730/1
ATGTCAGATAAATTACTGGTTCCTAAGTTACGAATTAGCTAACCTACTTTTTTCACGTGTTAAAATACAACAATAACATTCATGTACTGCCATTTGCGTCGACCGGCAACGCTAATGTCC
+
FE@EFFFE8FFBE;EFFFFD:FFEF:FEEFFFFEEBFBCFFCEFCFCFCF4EFFFFFFFEFFF>FAFFFFFFBDFFFFFCEFEFDFFFFFFFF:BFFFEFEEFFF=F@FFFF>EFCFCFF
```

**Reverse read** (`sampleID_2.fq.gz`):

```text
@E250063195L1C001R00300000730/2
ATACCAGAAGACGGTCAGTGTCGATATTAAAACCTCTTCTGCCCATCCTATCACGTGATTTGACAGTATAACGGACATTAGCGTTGCCGGTCGACGCAAATGGCAGTACATGAATGTTAT
+
FFFFFFEFCEEGEFEDFF@FDFFFFF?FEFEFFFFFF<BFFFDFFFFF8GFFFFFFDF6;EFFFFF=F;GFEFAFFF9FFFFFFBGFGEFFFEFFFFFF2AGFFFF9GFF%GFF6E(;F?
```

### Genome (fasta) 

The genome sequence, stored in the `./genome` directory, generally comes in regular fasta format with each genomic sequence (contig, scaffold or chromosome) coming in two lines. The first line, starting with “>”, gives the name of the sequence, while the second line gives the actual nucleotide sequence.  

### Samples

The `./samples` folder should contain a simple text file, called `samples.txt`, that lists all sample names. Important, the sample names should match the basename of the raw (or cleaned) fastq files. Hence, if the fastq files of the first sample are `sample01_1.fq.gz` and `sample01_2.fq.gz`, the sample name in the `samples.txt`file should be **sample01** . 
