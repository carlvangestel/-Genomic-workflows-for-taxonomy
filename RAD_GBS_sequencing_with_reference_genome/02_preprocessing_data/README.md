## Pre-processing of raw sequence reads

Before aligning sequencing data to a reference genome and performing SNP identification, it is essential to assess the quality of the raw reads and carry out any necessary preprocessing. This ensures that only high-quality, reliable data are used in downstream analyses.

#### Step 1. Demultiplex libraries
During library construction each sample was provided a unique barcode and pooled into a common library that is being sequenced. Therefore, as a first step we will demultiplex our sequenced library to separate the pooled NGS data back into individual sample FASTQ files based on these unique barcodes. We use the function process_radtags of the software package STACKS, which is specifically built to handle data derived from restriction-enzyme based protocols, to simultaneously demultiplex and perfrom a first preliminary clean-up of raw sequencing reads.    

```bash
process_radtags -P -1 ./library1_R1.fastq.gz -2 ./library1_R2.fastq.gz -o ../samples/  -b ./barcode_lib1.txt -e sbfI -r -c -q --inline_index
-i gzfastq
```
The process_radtags command is configured to handle paired-end sequencing data, as indicated by the -P flag. The forward and reverse reads are provided via the -1 and -2 options, respectively, pointing to gzipped FASTQ files (library1_R1.fastq.gz and library1_R2.fastq.gz). Processed output files are written to the directory specified by -o (../samples/). Several options are used to improve data quality and retention. The -r flag enables the rescue of barcodes and RAD-tags, allowing reads with minor sequencing errors in the barcode or restriction site to still be retained. The -c option removes reads containing uncalled bases (Ns), ensuring cleaner data, while -q applies a sliding window quality filter to discard low-quality reads. The --inline_index option specifies that one barcode is located within the read sequences themselves and the other in a separate index reads<sup>(*)</sup>. Finally, the -i gzfastq flag defines the input format as gzipped FASTQ files.

_(*) Note that the flag '--inline_index is entirely determined on how you construct your library in the wetlab._

The barcode file ('barcode_lib1.txt') contains a list of all barcodes needed to assign reads to each sample. We constructed our libraries in such a way the first index is sample-specific and actually part of the DNA RADtag sequence itself, while the other barcode is library specific and sequenced via a separate index primer. For example, the list below refers to a specific library (barcode 'ACAGTG') and contains 5 samples, each with a unique barcode ('ACACTGAC', 'ACGTAGCA', etc.)

```bash
ACACTGAC	ACAGTG
ACGTAGCA	ACAGTG
CACACAGT	ACAGTG
CAGTCTCA	ACAGTG
GTACTCGT	ACAGTG
```
A typical read in a fastq file of such a library will look like: 

```bash
@LH00478:285:2273MCLT1:2:1103:10477:28895 1:N:0:ACAGTG 
ACACTGACTGCAGGTACATGGCAGACCATCGTAAGAGTTGTAAAACGTTTAAGGGAGACGGACTGTGTCAGCCGACCTCGAGCACGTAGACCTCGTAATGTAGGACGCAAAGTGCAACCGGAAGATGTGCTAGCATACGCTC 
```
ACAGTG (in header) refers to the index barcode and all reads of all samples belonging to this library will have this barcode in the header of the sequence. ACACTGAC (at start of the read) refers to the inline barcode and is unique for this sample (all reads belonging to this sample will contain this inline barcode).


Overall, this command takes paired-end gzipped FASTQ data, demultiplexes reads into individual samples based on inline barcodes, performs cleaning and quality filtering, attempts to rescue slightly imperfect reads, and outputs high-quality, sample-specific read files to the designated directory.

#### Step 2. 'T' trimming 
This step might not be necessary for all protocols as it depends on how your libraries were constructed. During library preparation, we ligated the P2 adaptor to our RADtag via 'A-tailing'. Hence, all our reverse reads start with a 'T', which is not part of the true DNA sequence and should be removed from our reads.

```bash
for f in *.2.fq.gz
do
    fastp -i "$f" -o "${f%fq.gz}.trimmed.fq.gz" --trim_front1 1
done
```

#### Step 3. Quality check of raw reads
The first step in this process is to evaluate read quality using the software package FastQC. Detailed documentation and guidance can be found in the official FastQC manual: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

FastQC provides a range of quality metrics that should be carefully examined. In particular, verify that per-base quality scores remain above an acceptable threshold across the entire read length, noting that quality often declines toward the ends of reads. Additionally, inspect the presence of overrepresented sequences, assess the extent of missing or ambiguous data, and check for potential adapter contamination. These indicators help determine whether trimming or filtering steps are required before proceeding. We will run 
 
```bash
lib="library1_R1.fastq.gz library1_R2.gz"
for i in $lib;
do fastqc ./$i -o ./FASTQC/FASTQC_raw
done
```

module load MultiQC
```bash
multiqc ./*_fastqc.zip
```

#### Step 3) trim T at end reverse read
explanation wet lab protocol
optional, depending on the library 

#### 3) Trimming and Remove Adapter sequences  
```bash
java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
```

#### 4) Demultiplex 
module load STACKS
```bash
process_radtags -P -1 ./library1_R1.fastq.gz -2 ./library1_R2.fastq.gz -o ../samples/  -b ./barcode_lib1.txt -e sbfI -r -c -q --inline_index
-i gzfastq
```

sample_ACGTAGCA-CAGATC.1.fq.gz      
sample_ACGTAGCA-CAGATC.2.fq.gz      
sample_TGTGTGAC-CAGATC.rem.1.fq.gz
sample_ACGTAGCA-CAGATC.rem.2.fq.gz

#### 5) Remove PCR duplicates
clone_filter -1 ./Sample1_ACACGACA-ACAGTG.1.fq.gz -2 ./Sample1_ACACGACA-ACAGTG.2.fq.gz -o ../dedup -i gzfastq -y gzfastq  
Explanation which methods allow this and which not, refer to Nature Reviews Genetics  

#### 6) Rename files
```bash
cp ./dedup/Sample1_CAGTGTGT-ATCACG.1.1.fq.gz ./names/[ID_samplename]_[location].1.fq.gz
cp ./dedup/Sample1_CAGTGTGT-ATCACG.2.2.fq.gz ./names/[ID_samplename]_[location].2.fq.gz
```
