## Pre-processing of raw sequence reads

Before aligning sequencing data to a reference genome and performing SNP identification, it is essential to assess the quality of the raw reads and carry out any necessary preprocessing. This ensures that only high-quality, reliable data are used in downstream analyses.

#### Step 1. Quality check of raw reads
Prior to processing the data, we evaluate read quality using the software package FastQC. Detailed documentation and guidance can be found in the official FastQC manual (https://www.bioinformatics.babraham.ac.uk/projects/fastqc).

FastQC provides a range of quality metrics that should be carefully examined. In particular, verify that per-base quality scores remain above an acceptable threshold across the entire read length, noting that quality often declines toward the ends of reads. Additionally, inspect the presence of overrepresented sequences, assess the extent of missing or ambiguous data, and check for potential adapter contamination. These indicators help determine whether trimming or filtering steps are required before proceeding.
 
```bash
lib="library1_R1.fq.gz -2 ./library1_R2.fq.gz"
for i in $lib;
do fastqc ./$i -o ./FASTQC/
done
```
Rather than exploring every single output separately, we will use MultiQC to compile this info into a single interactive htlm report.
```bash
multiqc ./FASTQC/*_fastqc.zip
```

#### Step 2. Demultiplex libraries
During library construction each sample was assigned a unique barcode and then pooled into a common library for sequencing, resulting in one fastq file provided by the sequencing facility. Therefore, we will first demultiplex our sequenced library to separate the pooled NGS data back into individual sample fastq files based on these unique barcodes. We use the function process_radtags of the software package Stacks, which is designed for restriction enzyme–based data and allows both demultiplexing and preliminary quality filtering of raw reads.    

```bash
process_radtags -P -1 ./library1_R1.fq.gz -2 ./library1_R2.fq.gz -o ../samples/  -b ./barcode_lib1.txt -e sbfI -r -c -q --inline_index
-i gzfastq
```
The process_radtags command is configured to handle paired-end sequencing data, as indicated by the -P flag. The forward and reverse reads are provided via the -1 and -2 options, respectively, pointing to gzipped fastq files (library1_R1.fastq.gz and library1_R2.fastq.gz). The -e option refers to the specific restriction enzyme used to digest the genome and the -b flag provides sample and library specific barcodes. Processed output files are written to the directory specified by -o (../samples/). Several options are used to improve data quality and retention. The -r flag enables the rescue of barcodes and RAD-tags, allowing reads with minor sequencing errors in the barcode or restriction site to still be retained. The -c option removes reads containing uncalled bases (Ns), ensuring cleaner data, while -q applies a sliding window quality filter to discard low-quality reads when the average quality score within a window drops below a certain threshold. The --inline_index option specifies that one barcode is located within the read sequences themselves and the other in a separate index reads<sup>(*)</sup>. Finally, the -i gzfastq flag defines the input format as gzipped fastq files.

_(*) Note that the flag '--inline_index' is entirely determined on how you construct your library in the wetlab. Read the manual to choose the appropriate format that suits your library._

The barcode file ('barcode_lib1.txt') contains a list of all barcodes needed to assign reads to each sample. We constructed our libraries in such a way the first index is sample-specific and actually part of the DNA RADtag sequence itself, while the other barcode is library specific and sequenced via a separate index primer on the Illumina sequencer. For example, the list below refers to a specific library (barcode 'ACAGTG') and contains 5 samples, each with a unique barcode ('ACACTGAC', 'ACGTAGCA', etc.)

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

#### Step 3. Remove PCR duplicates
Many restriction enzyme–based protocols include a PCR amplification step, which can introduce biases and errors in downstream analyses. A common issue associated with PCR duplicates is the miscalling of heterozygotes as homozygotes during SNP detection. Only wetlab protocols that incorporate random shearing of digested DNA (such as traditional RADseq) allow to identify and remove PCR duplicates. For a nice and detailed overview of different restriction enzyme-based digestion protocols we refer to Andrews et al. (2016) and Davey et al. (2011). To remove PCR duplicates we use the 'clone_filter' function implemented in the Stacks software package. 'Clone_filter' operates by comparing paired-end reads and grouping those that have similar sequences and share identical start and end positions. When multiple read pairs are found to be identical, they are considered PCR duplicates—artifacts of amplification rather than independent DNA molecules. The program retains only one representative read pair from each group and discards the others. This approach relies on the assumption that randomly sheared fragments are unlikely to produce identical start and end coordinates unless they are duplicates generated during PCR.

```bash
clone_filter -1 ./Sample1_ACACGACA-ACAGTG.1.fq.gz -2 ./Sample1_ACACGACA-ACAGTG.2.fq.gz -o ../dedup -i gzfastq -y gzfastq  
```
In this command, the -1 and -2 flags specify again the input files corresponding to the forward and reverse reads respectively. The -o flag defines the output directory where the filtered (de-duplicated) reads will be written. The -i option indicates the format of the input files, while the -y option specifies the format of the output files. 

#### Step 4. 'T' trimming 
This step might not be necessary for all protocols as it depends on how your libraries were constructed. During library preparation, we ligated the P2 adaptor to our RADtag via 'A-tailing'. Hence, all our reverse reads start with a 'T', which is not part of the true DNA sequence and should be removed from our reads.

```bash
for f in *.2.fq.gz
do
    fastp -i "$f" -o "${f%.fq.gz}.trimmed.fq.gz" --trim_front1 1
done
```

#### Step 5. Trimming and Remove Adapter sequences  
Depending on the output of the quality report you may need to trim reads for low quality bases or remove adapter sequences using Trimmomatic. 

```bash
java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
```
The command runs Trimmomatic in paired-end (PE) mode to perform adapter trimming on sequencing reads. The file names refer respectively to the forward and reverse input reads. The four output files define how reads are written after processing: output_forward_paired.fq.gz and output_reverse_paired.fq.gz contain read pairs where both mates survive filtering, whereas output_forward_unpaired.fq.gz and output_reverse_unpaired.fq.gz contain reads whose mate was discarded during trimming.

The ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 step specifies adapter removal. Here, TruSeq3-PE.fa is the FASTA file containing adapter sequences to be clipped. The value 2 sets the maximum number of mismatches allowed in the initial seed match between read and adapter. The value 30 is the palindrome clip threshold, used when detecting adapter contamination in paired-end reads by aligning read pairs against each other; higher values make clipping more stringent. The value 10 is the simple clip threshold, which applies to matches between the read and adapter sequences in a straightforward alignment; again, higher values increase stringency. Together, these parameters control how sensitively and strictly adapter sequences are detected and removed from the reads.

#### Step 6. Rename files
At this stage, you may want to rename your files into something more recognizable and conventient for you in downstream analyses, like only an ID name and location name for example.

```bash
cp ./dedup/Sample1_CAGTGTGT-ATCACG.1.1.fq.gz ./names/[ID_samplename]_[location].1.fq.gz
cp ./dedup/Sample1_CAGTGTGT-ATCACG.2.2.fq.gz ./names/[ID_samplename]_[location].2.fq.gz
```
<br><br><br>
_References.  
Andrews, K. R., Good, J. M., Miller, M. R., Luikart, G., & Hohenlohe, P. A. (2016). Harnessing the power of RADseq for ecological and evolutionary genomics. Nature Reviews Genetics, 17(2), 81–92. https://doi.org/10.1038/nrg.2015.28  
Davey, J. W., Hohenlohe, P. A., Etter, P. D., Boone, J. Q., Catchen, J. M., & Blaxter, M. L. (2011). Genome-wide genetic marker discovery and genotyping using next-generation sequencing. Nature Reviews Genetics, 12(7), 499–510. https://doi.org/10.1038/nrg3012_

