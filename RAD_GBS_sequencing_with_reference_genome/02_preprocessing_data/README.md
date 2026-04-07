Before aligning sequencing data to a reference genome and performing SNP identification, it is essential to assess the quality of the raw reads and carry out any necessary preprocessing. This ensures that only high-quality, reliable data are used in downstream analyses.

The first step in this process is to evaluate read quality using the software package FastQC. Detailed documentation and guidance can be found in the official FastQC manual: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

FastQC provides a range of quality metrics that should be carefully examined. In particular, verify that per-base quality scores remain above an acceptable threshold across the entire read length, noting that quality often declines toward the ends of reads. Additionally, inspect the presence of overrepresented sequences, assess the extent of missing or ambiguous data, and check for potential adapter contamination. These indicators help determine whether trimming or filtering steps are required before proceeding.
#### 1) Quality check reads 
module load FastQC
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

#### 3) trim T at end reverse read
explanation wet lab protocol


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
