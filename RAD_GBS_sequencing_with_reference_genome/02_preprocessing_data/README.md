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

#### 2) Trimming  
#### 3) Remove Adapter sequences  
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

#### 6) Rename files
```bash
cp ./dedup/Sample1_CAGTGTGT-ATCACG.1.1.fq.gz ./names/[ID_samplename]_[location].1.fq.gz
cp ./dedup/Sample1_CAGTGTGT-ATCACG.2.2.fq.gz ./names/[ID_samplename]_[location].2.fq.gz
```
