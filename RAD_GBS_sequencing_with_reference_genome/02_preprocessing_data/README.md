#### 1) FASTQC  
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
#### 2) Demultiplex 
```bash
process_radtags -P -1 ./library1_R1.fastq.gz -2 ./library1_R2.fastq.gz -o ../samples/  -b ./barcode_lib1.txt -e sbfI -r -c -q --inline_index
-i gzfastq
```
#### 3) Trimming  

#### 4) Remove Adapter sequences  

#### 5) Remove PCR duplicates

#### 6) Rename files
