#### 1) FASTQC  

```bash
module load FastQC
sample="sample_1_R1.fastq.gz sample_1_R2.gz"
for i in $sample;
do fastqc ./$i -o ./FASTQC/FASTQC_raw
done
```

```bash
module load MultiQC
multiqc ./*_fastqc.zip
```
#### 2) Demultiplex 
#### 3) Trimming  
#### 4) Remove Adapter sequences  
#### 5) Remove PCR duplicates
#### 6) Rename files
