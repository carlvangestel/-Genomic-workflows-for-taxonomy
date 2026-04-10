##TEST##



#### Inspect MAPQ distribution before deciding threshold
view MAPQ stats: samtools view input.bam | awk '{print $5}' | sort -n | uniq -c
prints a histogram of MAPQ values to assist in picking a sensible cutoff (e.g., if 90% are ≥30, that’s safe, DISCARD reads MAPQ=0 before calculating percentage).
run_MAPQthreshold.sh

#### Remove unmapped, secondary and supplementary alignments (-F 2308);  not properly paired reads -f 2; and retain only with a mapping quality (MAPQ) >30 (-q 30)
Sort and index bam file
run_filter_bam.sh

#### count how many reads (percentage) retained after filtering
run_readsRetained.sh


## Mapping Reads to the Reference Genome

To gain meaningful insights from sequencing data, we first need to determine where each read originates from in the genome. This process is called **mapping** or **alignment**, where sequencing reads are matched to a reference genome. There are many mapping tools available, often tailored for specific data types - for example, STAR or HISAT for RNA-seq, minimap2 for long reads, miniprot for protein sequences.   
In this workflow, we use **BWA-MEM**, a widely used algorithm within the BWA software suite, to align paired-end reads to a reference genome.


### 1. Index the Reference Genome

Before mapping, the reference genome must be indexed with BWA. Assuming `$GENOME` points to the reference FASTA file (e.g., `GENOME="./genome/refgenome.fasta"`), indexing is performed using:

```bash
bwa index $GENOME
```

For downstream applications, it is often useful to also generate a .fai index using SAMtools:

```bash
samtools faidx $GENOME
```

### 2. Map Paired-End Reads

Assuming that `$READ1` and `$READ2` point to the the forward and reverse fastq files (optionally gzipped) of an individual (e.g. `READ1="./reads/sample01_1.fq.gz"` and `READ2="./reads/sample01_2.fq.gz"`), `$SAMPLE` refers to the sampleID of an individual (`SAMPLE="sample01"`) and `$BAM_OUT` to the output file (`BAM_OUT="./bam/sample01.bam"`) the full command to map paired reads to a reference genome is:

```bash
bwa mem -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$GENOME" "$READ1" "$READ2" | samtools view -bS | samtools sort -o "$BAM_OUT"
```

Let's break down the command. The first part, `bwa mem "$GENOME" "$READ1" "$READ2"` does the actual alignment of the reads to the reference genome and outputs a SAM file—a text-based format containing alignment information. Since SAM files are large, they're typically converted to the more compact, binary BAM format using the program **samtools**. This is done on-the-fly using a pipe (`|`), which passes the SAM output directly into **samtools view**. The `-bS` flags tell **samtools** view to convert from SAM (`-S`) to BAM (`-b`). The `-R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA"`  option in **bwa mem** adds read group information to the header of the resulting BAM file. This includes:

- ID: read group ID (usually the same as the sampleID)
- SM: sampleID
- PL: sequencing platform (e.g., ILLUMINA)

Adding read group metadata is important for downstream tools especially when working with multiple samples. 

Finally, since most downstream tools require BAM files sorted by genomic position, we sort the output using **samtools sort** and save it as the final BAM file with `-o "$BAM_OUT"`. It's best to name the output file using the individual's ID, matching entries in the `samples.txt` file. 
Once the mapping is completed, an index file is created (**samtools index**) to enable fast and efficient access to specific regions within the BAM without reading the entire file using the command `samtools index "$BAM_OUT"`

While mapping can be done sequentially for each individual, it is more efficient on a computing cluster to run mappings in parallel. The mapping and sample processing commands are then specified in a script that is submitted to the cluster. The following script (mapping.sh) automates this by reading paired fastq.gz files from the file samples.txt, mapping them, and outputting sorted BAM files named after each sample ID.

```bash

#!/bin/bash

cd ~/project

# Load modules
module load BWA
module load SAMtools

# Read sample name from file based on array task ID
SAMPLE=$(sed -n "${PBS_ARRAYID}p" ./samples/samples.txt)

# Define genome file
GENOME="./genome/refgenome.fasta"

# Define input and output files
READ1="./reads/${SAMPLE}_1.fq.gz"
READ2="./reads/${SAMPLE}_2.fq.gz"
BAM_OUT="./bam/${SAMPLE}.bam"

# Run BWA MEM, convert to BAM and sort
bwa mem -t 8 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$GENOME" "$READ1" "$READ2" | samtools view -bS | samtools sort -o "$BAM_OUT"
samtools index "$BAM_OUT"
```

### 3. Remove PCR duplicates

Before calling single nucleotide polymorphisms (SNPs) a last preprocessing step is required being removal of PCR duplicates. During library preparation, PCR amplification can generate multiple reads from the same original DNA fragment. These are not independent observations and may bias variant calling, so only one copy should be retained. PCR duplicates can be identified when both reads in a pair have the same alignment start positions. They can be removed using **picard’s MarkDuplicates** tool (https://broadinstitute.github.io/picard/) using the following script:

```bash
java -Xmx8G -jar picard-tools-2.9.0/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT INPUT="$BAM_IN" OUTPUT="$BAM_RMD_OUT" METRICS_FILE="$BAM_RMD_METRICS"
```

This command removes duplicate reads from `$BAM_IN` (e.g., sampleID.bam) and writes the cleaned output to `$BAM_RMD_OUT` (e.g. `sampleID.rmd.bam`). It also produces a metrics file (`$BAM_RMD_METRICS`, e.g., `sampleID.rmd.metrics.txt`) reporting the proportion of duplicates found. Once completed, we again index the new duplicate-removed bam files as before `samtools index "$BAM_RMD_OUT"`. To save space, the original bam files can be removed. 

If you want to run this in parallel on a cluster, the following example script (picard_rmd.sh) can be used:

```bash
#!/bin/bash

cd ~/project

# Load modules
module load SAMtools

# Read sample name from file based on array task ID
SAMPLE=$(sed -n "${PBS_ARRAYID}p" ./samples/samples.txt)

# Define input and output files
BAM_IN="./bam/${SAMPLE}.bam"
BAM_RMD_OUT="./bam/${SAMPLE}.rmd.bam"
BAM_RMD_METRICS="./bam/${SAMPLE}.metrics.rmd.txt"

# Run Picard's MarkDuplicates tool
java -Xmx8G -jar picard-tools-2.9.0/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT INPUT="$BAM_IN" OUTPUT="$BAM_RMD_OUT" METRICS_FILE="$BAM_RMD_METRICS"

# Index the new duplicate-removed bam files
samtools index "$BAM_RMD_OUT"
```
