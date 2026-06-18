# Restriction-Associated DNA sequencing and Genotyping-by-Sequencing

## Intro
This section describes the workflow to infer patterns of species divergence using **whole-genome sequencing data of multiple individuals**. The analysis assumes the availability of Illumina or DNB sequencing reads covering the entire genome of multiple individuals of the same or related species in fastq format, as well as a reference genome of the focal or a closely related species. The workflow describes how to preprocess (clean) the raw sequencing reads, map the reads to the reference genome, identify single nucleotide polymorphisms (SNPs), perform different types analysis on the genetic relationship between the individuals and visualize the results. The workflow is provided in different components that are generally run sequentially, but can be performed individually given the correct input file formats.Description of this workflow step, tools, commands, inputs, and outputs.

This section describes the workflow to infer patterns of species divergence using **reduced-representation sequencing data**, such as **Restriction site-Associated DNA sequencing (RADseq)** or **Genotyping-by-Sequencing (GBS)**. The analysis assumes the availability of Illumina short sequencing reads generated from restriction enzyme–based libraries in FASTQ format for multiple individuals of the same or related species. Although de novo approaches are also possible, this pipeline uses a reference genome of the focal species or a closely related species. The workflow describes how to preprocess (clean) the raw sequencing reads, demultiplex samples, map reads to a reference genome, and call single nucleotide polymorphisms (SNPs) from the sequenced restriction-associated regions. These SNPs are then used to perform analyses of genetic variation and relationships among individuals, such as population structure, genetic distance estimation, and phylogenetic inference, followed by visualization of the results.

### Workflow components
- [1. Directory structure, data format and tools](./01_dataformat_and_tools/)
- [2. Preprocess reads](./02_preprocessing_data/)
- [3. Map reads](./03_read_mapping)
- [4. Calling Single Nucleotide Variants (SNPs)](./04_SNP_calling)
- [5. Filter SNPs](./05_SNP_filtering)
- [6. SNP-based species tree using SVDQuartets](./06_SpeciesTree)
