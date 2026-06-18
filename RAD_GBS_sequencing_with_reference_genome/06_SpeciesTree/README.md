## A SNP-based species tree

The VCF
Most methods for species-tree inference with the multi-species coalescent model use gene trees or gene sequences as input. In both cases it is assumed that within-gene recombination is absent and that therefore the entire gene has only a single true phylogeny. This assumption may often not be met, particularly when population sizes are large or speciation is frequent, or when the gene sequences used for inference are in fact concatenated sequences of exons that are distant to each other on the chromosome. In those cases, species-tree inference based on genes may potentially be affected by the same issues that affect concatenation methods more generally, leading to overconfidence in node support and unreliable tree inference. An approach that completely avoids these issues is species-tree inference based on single-nucleotide polymorphisms (SNPs), where single SNPs instead of genes are considered as phylogenetic markers, and thus each SNP is explicitely allowed have its own phylogenetic history. One tool that implements SNP-based phylogenetic inference is SVDQuartets (Chifman and Kubatko 2014), which estimates the topology of the species tree under the multi-species-coalescent model.

### 1. SVDquartets: Singular Value Decomposition Scores for Species Quartets

https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/
