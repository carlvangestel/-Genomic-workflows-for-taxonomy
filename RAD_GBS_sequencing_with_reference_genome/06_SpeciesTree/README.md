## A SNP-based species tree

Most species-tree inference methods use either gene trees or genetic sequences as input. These approaches generally assume that recombination is absent within these regions, such that each region represents a single underlying phylogenetic history. Because RADtags are typically shorter than 1200 bp, they often contain limited phylogenetic information, which can result in poorly resolved or unreliable gene-tree estimates. Conversely, concatenating all SNPs into a single alignment may lead to inflated node support and fails to account for potential variation in phylogenetic histories across different regions of the genome. An alternative approach that avoids these issues is Singular Value Decomposition Scores for Species Quartets (SVDQuartets) (Chifman and Kubatko, 2014), a species-tree inference method allowing each SNP to have its own evolutionary history, and estimates species-tree topology under the multispecies coalescent model.


https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/
