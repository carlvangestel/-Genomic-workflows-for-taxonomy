## A SNP-based species tree

Most species-tree inference methods use either gene trees or genetic sequences as input. These approaches generally assume that recombination is absent within these regions, such that each region represents a single underlying phylogenetic history. Because RADtags are typically shorter than 1200 bp, they often contain limited phylogenetic information, which can result in poorly resolved or unreliable gene-tree estimates. Conversely, concatenating all SNPs into a single alignment may lead to inflated node support and fails to account for potential variation in phylogenetic histories across different regions of the genome. An alternative approach that avoids these issues is Singular Value Decomposition Scores for Species Quartets (SVDQuartets) (Chifman and Kubatko, 2014), a species-tree inference method allowing each SNP to have its own evolutionary history, and estimates species-tree topology under the multispecies coalescent model. The latter however means no branch lengths can be inferred.
We will use PAUP* to estimate the SNP-based species tree using the SVDQuartets tool.

PAUP* does not accept vcf files and therefore we convert it to a NEXUS format using the script vcf2phylip.

python vcf2phylip.py --input myfile.vcf --phylip-disable --nexus

Example 4: If, for example, you wish to disable the creation of the PHYLIP matrix and only create a NEXUS matrix:

python vcf2phylip.py --input ./vcf/project.HQ.minDP15.80shared.vcf --phylip-disable --nexus
# Which is equivalent to:
python vcf2phylip.py -i myfile.vcf -p -n
# This command will create only a NEXUS matrix called myfile_min4.nexus




Reference.
Chifman, J. and L. Kubatko. 2014. Quartet inference from SNP data under the coalescent, Bioinformatics, 30(23): 3317-3324. 


