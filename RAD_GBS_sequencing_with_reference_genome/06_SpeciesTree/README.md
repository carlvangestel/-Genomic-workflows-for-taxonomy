## A SNP-based species tree

Most species-tree inference methods use either gene trees or genetic sequences as input. These approaches generally assume that recombination is absent within these regions, such that each region represents a single underlying phylogenetic history. Because RADtags are typically shorter than 1200 bp, they often contain limited phylogenetic information, which can result in poorly resolved or unreliable gene-tree estimates. Conversely, concatenating all SNPs into a single alignment may lead to inflated node support and fails to account for potential variation in phylogenetic histories across different regions of the genome. An alternative approach that avoids these issues is Singular Value Decomposition Scores for Species Quartets (SVDQuartets) (Chifman and Kubatko, 2014), a species-tree inference method allowing each SNP to have its own evolutionary history, and estimates species-tree topology under the multispecies coalescent model. The latter however means no branch lengths can be inferred.
We will use PAUP* to estimate the SNP-based species tree using the SVDQuartets tool.

PAUP* does not accept vcf files and therefore we convert it to a NEXUS format using the script vcf2phylip.

```bash
python vcf2phylip.py --input ./vcf/project.HQ.minDP15.80shared.vcf --phylip-disable --nexus
```
This code will generate a NEXUS matrix named project.HQ.minDP15.80shared.nex. --phylipl-disable prevents the creation of the PHYLIP matrix (default) and the --nexus flag points out you wish to vcf to be transformed into a nexus file.

Next, inform PAUP* which samples belong to which taxa by adding the following text to the end of project.HQ.minDP15.80shared.nex file:
```
BEGIN SETS;
        TAXPARTITION SPECIES =
                ARC:1-5,
                BDA:6-10,
                BON:11-17,
                CC:18-22,
                CDC:23-27 78,
                CM:29-33,
                COR:34-38,
                GAR:39-43,
                GSL:44-48,
                MAR:49-50,
                PSA:51-55,
                RDA:56-59,
                SL:60-62,
                TRM:63-67,
                VAG:68-72,
                VDC:28 73-77;

  END;
```
Here, the first 5 sequences (1-5) of the nexus file belong to the taxa 'ARC', the next 5 sequences (6-10) to taxa 'BDA', and so on. 


1. Open the Nexus file NC_031969.f5.sub4.nex in PAUP*, and make sure that the option "Execute" is set in the opening dialog
2. If you have an outgroup: specify the two samples of Astatotilapia burtoni ("IZA1" and "IZC5") as the outgroup, click on "Define Outgroup..." in PAUP*'s "Data" menu, as shown below.
3. Click "SVDQuartets..." again in the "Analysis" menu. In the settings for the SVDQuartets analysis, again use "Evaluate all possible quartets" and make sure that "Distribute" is selected for "Handling of ambiguities". Unlike in the first analysis, now assign the samples to species by setting a tick in the last checkbox next to "Assign tips to species using taxon partition". The "SPECIES" taxon partition should already be selected in the drop-down menu to the right of it; this is the taxon set definition that we had added in Nexus format to the end of the alignment file. To also perform a bootstrapping analysis this time to assess node support for the species tree, set a tick in the checkbox next to "Perform bootstrapping". With 100 bootstrap replicates, the bootstrapping analysis might take around 20 minutes; if you would prefer not to wait that long you could set the number of bootstrap replicates to 50 instead. Another way to speed up the analysis is to use all CPUs available on your machine. To do so, simply click the "#CPUs" button in the bottom right of the settings window. This window should then look as shown in the next screenshot; click "OK" if it does.
4. You'll see that PAUP* has already generated a species tree, and that it reports the number of quartets that are comparable or incomparable with this tree. If all quartets would be comparable with this tree, this would indicate complete absence of incomplete lineage sorting and other processes that might violate the multi-species coalescent model, such as introgression or paralogous sequences in the alignment. However, these processes can not be excluded, and introgression is in fact very likely among the species included in this dataset.
5. Once bootstrapping has completed, PAUP* should print the same species tree once again, this time with bootstrap node support values as shown in the screenshot below.

Note: you can run the same analysis with all samples assigned as separate taxa to explore to what extent all samples of the same putative taxa form a monophyletic group.   

Reference.
Chifman, J. and L. Kubatko. 2014. Quartet inference from SNP data under the coalescent, Bioinformatics, 30(23): 3317-3324. 


