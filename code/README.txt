Caudata phylogeny

version 0.0 (This was a test run only!)

- PyPHLAWD -

## The configuration file used as is, with no changes to the parameters.
+ spiritusanti$ python3 src/setup_clade.py "Caudata" vrt.05102018.db SALAS/ src/conf.py
+ spiritusanti$ python3 src/find_good_clusters_for_concat.py --dir SALAS/Caudata_8293/ --database vrt.05102018.db
	+ Do you want to rename these clusters? y
	+ Do you want to make trees and trim tips for these gene regions y/n n
	+ Do you want to concat? y/n y
	+ Do you want to make a constraint? y/n n

+ spiritusanti$ python3 src/add_outgroup_to_matrix.py -b vrt.05102018.db -m SALAS/Caudata_8293/Caudata_8293_outaln -p SALAS/Caudata_8293/Caudata_8293_outpart -t 8353 -mn 2 -mx 5 -o Caudata_outgrp
## This adds five species as outgroup, overlapping with at least two genes. 
## The outgroup here is 8353 (NCBI Taxonomy ID), which correspondes to Xenopus.

- RAxML -

+ spiritusanti$ python3 src/change_ncbi_to_name_tre_fromurl.py http://141.211.236.35:10999 ../Documents/1.PROYECTOS/9.SALAMANDRAS/RAxML_bestTree.Caudata_v0.tre > Caudata_v0_names.tre

version 1.0

- PyPHLAWD -

## This builds an updated database of vertebrate (VRT) sequences from NCBI (created 28-06-2020).
+ spiritusanti$ ./phlawd_db_maker vrt ~/PyPHLAWD/vrt.06302020.db
 
## The configuration file used as is, with no changes to the parameters.
## Exclude patterns: 'sp.', 'var.', 'cf.', 'aff.'
+ spiritusanti$ python3 src/setup_clade.py "Caudata" vrt.06302020.db New_SALAS/ src/conf.py

## This finds good clusters (orthologs) and creates a concatenated alignement.
+ smallest_cluster = 20; cluster_prop = 0.2
+ Found 75 good clusters!
## NOTE: not sure how to get the names of the genes for each cluster, but from a quick look there are both mtDNA and nDNA.

+ spiritusanti$ python3 src/find_good_clusters_for_concat.py --dir New_SALAS/Caudata_8293/ --database vrt.06302020.db
 	+	Do you want to rename these clusters? y
	+	Do you want to make trees and trim tips for these gene regions y/n n
	+	Do you want to concat? y/n y
	+	Do you want to make a constraint? y/n n

## I added an outgroup to the concatenated alignment.
## Added five species as outgroup, overlapping with at least two genes. 
## The outgroup here is 8353 (NCBI Taxonomy ID), which correspondes to Xenopus.
+ spiritusanti$ python3 src/add_outgroup_to_matrix.py -b vrt.05102018.db -m New_SALAS/Caudata_8293/Caudata_8293_outaln -p New_SALAS/Caudata_8293/Caudata_8293_outpart -t 8353 -mn 2 -mx 5 -o Caudata_outgrp

- RAxML -
## The ML inference was performed using v.8 as implemented in CIPRES.
raxmlHPC-HYBRID -T 4 -f a -N autoMRE -n Caudata_v1 -s infile.txt -c 25 -m GTRCAT -p 12345 -q part.txt -k -x 12345 -o 191247,224340,8364,8355,8354

## Relevant options:
+ -f a : rapid bootstrap analysis / best-scoring ML tree
+ -o : five species as outgroup (NCBI Taxonomy ID)
+ -s : infile.txt (== concatenated alignment of 93 clusters)
+ -m : GTRCAT model
+ -q : partition file for the concatenated alignment
+ -N autoMRE : bootstrapping is halted automatically by RAxML

## The phylogeny has names corresponding to the NCBI Taxonomy IDs (it is easier to change these until the very end).
## The NCBI Taxonomy IDs were transformed into taxa names in R (not with python as implemented in PyPHLAWD). 
## For this I used the 'taxize' package and assigned family and species names to each taxa.
## NOTE: I used my personal ENTREZ API KEY to increase the rate of requests/second to process all the IDs (otherwise it crashes).
+ usethis::edit_r_environ(); ENTREZ_KEY=myKEY

## The monophyly check basically consists of the following:
## 1. For a given taxon, get all the species belonging to that taxon and then find their MRCA in the tree.
## 2. Extract a sub-tree from the taxon's MRCA and count the number of observed tips in the sub-tree.
## 3. If the taxon is monophyletic, the number of observed tips in the sub-tree is exactly the same as the number of species for the taxon.
+     If the taxon is not-monophyletic, the number of observed tips in the sub-tree is greater than the number of species for the taxon.

## Results of monophyly check for genera, along with decision (discussed on 20 August, 2020)
+ Cynops		OK!
+ Pachytriton		Constrain
+ Paramesotriton 	Constrain
+ Pseudotriton		Constrain
+ Eurycea		Constrain
+ Pseudoeurycea		Constrain
+ Bolitoglossa		Constrain
+ Onychodactylus	Constrain
+ Hynobius		Constrain


version 2.0-5.0

- R -

## The alignment now has the taxa names (instead of NCBI's taxonomic IDs). 
## This way it is easier to set the monophyletic constraints on the tree inference.

## Get extra species and align to profile.
+ Aquiloeurycea cafetalera: 16S - MT524627.1.
+ Chiropterotriton nubilus: 16S - MK335405.1; COI - MK335251.1

## Also merge some clusters for two markers:
+ 16S - two clusters
+ COI - two clusters

- AliView -

## Trim the alignment in AliView to preserve partitions -> generate new partitions file.
## Align second_* into profile.
## Align sequences for extra species into profile.

## Set constraints tree for RAxML
## Build the tree manually with the monophyletic constraints mentioned above.

- RAxML -

## The ML inference was performed using v.8 as implemented in CIPRES.
+ raxmlHPC-HYBRID -T 4 -f a -N autoMRE -n Caudata_v2 -s infile.txt -c 25 -m GTRCAT -p 12345 -q part.txt -k -x 12345 -o 191247,224340,8364,8355,8354


version PIS

