# Allelic diversity

This repository contains the set of scripts used to **explore the allelic diversity in a large collection of isolates** based on a cg/wgMLST matrix and a metadata table.

## Input/Output of _allelic_heatmap.py_
The script _allelic_heatmap.py_ contains the set of functions necessary to build a heatmap with the allele distribution across a species tree. This script generates an **interactive heatmap** that can be easily explored by the user to detect potential recombination events between different groups of samples.

#### Input
- allelic matrix (cg/wgMLST)
- metadata table
- dendrogram

#### Output
- HTML with the tree and the heatmap
- *_allelic_detailed_matrix.tsv* - a table similar to the allele matrix where, besides the allele information, the user can find the information whether it is exclusive of a certain group (according to a user-selected metadata variable, e.g. lineage) or shared by different groups. 
- *_allele_summary.tsv* - a summary table indicating, for each group of samples, which alleles are shared between groups and which are exclusive of a certain group.
- *_allele_summary_counts.tsv* - a summary table with the counts of loci and alleles that are exclusive of each group or shared by some groups.
- *_allele_dominance.tsv* - a table indicating the absolute and relative counts of each allele in each group of samples. 

## Input/Output of _allelic_distribution.py_
The script _allelic_distribution.py_ contains the set of functions necessary to obtain information about the proportion of isolates from a given group (e.g. lineage) that harbor a certain allele.

#### Input
- allelic matrix (cg/wgMLST)
- metadata table

#### Output
- *_summary.tsv* - a table where each row corresponds to an allele and each column to a group of samples (according to a user-selected metadata variable, e.g. lineage), indicating what is the proportion of isolates of the group that contains the allele.
- *_distribution.tsv* - a table where each row corresponds to a locus and each column to a group, indicating the allele distribution per group.

## Installation

These scripts are written in python 3.10. To facilitate their installation, we provide the conda installation file _allelic_diversity.yml_.

```
cd /PATH/TO/allelic_diversity/
conda env create -f allelic_diversity.yml
conda activate allelic_diversity
python allelic_heatmap.py -h
python allelic_distribution.py -h
```

## Citation
If you use any of these scripts, please do not forget to cite:   

Zohra Lodhia, Verónica Mixão, Joana Isidro, Rita Ferreira, Dora Cordeiro, Cristina Correia, Inês João, João Paulo Gomes, Maria José Borrego and Vítor Borges (2025) Exploring _Chlamydia trachomatis_ global genomic diversity with a novel core-genome MLST (cgMLST) approach. bioRxiv.

## Funding
This work was funded by national funds through FCT - Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001 (doi: 10.54499/2022.00851.CEECIND/CP1748/CT0001), and by the European Union project “Sustainable use and integration of enhanced infrastructure into routine genome-based surveillance and outbreak investigation activities in Portugal” - GENEO [101113460] on behalf of the EU4H programme [EU4H-2022-DGA-MS-IBA-01-02].
