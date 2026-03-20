# Allelic diversity

This repository contains the set of scripts used to **explore the allelic diversity in a large collection of isolates** based on a cg/wgMLST matrix and a metadata table.

##Scripts

### _allelic_heatmap.py_
The script _allelic_heatmap.py_ contains the set of functions necessary to build a heatmap with the allele distribution across a species tree. This script generates an **interactive heatmap** that can be easily explored by the user to detect potential recombination events between different groups of samples.

#### Input
- cg/wgMLST allelic matrix (*tsv)
- metadata table (*tsv)
- dendrogram (*nwk)

#### Output
- HTML with the tree and the heatmap
- *_allelic_detailed_matrix.tsv* - a table similar to the allele matrix where, besides the allele information, the user can find the information whether it is exclusive of a certain group (according to a user-selected metadata variable, e.g. lineage) or shared by different groups. 
- *_allele_summary.tsv* - a summary table indicating, for each group of samples, which alleles are shared between groups and which are exclusive of a certain group.
- *_allele_summary_counts.tsv* - a summary table with the counts of loci and alleles that are exclusive of each group or shared by some groups.
- *_allele_dominance.tsv* - a table indicating the absolute and relative counts of each allele in each group of samples. 

#### Usage
```
  -h, --help            show this help message and exit
  -a ALLELIC_MATRIX, --allelic_matrix ALLELIC_MATRIX
                        TSV allelic matrix
  -m METADATA, --metadata METADATA
                        Metadata TSV
  -c COLORS, --colors COLORS
                        (Optional) Colors TSV (column, category, color). If omitted, colors are auto-assigned and saved to <prefix>_colors.tsv
  -hcol HEADER, --header HEADER
                        Metadata column to use for categories
  -t TREE, --tree TREE  Newick tree file
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Output prefix
  -mode {frequency,count}, --mode {frequency,count}
                        Dominant category mode: 'frequency' (default) or 'count'
  -ncat NCAT, --ncat NCAT
                        Optional: if an allele is present in at least N categories, color it gray (nearly conserved)
  -ncat_dominant NCAT_DOMINANT, --ncat_dominant NCAT_DOMINANT
                        Optional: if an allele is dominant in at least N categories, color it gray (nearly conserved)
```

### _allelic_distribution.py_
The script _allelic_distribution.py_ contains the set of functions necessary to obtain information about the proportion of isolates from a given group (e.g. lineage) that harbor a certain allele.

#### Input
- cg/wgMLST allelic matrix (*tsv)
- metadata table (*tsv)

#### Output
- *_summary.tsv* - a table where each row corresponds to an allele and each column to a group of samples (according to a user-selected metadata variable, e.g. lineage), indicating what is the proportion of isolates of the group that contains the allele.
- *_distribution.tsv* - a table where each row corresponds to a locus and each column to a group, indicating the allele distribution per group.

#### Usage
```
  -h, --help            show this help message and exit
  -m METADATA, --metadata METADATA
                        [MANDATORY] Metadata file in .tsv format.
  -a ALLELES, --allele ALLELES
                        [MANDATORY] Allele matrix in .tsv format.
  -c GROUP_COLUMN, --group-column GROUP_COLUMN
                        [MANDATORY] Name of the metadata column with the groups of interest.
  -g GROUP_INTEREST, --group-interest GROUP_INTEREST
                        [OPTIONAL] Comma-separated list of groups of interest to select the alleles under analysis. If nothing is indicated, all goups
                        will be considered as of interest and all alleles will be analysed.
  -o OUTPUT, --output OUTPUT
                        [OPTIONAL] Tag for output files. default: Allele_distribution
```

## Command line examples


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

[Zohra Lodhia, Verónica Mixão, Joana Isidro, Rita Ferreira, Dora Cordeiro, Cristina Correia, Inês João, João Paulo Gomes, Maria José Borrego and Vítor Borges (2025) Advancing _Chlamydia trachomatis_ genomic surveillance and research with a novel core-genome MLST (cgMLST) approach. Research Square. doi: 10.21203/rs.3.rs-7743240/v1](https://www.researchsquare.com/article/rs-7743240/v1)

## Funding
This work was funded by national funds through FCT - Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001 (doi: 10.54499/2022.00851.CEECIND/CP1748/CT0001) and of the doctoral fellowship SFRH/BD/147446/2019 (doi: 10.54499/SFRH/BD/147446/2019), and by the European Union project “Sustainable use and integration of enhanced infrastructure into routine genome-based surveillance and outbreak investigation activities in Portugal” - GENEO [101113460] on behalf of the EU4H programme [EU4H-2022-DGA-MS-IBA-01-02].
