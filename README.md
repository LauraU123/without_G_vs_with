# Comparing RSV Phylogeny with and without the G gene


This workflow constructs phylogenetic analyses for RSV genomes with G gene excluded and included.
Thus, the impact of the G gene on the tree structure can be visualised. 

### Input Data

The worklfow takes as input RSV sequences (fasta), metadata (tsv), as well as references (gbk and fasta).

### Output Data

The output are annotated RSV tree files constructed based on full genomes with and without the G gene respectively.
These two trees are constructed based on the same data, and thus tanglegrams comparing the trees can be constructed.


### Running the workflow

This workflow can be run from the command line: snakemake --cores all.

Builds of interest can be specified in the configfile.yaml in the config folder. These are available for RSV-A and RSV-B.
