Generate bacterial species trees through alignments of 92 core genes.

Dependencies:
* ASTRAL (https://github.com/smirarab/ASTRAL)
* RAxML (https://github.com/stamatak/standard-RAxML)
* MAFFT (http://mafft.cbrc.jp/alignment/software/)
* HMMER3 (http://hmmer.org/)
* BioPerl (http://bioperl.org/index.html)
* prodigal (http://prodigal.ornl.gov/)

Install:
No formal install is required, but paths will need to be updated in core_species_tree.pl

Usage:
perl core_species_tree.pl genomes/*.fna

Output:
* astral.species.tre -- species tree
* astral.species.bs -- astral bootstrap trees
* astral.species.consensus -- consensus bootstrap tree

Workflow:
* Call genes with prodigal
* hmmsearch genes for TIGRFAM hits
* Gene family protein alignment
* Convert to nucleotide
* RAxML gene trees
* ASTRAL species tree