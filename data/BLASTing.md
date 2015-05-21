This folder contains subfolders with BLAST output results.

The folder names follow the format: ``query-database``

``query``: Denotes the data that was BLASTed
	- ``all_orfs``: all predicted open reading frames
	- ``all_genes``: all genes (subset of ORFs that have gene calling information).
	- ``has_cog``: subset of genes that have COG data
	- ``no_cog``: subset of genes without COG data

``database``: Database that the sequences were BLASTed against.
	- ``nr``: the NCBI non-redundant protein database
	- ``eggNOG``: the eggNOGv4 database
	- ``COG_NCBI``: the NCBI COG database

See ``BLAST_template.sh`` for a template bash script to automate BLASTing.
