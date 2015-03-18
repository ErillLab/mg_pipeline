mg_pipeline
===========

This repository provides a set of tools to help analyze transcriptional regulatory networks in metagenomic data.

The main approach of this project is to predict transcription factor binding sites in gene promoters to infer their regulation. From those data we can then proceed to test hypotheses about regulatory interactions between genes and across taxa.

The mg-rast directory contains a pipeline for processing metagenomic data from [MG-RAST](https://metagenomics.anl.gov).

Currently we are using the [Integrated Gene Catalog (IGC)](http://meta.genomics.cn/metagene/meta/home) (see [Li et al., 2014](http://www.nature.com/nbt/journal/v32/n8/full/nbt.2942.html)) as our main dataset.

See `igc_pipeline.py` for the main processing code.