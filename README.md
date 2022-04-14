# eqtls.model-nf

A Nextflow pipeline for predicting tissue-active eQTLs from chromatin features.

The pipeline performs the following analysis steps:

*
*
*

The pipeline uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.

## Requirements

- Unix-like operating system (Linux, MacOS, etc.)
- Java 8 or later 

## Pipeline usage

Launching the pipeline with the `--help` parameter shows the help message:

```
nextflow run eQTLs.model.nf --help
```

```
N E X T F L O W  ~  version 20.10.0
Launching `eQTLs.model.nf` [reverent_heisenberg] - revision: c5d78a089d

eQTLs.model-nf: A pipeline for predicting tissue-active eQTLs.
==============================================================================================
The pipeline takes as input one-gene eQTLs from a donor tissue and predicts which of them are active in a target tissue.

Usage:
nextflow run eQTLs.model.nf [options]

Parameters:
--index			    	Index file containing target tissue info.
--outFolder	            	Output directory (default: results).
--dt			    	Donor tissue providing the list of eQTLs to be predicted in a target tissue.
--eqtls_dt		    	BED file containing all eQTLs in the donor tissue (chrom, start, end, tissue). NOTE: coordinates are 0-based.
--eqtls_slope_distance_dt    	File containing SNP, gene_id, tss_distance, slope for every eQTL-gene pair in donor tissue.
--eqtls_oneGene_dt           	File containing one-gene eQTLs in the donor tissue (chrom_start_end). NOTE: coordinates are 0-based.
--exp_list		    	File containing a list of all EN-TEx functional genomics experiments used for the predictions (bigBed file_id, target, tissue). bigBed = peaks file; target = histone mark / TF / ATAC / DNase.
--bigbed_folder              	Folder containing the files listed in "exp_list".
--entex_rnaseq_m             	Matrix of TPM values for genes (rows) across samples (columns). NOTE: the file is gzipped.
--TSSs			    	BED file containing all non-redundant TSSs +/- 2 Kb (if two isoforms have the same TSS, it will be counted only once). Chrom, start, end, transcript_id, placeholder, strand, gene_id.
--cCREs			    	BED file containing GRCh38 cCREs from ENCODE3.
--repeats             	    	BED file containing repeated elements in GRCh38. NOTE: the file is gzipped.
--keep_only_tested_snps      	Whether the prediction should be restrited to SNPs tested in the taregt tissue by GTEx (default: false).
```

## Input files and format


## Pipeline results
