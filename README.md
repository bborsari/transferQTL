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
The pipeline takes as input one-gene eQTLs from a donor tissue and predicts which of them are active in one or more target tissues.

Usage:
nextflow run eQTLs.model.nf [options]

Parameters:
--dt			    	Donor tissue providing the catalog of eQTLs to be predicted in the target tissue(s).
--eqtls_dt		    	Donor-tissue eQTLs.
--eqtls_slope_distance_dt    	Slope and TSS-distance of eQTLs in donor tissue.
--eqtls_oneGene_dt           	Donor-tissue eQTLs linked to only one gene.
--index                         Index file containing target tissue info.
--exp_list		    	File containing a list of all EN-TEx functional genomics experiments used for the predictions (bigBed file_id, target, tissue). bigBed = peaks file; target = histone mark / TF / ATAC / DNase.
--bigbed_folder              	Folder containing the files listed in "exp_list".
--entex_rnaseq_m             	Matrix of TPM values for genes (rows) across samples (columns). NOTE: the file is gzipped.
--TSSs			    	BED file containing all non-redundant TSSs +/- 2 Kb (if two isoforms have the same TSS, it will be counted only once). Chrom, start, end, transcript_id, placeholder, strand, gene_id.
--cCREs			    	BED file containing GRCh38 cCREs from ENCODE3.
--repeats             	    	BED file containing repeated elements in GRCh38. NOTE: the file is gzipped.
--keep_only_tested_snps      	Whether the prediction should be restrited to SNPs tested in the taregt tissue by GTEx (default: false).
--outFolder                     Output directory (default: results).
```

## Input data and files format

`eQTLs.model-nf` requires the following input data:

* **Donor tissue** (option `--dt`). The name of the tissue used as donor tissue. 

* **eQTLs from donor tissue** (option `--eqtls_dt`). BED file containing all eQTLs reported in the donor tissue. See example below (coordinates are zero-based):

```
chr1	64763	64764	Lung
```

* **Slope and TSS-distance of eQTLs in donor tissue** (option `--eqtls_slope_distance_dt`). This info is, for instance, reported by GTEx. Z-score can be an alternative measure to the slope. See example below:

```
SNP			gene_id			tss_distance	slope
chr1_64763_64764	ENSG00000227232.5	35211		0.370865
```

* **One-gene eQTLs from donor tissue** (option `--eqtls_oneGene_dt`). One-column file contaning SNP_id (chrom_start_end) of donor-tissue eQTLs linked to only one gene. Our model predicts only one-gene eQTLs, for the time being. 

* **Index file** (option `--index`). tsv file containing relevant info for the target tissue(s). Here is an example of the file format:

```
target_tissue1	/path/to/target_tissue1.eQTLs.bed	/path/to/target_tissue1.TestedSNPs.bed	/path/to/target_tissue1.ChromatinSignalTable.tsv
target_tissue2	/path/to/target_tissue2.eQTLs.bed	/path/to/target_tissue2.TestedSNPs.bed	/path/to/target_tissue2.ChromatinSignalTable.tsv
target_tissue3  /path/to/target_tissue3.eQTLs.bed       /path/to/target_tissue3.TestedSNPs.bed  /path/to/target_tissue3.ChromatinSignalTable.tsv
```

The fields in the file correspond to:

1. Target tissue. We will predict which donor-tissue eQTLs are active in a given target tissue. Multiple models can be trained in parallel on multiple target tissues. 
2. Path to BED file containing eQTLs in the target tissue (in our case we use GTEx eQTLs). This file contains chrom, start, end, tissue. NOTE: coordinates are 0-based.
3. Path to BED file containing all SNPs tested in the target tissue (in our case we consider SNPs tested by GTEx). This file contains chrom, start, end, tissue. NOTE: coordinates are 0-based.
4. Path to tsv file containing, for every chromatin feature available in the target tissue, the fold-change signal around the SNPs being predicted. This corresponds to the average fold-change signal in a +/- 5 Kb window centered at the SNP. 

* **Chromatin experiments in target tissue(s)** (option `--exp_list`). tsv file containing a list of peak-calling files used for the predictions. These files correspond to experiments (ChIP/ATAC/DNase-seq) performed in the target tissue(s). The content of this file is: bigBed peaks file_id, target 



## Pipeline results
