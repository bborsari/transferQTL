# transferQTL-nf

A Nextflow pipeline for predicting tissue-active eQTLs from chromatin features.

The pipeline performs the following analysis steps:

* Feature extraction based on EN-TEx functional genomics assays. 
* Prediction of donor-tissue eQTLs active in a target tissue (using Random Forest).
For more details see [pipeline results](https://github.com/bborsari/transferQTL#pipeline-results).

![](https://github.com/bborsari/transferQTL/blob/main/wiki.images/schema.png)

The pipeline uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.

## Requirements

- Unix-like operating system (Linux, MacOS, etc.)
- Java 8 or later 

## Pipeline usage

Launching the pipeline with the `--help` parameter shows the help message:

```
nextflow run transferQTL.nf --help
```

```
N E X T F L O W  ~  version 20.10.0
Launching `transferQTL.nf` [reverent_heisenberg] - revision: c5d78a089d

eQTLs.model-nf: A Nextflow pipeline for predicting tissue-active eQTLs.
==============================================================================================
The pipeline takes as input eQTLs from a donor tissue and predicts which of them are active in one or more target tissues.

Usage:
nextflow run transferQTL.nf [options]

Parameters:
--dt			    	Donor tissue.
--eqtls_dt		    	Donor-tissue eQTLs.
--eqtls_slope_distance_dt    	Slope and TSS-distance of donor-tissue eQTLs.
--eqtls_oneGene_dt           	Donor-tissue eQTLs linked to only one gene (eQTLs with multiple eGenes are filtered out).
--exp_list		    	Functional genomics experiments used for feature extraction.
--bigbed_folder              	Directory containing peak-calling files.
--entex_rnaseq_m             	EN-TEx gene expression matrix.
--TSSs			    	List of annotated TSSs.
--cCREs			    	GRCh38 cCREs from ENCODE3.
--repeats             	    	GRCh38 repeats.
--keep_only_tested_snps      	Whether to restrict the prediction to donor-tissue eQTLs tested in the target tissue (default: false).
--index                         Index file with target tissue info.
--outFolder                     Output directory (default: results).
```

## Input data and files format

`transferQTL.nf` requires the following input data:

* **Donor tissue** (`--dt`). The tissue providing the source catalog of eQTLs. The pipeline will predict whether these eQTLs are active in one or more target tissues. Example: `Lung`

* **Donor-tissue eQTLs** (`--eqtls_dt`). BED file containing the source catalog of donor-tissue eQTLs. See example below (coordinates are zero-based):

```
chr1	64763	64764	Lung
```

* **Slope and TSS-distance of donor-tissue eQTLs** (`--eqtls_slope_distance_dt`). The info in this tsv file is, for instance, reported by GTEx. Z-score can be an alternative measure to the slope. See example below:

```
SNP			gene_id			tss_distance	slope
chr1_64763_64764	ENSG00000227232.5	35211		0.370865
```

* **Donor-tissue eQTLs linked to only one gene** (`--eqtls_oneGene_dt`). One-column file containing a list of eQTLs linked to only one gene in the donor tissue. See example below:

```
chr1_64763_64764
```

* **Functional genomics experiments used for feature extraction** (`--exp_list`). tsv file containing a list of peak-calling bigBed files used to extract chromatin features of donor-tissue eQTLs in the relevant target tissue. These files correspond to EN-TEx functional genomics experiments (ChIP/ATAC/DNase-seq). See example below:

```
ENCFF821QBE     CTCF            Adrenal_Gland
ENCFF945HXI     H3K27ac         Adrenal_Gland
ENCFF459ANZ     POLR2A          Adrenal_Gland
ENCFF036JMQ     H3K4me1         Adrenal_Gland
ENCFF158ICO     ATAC            Adrenal_Gland
```

* **Directory containing peak-calling files** (`--bigbed_folder`). Folder containing the files listed in `--exp_list`.

* **EN-TEx gene expression matrix** (`--entex_rnaseq_m`). [bgzip](http://www.htslib.org/doc/bgzip.html)-compressed matrix of TPM values for genes (rows) across donor_tissue samples (columns). 

* **List of annotated TSSs** (`--TSSs`). BED file containing coordinates of all annotated TSSs +/- 2 Kb. See example below:

```
chrX	100634688	100638689	ENST00000496771.5	0	-	ENSG00000000003.14
```

* **GRCh38 cCREs from ENCODE3** (`--cCREs`). BED file containing GRCh38 cCREs from ENCODE3 (downloaded from http://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.bed).

* **GRCh38 repeats** (`--repeats`). [bgzip](http://www.htslib.org/doc/bgzip.html)-BED file containing repeated elements annotated in GRCh38. The file was downloaded from http://genome.ucsc.edu/cgibin/hgTables, after setting `group` = `repeats` and `track` = `Repeatmasker`.

* **Keep only tested SNPs** (`--keep_only_tested_snps`). Whether to restrict the predictions only to donor-tissue eQTLs tested by GTEx for eQTL-gene association in the target tissue (`default` = `false`). The model performance can only be evaluated with SNPs that were tested for being eQTLs in the target tissue. 

* **Index file with target tissue info** (`--index`). tsv file containing relevant info for the target tissue(s). Here is an example of the file format:

```
Adrenal_Gland	/path/to/Adrenal_Gland.eQTLs.bed	/path/to/Adrenal_Gland.TestedSNPs.bed	/path/to/Adrenal_Glaand.ChromatinSignalTable.tsv
Artery_Aorta	/path/to/Artery_Aorta.eQTLs.bed		/path/to/Artery_Aorta.TestedSNPs.bed	/path/to/Artery_Aorta.ChromatinSignalTable.tsv
Colon_Sigmoid   /path/to/Colon_Sigmoid.eQTLs.bed        /path/to/Colon_Sigmoid.TestedSNPs.bed   /path/to/Colon_Sigmoid.ChromatinSignalTable.tsv
```

The fields in the file correspond to:

1. Target tissue. We will predict which donor-tissue eQTLs are active in a given target tissue. Multiple models can be trained in parallel on multiple target tissues. 
2. Path to BED file containing eQTLs in the target tissue (in our case we use GTEx eQTLs). See example below (coordinates are zero-based):

```
chr1	10000448	10000449	Adrenal_Gland
```

3. Path to BED file containing all SNPs tested in the target tissue (in our case we consider SNPs tested by GTEx). See example below (coordinates are 0-based):

```
chr1	10000043	10000044	Adrenal_Gland
```


4. Path to tsv file containing, for every chromatin feature available in the target tissue, the fold-change signal around the SNPs being predicted. This corresponds to the average fold-change signal in a +/- 5 Kb window centered at the SNP. See example below:

```
SNP			H3K27ac		H3K4me1		H3K9me3		CTCF		POLR2A		ATAC		DNase
chr1_10000043_10000044	0.0502729	0.580827	0.242431	0.0718333	0.136759	0.021416	0.0263053
``` 


## Pipeline results

The output folder can be specified with `--outFolder` (default: `results`).  

* **Feature extraction**. We employ [39 features](https://github.com/bborsari/transferQTL/blob/main/wiki.images/features.pdf) in the predictions.
  * Feature 1 is computed by `Process #7`.
  * Features 2-3 are provided by the user (`--eqtls_slope_distance_dt`).
  * Features 4-15 are computed by `Processes #1,2,5`.
  * Features 17-24 are provided by the user (`--index`).
  * Features 25-35 are computed by `Processes #1-4`.
  * Features 37-39 are computed by `Process #6`.

The pipeline produces a final `SNPs x features` table for every pair of donor-target tissues. These tables are stored inside `results/input.tables`. See example below:

```
SNP                        is_eQTL  ATAC  CTCF  DNase  H3K27ac  H3K27me3  H3K36me3  H3K4me1  H3K4me3  H3K9me3  POLR2A  POLR2AphosphoS5  sum  tss_distance  slope     cv                is_out_repeat  is_cCRE  H3K27ac_p  H3K4me3_p  H3K4me1_p  H3K27me3_p  H3K36me3_p  H3K9me3_p  CTCF_p  POLR2A_p  POLR2AphosphoS5_p  EP300_p  ATAC_p  DNase_p  H3K27ac_k  H3K4me1_k  H3K27me3_k  H3K9me3_k  CTCF_k     POLR2A_k   ATAC_k     DNase_k      is_proximal
chr10_100004826_100004827  y        0     0     0      0        0         0         1        0        0        0       0                1    5092          0.314238  1.27037111031612  1              0        0.296296   0          0.52       0           0           0          0       0         0                  0        0       0        0.137448   1.57745    0.438172    0.568175   0.0815184  0.0446231  0.0269283  3.48443e-05  0
```
* **Model prediction**.
The pipeline trains a random forest model for every pair of donor-target tissues using the [R package caret](https://topepo.github.io/caret/) and by implementing a 5-fold cross-validation schema (`Process #8`). The model predictions are validated against the GTEx eQTL catalog for a given target tissue.
Results of the prediction are stored inside `results/output.objs`.
