/*
 * Copyright (c) 2021, Beatrice Borsari
 *
 * This file is part of 'eQTLs.model-nf':
 * A Nextflow pipeline for predicting tissue-active eQTLs from chromatin features.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */



// Print usage

if (params.help) {
    log.info ''
    log.info 'eQTLs.model-nf: A Nextflow pipeline for predicting tissue-active eQTLs.'
    log.info '=============================================================================================='
    log.info 'The pipeline takes as input eQTLs from a donor tissue and predicts which of them are active in one or more target tissues.'
    log.info ' '
    log.info 'Usage: '
    log.info '    nextflow run eQTLs.model.nf [options]'
    log.info ''
    log.info 'Parameters:'
    log.info ' --dt			    Donor tissue.'
    log.info ' --eqtls_dt		    Donor-tissue eQTLs.'
    log.info ' --eqtls_slope_distance_dt    Slope and TSS-distance of donor-tissue eQTLs.'
    log.info ' --eqtls_oneGene_dt           Donor-tissue eQTLs linked to only one gene.'
    log.info ' --exp_list		    Functional genomics experiments used for feature extraction.'
    log.info ' --bigbed_folder              Directory containing peak-calling files.'
    log.info ' --entex_rnaseq_m             EN-TEx gene expression matrix.'
    log.info ' --TSSs			    List of annotated TSSs.'
    log.info ' --cCREs			    GRCh38 cCREs from ENCODE3.'
    log.info ' --repeats             	    GRCh38 repeats.'
    log.info ' --keep_only_tested_snps      Whether to restrict the prediction to donor-tissue eQTLs tested in the target tissue (default: false).'
    log.info ' --index			    Index file with target tissue info.'
    log.info ' --outFolder		    Output directory (default: results).'
    log.info ''
    exit(1)
}


/* ~~~~~~~~~~~~~~~
 * BEGIN
 * ~~~~~~~~~~~~~~~
 */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains:
 * 1. donor_tissue eQTLs
 * 2. folder w/ bigBed peak files for feature extraction
 */
Channel.of(file(params.eqtls_dt), params.bigbed_folder)
.collect()
.set{start_ch}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create channels that contain: 
 * (a): target_tissue, target_tissue eQTLs (eqtls_tt), 
 * (b): target_tissue, target_tissue tested snps (tested_tt)
 * (c): target_tissue, table of chromatin features' signal around the SNV (signal_table_tt)
 */
 
index = file(params.index)

Channel.from(index.readLines())
.map { line ->
  def list = line.tokenize()
  def tt = list[0]
  def eqtls_tt = resolveFile(list[1], index)
  def tested_tt = resolveFile(list[2], index)
  def signal_table_tt = resolveFile(list[3], index)
  [ tt, eqtls_tt, tested_tt, signal_table_tt ]
}
.multiMap { it ->
     a: it[0, 1]
     b: it[0, 2]
     c: it[0, 3]
}.set{tt_ch}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains,
 * for every bigBed file,
 * 1. file_id
 * 2. assay_name (either histone mark, TF, ATAC, DNase)
 * 3. tissue in which the experiment was performed (can include also dt)
 */
exps = file(params.exp_list)

Channel.from(exps.readLines())
.map { line ->
  def list = line.tokenize()
  def file_id = list[0]
  def assay_name = list[1]
  def tissue = list[2]
  [ file_id, assay_name, tissue ]
}.set{bigbed_ch}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains,
 * for every bigBed file,
 * 1. donor_tissue eQTLs
 * 2. bigBed file
 * 3. file_id
 * 4. assay_name (either histone mark, TF, ATAC, DNase)
 * 5. tissue
 */
start_ch
.combine(bigbed_ch)
.map { it ->
  [ it[0], file("$baseDir/"+it[1]+"/"+it[2]+".bigBed"), it[2], it[3], it[4] ]
}
.set{bedtools_ch}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #1: Intersect donor-tissue eQTLs
 * w/ peaks of histone marks / TFs / ATAC and DNase
 * in every specified target tissue
 */
process bedtools_intersect {

  publishDir "${params.outFolder}/bedfiles", overwrite: false

  conda '/users/rg/bborsari/.conda/envs/ENCODE_RC/'

  input:
  tuple file(eqtls_dt), file(bigbed), file_id, assay_name, tissue from bedtools_ch

  output:
  tuple assay_name, tissue, file("${file_id}.bed") into aggregate_ch

  script:

  """
  # bedtools intersect
  intersect.peaks.cosi.regions.sh --regions $eqtls_dt --peaks $bigbed --outFile ${file_id}.bed

  # keep only eQTLs with an intersection
  awk '\$NF>0{print \$1"_"\$2"_"\$3}' ${file_id}.bed | sort -u > ${file_id}.bed.tmp 
  mv ${file_id}.bed.tmp ${file_id}.bed
  """

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains,
 * for a given assay and tissue,
 * 1. donor_tissue eQTLs
 * 2. assay_name
 * 3. tissue
 * 4. tables of presence/absence of peaks
 */
Channel.of(file(params.eqtls_dt))
.collect()
.set{seed1_ch}

seed1_ch
.combine(aggregate_ch)
.groupTuple(by: [0, 1, 2])
.set{ready2aggregate_ch}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #2: Aggregate experiments performed across multiple donors
 * for the same assay (histone mark, TF, ATAC or DNase) and tissue 
 * aka: for a given snp, 
 * if there is a peak in >= 1 donor: 1
 * else: 0
 */
process aggregate_exp {

  publishDir "${params.outFolder}/agg.tables", overwrite: false

  input:
  tuple file(eqtls_dt), assay, tissue, file(bed) from ready2aggregate_ch

  output:
  tuple assay, file("${tissue}.${assay}.tsv") into tissue_assay_agg_ch_a
  tuple tissue, file("${tissue}.${assay}.tsv") into tissue_assay_agg_ch_b

  script:
  """
  awk 'BEGIN{FS="\t";OFS="_"}{print \$1,\$2,\$3}' $eqtls_dt > ${tissue}.${assay}.tsv

  for f in ${bed}; do
	if [ -s \$f ];
	then
		join.py -b ${tissue}.${assay}.tsv -a <(awk '{print \$0"\t"1}' \$f) -u -p 0 > ${tissue}.${assay}.tmp;
		mv ${tissue}.${assay}.tmp ${tissue}.${assay}.tsv
	fi 
  done

  awk 'BEGIN{FS=OFS="\t"}{n=0; for (i=2;i<=NF;i++){if (\$i>0){n+=1}}; if (n>0){print \$1, 1} else {print \$1, 0}}' ${tissue}.${assay}.tsv > ${tissue}.${assay}.tmp
  mv ${tissue}.${assay}.tmp ${tissue}.${assay}.tsv
  """

}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains,
 * for a given assay,
 * 1. donor_tissue eQTLs
 * 2. tables of presence/absence of peaks for every tissue
 */
Channel.of(file(params.eqtls_dt))
.collect()
.set{seed2_ch}

seed2_ch
.combine(tissue_assay_agg_ch_a)
.groupTuple(by: [0, 1])
.set{ready2build_byAssay_summary_table_ch}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #3: Prepare a summary table by assay
 * aka for each assay register presence/absence 
 * of peaks in every profiled tissue
 */
process build_table_byTarget {

  publishDir "${params.outFolder}/assay.summary.tables", overwrite: false

  input:
  tuple file(eqtls_dt), assay, file(tsv) from ready2build_byAssay_summary_table_ch

  output:
  file("${assay}.summary.table.tsv") into byAssay_summary_table_ch

  script:
  """
  awk 'BEGIN{FS="\t";OFS="_"}{print \$1,\$2,\$3}' $eqtls_dt > ${assay}.summary.table.tsv

  for f in ${tsv}; do
	if [ -s \$f ];
	then
		tissue="\$(basename \$f | awk '{split(\$1, a, "."); print a[1]}')"
		join.py -b ${assay}.summary.table.tsv -a \$f |\
		sed "1i\${tissue}" > ${assay}.summary.table.tmp; 
		mv ${assay}.summary.table.tmp ${assay}.summary.table.tsv
	fi
  done

  awk 'BEGIN{FS=OFS="\t"}NR>1{a=0;n=(NF-1);for (i=2;i<=NF;i++){if (\$i>0){a+=1}}; print \$1, a/n}' ${assay}.summary.table.tsv |\
  sed "1i$assay" > ${assay}.summary.table.tmp
  mv ${assay}.summary.table.tmp ${assay}.summary.table.tsv
  """
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains
 * all summary tables by assay
 */
byAssay_summary_table_ch
.collect()
.groupTuple()
.flatten()
.collect()
.set{merge_byAssay_summary_table_ch}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #4: Merge all summary tables by 
 * assay into a single file
 */
process merge_byAssay_summary_table {

  publishDir "${params.outFolder}/assay.summary.tables", overwrite: false

  input:
  file merge_byAssay_summary_table_ch

  output:
  file("assays.summary.table.tsv") into assays_table_ch

  script:
  """
  paste <(sed 's/H3K27ac/SNP\tH3K27ac/' H3K27ac.summary.table.tsv) <(cut -f2 H3K4me3.summary.table.tsv) <(cut -f2 H3K4me1.summary.table.tsv) <(cut -f2 H3K27me3.summary.table.tsv) <(cut -f2 H3K36me3.summary.table.tsv) <(cut -f2 H3K9me3.summary.table.tsv) <(cut -f2 CTCF.summary.table.tsv) <(cut -f2 POLR2A.summary.table.tsv) <(cut -f2 POLR2AphosphoS5.summary.table.tsv) <(cut -f2 EP300.summary.table.tsv) <(cut -f2 ATAC.summary.table.tsv) <(cut -f2 DNase.summary.table.tsv) > assays.summary.table.tsv
  """

}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains,
 * for every tissue that is not the donor_tissue,
 * 1. target_tissue
 * 2. target_tissue eQTLs
 * 3. target_tissue tables of presence/absence peaks for every assay_name
 * 4. donor_tissue eQTLs
 */
Channel.of(file(params.eqtls_dt))
.collect()
.set{seed3_ch}

tissue_assay_agg_ch_b
.groupTuple(by: [0], sort: {it.name})
.filter({it[0] != params.dt})
.concat(tt_ch.a)
.groupTuple(by: [0])
.map { line ->
  def list = line
  def tt = list[0]
  def tables = list[1][0]
  def eqtls = list[1][1]
  [ tt, eqtls, tables ]
}
.combine(seed3_ch)
.set{ready2build_byTissue_summary_table_ch}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #5: Prepare a summary table by tissue
 * aka for each tissue register presence/absence
 * of peaks for every assay
 */
process build_table_byTissue {

  publishDir "${params.outFolder}/tissue.summary.tables", overwrite: false

  input:
  tuple tissue, file(eqtls_tt), file(tsv), file(eqtls_dt) from ready2build_byTissue_summary_table_ch

  output:
  tuple tissue, file("${tissue}.summary.table.peaks.tsv") into byTissue_summary_table_ch

  script:
  """
  # annotate whether an eQTL in the donor tissue is also eQTL in the target tissue
  join.py -b <(awk 'BEGIN{FS="\t";OFS="_"}{print \$1,\$2,\$3}' $eqtls_dt) -a <(awk '{print \$1"_"\$2"_"\$3"\t"1}' $eqtls_tt) -u -p 0 |\
  sed '1iSNP\tis_eQTL' > ${tissue}.summary.table.peaks.tsv

  # for each assay performed in a given tissue
  # concatenate the binary table of presence/absence of peaks
  for f in ${tsv}; do
        if [ -s \$f ];
        then
		assay="\$(basename \$f | awk '{split(\$1, a, "."); print a[2]}')"
		~/bin/join.py -b ${tissue}.summary.table.peaks.tsv -a \$f --b_header |\
		sed "1s|V1|\${assay}|" > ${tissue}.summary.table.peaks.tmp; 
		mv ${tissue}.summary.table.peaks.tmp ${tissue}.summary.table.peaks.tsv;			
	fi
  done
  """
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains:
 * 1. donor-tissue eQTLs that are associated to only one gene
 * 2. slope and tss_distance for donor-tissue eQTLs
 * 3. entex rnaseq matrix 
 */
Channel.of(file(params.eqtls_oneGene_dt), file(params.eqtls_slope_distance_dt), file(params.entex_rnaseq_m))
.collect()
.set{additional_features_ch_a}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains:
 * 1. donor-tissue eQTLs
 * 2. bed file of repeated regions
 * 3. bed file of ENCODE cCREs
 * 4. bed file of annotated TSSs (gencode v24)
 */
Channel.of(file(params.eqtls_dt), file(params.repeats), file(params.cCREs), file(params.TSSs))
.collect()
.set{ready2build_additional_features_ch_b}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #6: Prepare table of donor_tissue eQTLs that
 * 1. are TSS-proximal (aka they intersect annotated TSSs)
 * 2. are outside repeated regions
 * 3. are inside ENCODE cCREs
 */
process get_additional_features_b {

  publishDir "${params.outFolder}/additional.features", overwrite: false

  conda '/users/rg/bborsari/.conda/envs/ENCODE_RC/'

  input:
  tuple file(eqtls_dt), file(repeats), file(cCREs), file(TSSs) from ready2build_additional_features_ch_b

  output:
  tuple file("out.of.repeat.eqtls.tsv"), file("in.cCREs.eqtls.tsv"), file("proximal.eqtls.tsv") into additional_features_ch_b

  script:
  """
  # donor_tissue eQTLs that are proximal to annotated TSSs
  bedtools intersect -a $eqtls_dt -b $TSSs -u |\
  awk '{print \$1"_"\$2"_"\$3"\t"1}' > proximal.eqtls.tsv

  # donor_tissue eQTLs that are outside of repeated regions
  bedtools intersect -a $eqtls_dt -b <(zcat $repeats) -v |\
  awk '{print \$1"_"\$2"_"\$3"\t"1}' > out.of.repeat.eqtls.tsv

  # donor_tissue eQTLs that are inside ENCODE cCREs
  bedtools intersect -a $eqtls_dt -b $cCREs -u |\
  awk '{print \$1"_"\$2"_"\$3"\t"1}' > in.cCREs.eqtls.tsv
  """

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Create a channel that contains:
 * 1. target tissue
 * 2. presence/absence of peaks summary table
 * 3. list of tested snps in target tissue
 * 4. kun's tables of chromatin signal around eqtls
 */
byTissue_summary_table_ch
.concat(tt_ch.b)
.concat(tt_ch.c)
.groupTuple(by: [0])
.map { line ->
  def list = line
  def tt = list[0]
  def summary_table_tt = list[1][0]
  def tested_snps_tt = list[1][1]
  def signal_table_tt = list[1][2]
  [ tt, summary_table_tt, tested_snps_tt, signal_table_tt ]
}
.set{ready2build_input4model_ch}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #7: Prepare input table to run predictive model
 */
if (params.keep_only_tested_snps) {

  process build_input4model_1 {

  publishDir "${params.outFolder}/input.tables", overwrite: false

  input:
  tuple tissue, file(summary_table_tt), file(tested_snps_tt), file(signal_table_tt) from ready2build_input4model_ch
  tuple file(oneGene_eqtls_dt), file(slope_distance_dt), file(entex_rnaseq_m) from additional_features_ch_a
  tuple file(out_repeats_eqtls_dt), file(in_ccres_eqtls_dt), file(proximal_eqtls_dt) from additional_features_ch_b
  file(assays_summary_table) from assays_table_ch

  output:
  tuple tissue, file("${tissue}.input4model.tsv") into input4model_ch

  script:
  """
  input4model.R --target_tissue $tissue --input_matrix $summary_table_tt --outFile ${tissue}.input4model.tsv --keep_only_tested TRUE --tested_snps_tt $tested_snps_tt --one_gene_eqtls $oneGene_eqtls_dt --entex_rnaseq_m $entex_rnaseq_m --slope_distance $slope_distance_dt --out_repeats_eqtls $out_repeats_eqtls_dt --in_cCREs_eqtls $in_ccres_eqtls_dt --proximal_eqtls $proximal_eqtls_dt --signal_tables $signal_table_tt --proportion_marked_tissues $assays_summary_table
  """
  }

}

else {
  
  process build_input4model_2 {

  publishDir "${params.outFolder}/input.tables", overwrite: false

  input:
  tuple tissue, file(summary_table_tt), file(tested_snps_tt), file(signal_table_tt) from ready2build_input4model_ch
  tuple file(oneGene_eqtls_dt), file(slope_distance_dt), file(entex_rnaseq_m) from additional_features_ch_a
  tuple file(out_repeats_eqtls_dt), file(in_ccres_eqtls_dt), file(proximal_eqtls_dt) from additional_features_ch_b
  file(assays_summary_table) from assays_table_ch

  output:
  tuple tissue, file("${tissue}.input4model.tsv") into input4model_ch
  
  script:
  """
  input4model.R --target_tissue $tissue --input_matrix $summary_table_tt --outFile ${tissue}.input4model.tsv --keep_only_tested FALSE --tested_snps_tt $tested_snps_tt --one_gene_eqtls $oneGene_eqtls_dt --entex_rnaseq_m $entex_rnaseq_m --slope_distance $slope_distance_dt --out_repeats_eqtls $out_repeats_eqtls_dt --in_cCREs_eqtls $in_ccres_eqtls_dt --proximal_eqtls $proximal_eqtls_dt --signal_tables $signal_table_tt --proportion_marked_tissues $assays_summary_table
  """

  }

}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Process #8: Run predictive model
 */
process run_model {

publishDir "${params.outFolder}/output.objs", overwrite: false

input:
tuple tissue, file("${tissue}.input4model.tsv") from input4model_ch

output:
tuple tissue, file("${tissue}.RData") into output4model_ch
  
script:
"""
rf.model.R -i ${tissue}.input4model.tsv -t $tissue -p 0.7 -m rf
"""

}




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FUNCTIONS
/*
 * Given a string path resolve it against the index file location.
 * Params:
 * - str: a string value represting the file path to be resolved
 * - index: path location against which relative paths need to be resolved
 */
def resolveFile( str, index ) {
  if( str.startsWith('/') || str =~ /^[\w\d]*:\// ) {
    return file(str)
  }
  else if( index instanceof Path ) {
    return index.parent.resolve(str)
  }
  else {
    return file(str)
  }
}
