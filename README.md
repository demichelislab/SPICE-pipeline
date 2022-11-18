<img src="/images/SPICE_COLORI.png" width=50%/>

# SPICE-pipeline
Framework for allele specific analysis of matched tumor and normal next generation sequencing data

## Pre-computed data

All precomputed TCGA allele-specific data, genomic variants data (SNVs and
indels) and association of expression to LOH status data are available [here](https://www.dropbox.com/sh/jzk9roi8c238cc0/AABTTlO7YTneLrjyVUIV8vJia?dl=0) and from the official Zenodo repository (https://zenodo.org/record/5266542#.Y3d-9-zMKAx).


### Files descriptions

The following sections describe the content of each file.

#### Genomic data

The genomic data table contain the corrected copy number values and the allele
specific copy number for each gene for each sample analyzed in the study.

In the table below each field of the results is described:

| Field                        | Description                                                           |
| --------                     | --------                                                              |
| dataset                      | The name of the dataset                                               |
| sample_id                    | The ID of the sample in this study                                    |
| hugo                         | The HUGO symbol of the gene                                           |
| log2                         | The raw log2 value                                                    |
| log2_int                     | The discretized log2 value                                            |
| log2_corr                    | The corrected log2 value                                              |
| log2_corr_int                | The corrected and discretized log2 value                              |
| as_cn_disc                   | The allele specific discretized copy number                           |
| count_snvs_deleterious       | The number of SNVs annotated as deleterious present in the gene       |
| count_insertions_deleterious | The number of insertions annotated as deleterious present in the gene |
| count_deletions_deleterious  | The number of deletions annotated as deleterious present in the gene  |
| count_snvs                   | The total number of SNVs present in the gene                          |
| count_insertions             | The total number of insertions present in the gene                    |
| count_deletions              | The total number of deletions present in the gene                     |
| segment_id                   | The ID of the segment containing the gene                             |

For a detailed description of log2_corr, log2_corr_int, as_cn_disc and allele specific analysis
please see: https://doi.org/10.1002/cpbi.81

#### SNVs data

The SNV data table contains the information about SNVs that were called in the
TCGA samples

In the table below the most important fields in the table are described. For
other fields see [VEP documentation](https://m.ensembl.org/info/docs/tools/vep/script/vep_options.html).

| Field              | Description                                                                                |
| --------           | --------                                                                                   |
| sample_id          | The ID of the sample in this study                                                         |
| uid                | The ID to identify the SNV                                                                 |
| uploaded_variation | The ID of the variant if present in dbSNP                                                  |
| chromosome         | The chromosome where the SNV is located                                                    |
| position           | The nucleotide position where the SNV is located along the chromosome                      |
| allele             | The alterative allele present at the SNV coordinate                                        |
| ...                | See [VEP documentation](https://m.ensembl.org/info/docs/tools/vep/script/vep_options.html) |
| rc_ref_normal      | The count of high quality reads supporting the reference allele in the normal sample       |
| rc_alt_normal      | The count of high quality reads supporting the alternative allele in the normal sample     |
| af_normal          | The allelic fraction of the alternative allele in the normal sample                        |
| cov_normal         | The coverage at the location of the SNV in the normal sample                               |
| rc_ref_tumor       | The count of high quality reads supporting the reference allele in the tumor sample        |
| rc_alt_tumor       | The count of high quality reads supporting the alternative allele in the tumor sample      |
| af_tumor           | The allelic fraction of the alternative allele in the tumor sample                         |
| cov_tumor          | The coverage at the location of the SNV in the tumor sample                                |
| t_cov              | The coverage at the location of the SNV in the tumor sample                                |
| t_af               | The allelic fraction of the alternative allele in the tumor sample                         |
| n_admreads         | The number of reads coming from the admixture at the SNV location                          |
| t_ref_count_corr   | The number of reads supporting the reference allele corrected for admixture                |
| t_af_corr          | The corrected allelic fraction of the SNV                                                  |
| cn_int             | The integer copy number of the segment containing the SNV                                  |
| cn_snvmut          | The copy number of the SNV inferred based on the corrected SNV AF                          |
| vafexp             | The allelic fraction of the SNV based on cn_snvmut                                         |
| snv_clonality      | The clonality of the SNV computed based on the copy number                                 |
| snv_clonality_int  | The clonality of the SNV computed based on the rounded (integer) copy number               |
| study              | The study that contains the variant                                                        |

#### Indels data

The indels data table contains the information about indels that were called in
the TCGA samples

In the table below the most important fields in the table are described. For
other fields see [VEP documentation](https://m.ensembl.org/info/docs/tools/vep/script/vep_options.html).

| Field              | Description                                                                                |
| --------           | --------                                                                                   |
| sample_id          | The ID of the sample in this study                                                         |
| uid                | The ID to identify the variant                                                             |
| uploaded_variation | The ID of the variant if present in dbSNP                                                  |
| chromosome         | The chromosome where the variant is located                                                |
| start              | The nucleotide position where the variant starts                                           |
| end                | The nucleotide position where the variant ends                                             |
| allele             | The alterative variant present at the locus                                                |
| ...                | See [VEP documentation](https://m.ensembl.org/info/docs/tools/vep/script/vep_options.html) |
| rc_ref_normal      | Only valid for SNVs                                                                        |
| rc_alt_normal      | Only valid for SNVs                                                                        |
| af_normal          | Only valid for SNVs                                                                        |
| cov_normal         | Only valid for SNVs                                                                        |
| rc_ref_tumor       | Only valid for SNVs                                                                        |
| rc_alt_tumor       | Only valid for SNVs                                                                        |
| af_tumor           | Only valid for SNVs                                                                        |
| cov_tumor          | Only valid for SNVs                                                                        |
| t_cov              | Only valid for SNVs                                                                        |
| t_af               | Only valid for SNVs                                                                        |
| n_admreads         | Only valid for SNVs                                                                        |
| t_ref_count_corr   | Only valid for SNVs                                                                        |
| t_af_corr          | Only valid for SNVs                                                                        |
| cn_int             | Only valid for SNVs                                                                        |
| cn_snvmut          | Only valid for SNVs                                                                        |
| vafexp             | Only valid for SNVs                                                                        |
| snv_clonality      | Only valid for SNVs                                                                        |
| snv_clonality_int  | Only valid for SNVs                                                                        |
| study              | The study that contains the variant                                                        |

#### LOH/CN linear model results

These tables contain the results from linear models associating copy number (CN) and loss of heterozygosity (LOH) to gene expression.

In the table below the field present in the table are described.

| Field         | Description                                                      |
| --------      | --------                                                         |
| model_p       | p value of the model                                             |
| beta_cntot    | beta value for the variable "CN"                                 |
| stde_cntot    | standard error for the variable "CN"                             |
| t_val_cntot   | t value for the variable "CN"                                    |
| pval_cntot    | p value for the variable "CN"                                    |
| beta_loh      | beta value for the variable "LOH"                                |
| stde_loh      | standard error for the variable "LOH"                            |
| tval_loh      | t value for the variable "LOH"                                   |
| pval_loh      | p value for the variable "LOH"                                   |
| gene          | gene symbol of the tested gene                                   |
| class         | class of the gene ("TSG","OG","ESSENTIAL","OTHER")               |
| fdr_model     | fdr of the model                                                 |
| fdr_loh       | fdr of the variable "LOH"                                        |
| fdr_cntot     | fdr of the variable "CN"                                         |
| ass_coeff_CN  | coefficient of association (beta/stde) for the variable "CN"     |
| ass_coeff_LOH | coefficient of association (beta/stde) for the variable "LOH"    |

## Code

The code folder contains the code for the pipeline.

### Pipeline

The pipeline folder contains the complete code of analysis pipeline described
in the paper.
A CWL version of the pipeline is available here: https://github.com/demichelislab/SPICE-pipeline-CWL

#### Running the pipeline

In order to run the pipeline some resources files and tools are needed. A
complete example folder that includes both the code and the resources can be
downloaded from [here](https://www.dropbox.com/s/3g7ss1xhhgivfg7/pipeline_example.tar.gz?dl=0).
Be aware that the file is big (12.5GB). The complete instructions needed to
run the pipeline are inside the archive.

## Fundings

This project is funded by (ERC-CoG-2014-648670)
