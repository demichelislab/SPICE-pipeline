<img src="/images/SPICE_COLORI.png" width=50%/>

# SPICE-pipeline
Framework for allele specific analysis of matched tumor and normal next generation sequencing data

## Pre-computed data

All precomputed TCGA allele-specific data, genomic variants data (SNVs and
indels) and haploinsufficiency data are available [here](https://www.dropbox.com/sh/jzk9roi8c238cc0/AABTTlO7YTneLrjyVUIV8vJia?dl=0)

The table below contains the links to all the output files for each dataset:

 | Dataset  | Genomic data                                                                                 | SNVs                                                                             | Indels                                                                               | Haploinsufficiency data                                                                                         |
| -------- | --------                                                                                     | --------                                                                         | --------                                                                             | --------                                                                                                        |
| ACC      | [ACC genomic data](https://www.dropbox.com/s/gpi8whbdj0z871o/acc_genomic_data.tsv.gz?dl=0)   | [ACC SNVs](https://www.dropbox.com/s/j0p1ojyszmi1k62/acc_snv_data.tsv.gz?dl=0)   | [ACC Indels](https://www.dropbox.com/s/g927som1vyatr8i/acc_indel_data.tsv.gz?dl=0)   | [ACC haploinsufficiency tables](https://www.dropbox.com/s/6zsw6igh7y64qo0/acc_haploins_cnnl_hemi.tsv.gz?dl=0)   |
| BLCA     | [BLCA genomic data](https://www.dropbox.com/s/yma2scxl4ov28az/blca_genomic_data.tsv.gz?dl=0) | [BLCA SNVs](https://www.dropbox.com/s/3375yj50mj1myae/blca_snv_data.tsv.gz?dl=0) | [BLCA Indels](https://www.dropbox.com/s/gepfzgh97dhj2vp/blca_indel_data.tsv.gz?dl=0) | [BLCA haploinsufficiency tables](https://www.dropbox.com/s/k2mueefh7a8qngc/blca_haploins_cnnl_hemi.tsv.gz?dl=0) |
| BRCA     | [BRCA genomic data](https://www.dropbox.com/s/lqfqgzkjn9qbcpg/brca_genomic_data.tsv.gz?dl=0) | [BRCA SNVs](https://www.dropbox.com/s/8rxcwvp1mld21kg/brca_snv_data.tsv.gz?dl=0) | [BRCA Indels](https://www.dropbox.com/s/0sxjaqi7a7jzea4/brca_indel_data.tsv.gz?dl=0) | [BRCA haploinsufficiency tables](https://www.dropbox.com/s/wmrtjkrmguikmo6/brca_haploins_cnnl_hemi.tsv.gz?dl=0) |
| COAD     | [COAD genomic data](https://www.dropbox.com/s/kmft5vgcxiu68el/coad_genomic_data.tsv.gz?dl=0) | [COAD SNVs](https://www.dropbox.com/s/2tiewf4xpu8iirk/coad_snv_data.tsv.gz?dl=0) | [COAD Indels](https://www.dropbox.com/s/kpb1y1k3u8ttllh/coad_indel_data.tsv.gz?dl=0) | [COAD haploinsufficiency tables](https://www.dropbox.com/s/4lfzh2f5a2z0t4x/coad_haploins_cnnl_hemi.tsv.gz?dl=0) |
| GBM      | [GBM genomic data](https://www.dropbox.com/s/86y12xacczwlms2/gbm_genomic_data.tsv.gz?dl=0)   | [GBM SNVs](https://www.dropbox.com/s/clollkam8m9pqb3/gbm_snv_data.tsv.gz?dl=0)   | [GBM Indels](https://www.dropbox.com/s/ej694nqh30y9jyt/gbm_indel_data.tsv.gz?dl=0)   | [GBM haploinsufficiency tables](https://www.dropbox.com/s/q5ub3e261aeyjmq/gbm_haploins_cnnl_hemi.tsv.gz?dl=0)   |
| HNSC     | [HNSC genomic data](https://www.dropbox.com/s/9cegc2pztgqtk25/hnsc_genomic_data.tsv.gz?dl=0) | [HNSC SNVs](https://www.dropbox.com/s/lmiicjlcq93iwx7/hnsc_snv_data.tsv.gz?dl=0) | [HNSC Indels](https://www.dropbox.com/s/ykfrnkm8a94c4cq/hnsc_indel_data.tsv.gz?dl=0) | [HNSC haploinsufficiency tables](https://www.dropbox.com/s/xcncm3x0kmv46j4/hnsc_haploins_cnnl_hemi.tsv.gz?dl=0) |
| KICH     | [KICH genomic data](https://www.dropbox.com/s/r0ih03oz5vg2igh/kich_genomic_data.tsv.gz?dl=0) | [KICH SNVs](https://www.dropbox.com/s/e974ybyjndt9j0m/kich_snv_data.tsv.gz?dl=0) | [KICH Indels](https://www.dropbox.com/s/byz06jpd0ll8u46/kich_indel_data.tsv.gz?dl=0) | [KICH haploinsufficiency tables](https://www.dropbox.com/s/hyb99m4im94gjvf/kich_haploins_cnnl_hemi.tsv.gz?dl=0) |
| KIRC     | [KIRC genomic data](https://www.dropbox.com/s/x73ftij4rrw12tk/kirc_genomic_data.tsv.gz?dl=0) | [KIRC SNVs](https://www.dropbox.com/s/ytpu506mpwn0fp2/kirc_snv_data.tsv.gz?dl=0) | [KIRC Indels](https://www.dropbox.com/s/zkjw217zgk9u334/kirc_indel_data.tsv.gz?dl=0) | [KIRC haploinsufficiency tables](https://www.dropbox.com/s/hm3fqzgedspr84r/kirc_haploins_cnnl_hemi.tsv.gz?dl=0) |
| KIRP     | [KIRP genomic data](https://www.dropbox.com/s/76bfzrle9xvt2v6/kirp_genomic_data.tsv.gz?dl=0) | [KIRP SNVs](https://www.dropbox.com/s/1k563c9lsgwcm4r/kirp_snv_data.tsv.gz?dl=0) | [KIRP Indels](https://www.dropbox.com/s/h2kfa53omdcag03/kirp_indel_data.tsv.gz?dl=0) | [KIRP haploinsufficiency tables](https://www.dropbox.com/s/qvgkwspr2coll1l/kirp_haploins_cnnl_hemi.tsv.gz?dl=0) |
| LAML     | [LAML genomic data](https://www.dropbox.com/s/pv6s2f58ncghz6r/laml_genomic_data.tsv.gz?dl=0) | [LAML SNVs](https://www.dropbox.com/s/2klvrqv38opqot0/laml_snv_data.tsv.gz?dl=0) | [LAML Indels](https://www.dropbox.com/s/9fap2kck94yesbn/laml_indel_data.tsv.gz?dl=0) | [LAML haploinsufficiency tables](https://www.dropbox.com/s/k3bq37fre4uqxef/laml_haploins_cnnl_hemi.tsv.gz?dl=0) |
| LGG      | [LGG genomic data](https://www.dropbox.com/s/vcgn37q1unfk8er/lgg_genomic_data.tsv.gz?dl=0)   | [LGG SNVs](https://www.dropbox.com/s/9l5k4e6pb31dp9f/lgg_snv_data.tsv.gz?dl=0)   | [LGG Indels](https://www.dropbox.com/s/ss5jg0yltjm8j69/lgg_indel_data.tsv.gz?dl=0)   | [LGG haploinsufficiency tables](https://www.dropbox.com/s/270etwmjmgvjagz/lgg_haploins_cnnl_hemi.tsv.gz?dl=0)   |
| LIHC     | [LIHC genomic data](https://www.dropbox.com/s/giduawtpnl7g8je/lihc_genomic_data.tsv.gz?dl=0) | [LIHC SNVs](https://www.dropbox.com/s/15kxsjzlxwrets0/lihc_snv_data.tsv.gz?dl=0) | [LIHC Indels](https://www.dropbox.com/s/ex3m4tpgjkb4vro/lihc_indel_data.tsv.gz?dl=0) | [LIHC haploinsufficiency tables](https://www.dropbox.com/s/012c7dt1xxoiiu7/lihc_haploins_cnnl_hemi.tsv.gz?dl=0) |
| LUAD     | [LUAD genomic data](https://www.dropbox.com/s/fz47h94w7jk8co0/luad_genomic_data.tsv.gz?dl=0) | [LUAD SNVs](https://www.dropbox.com/s/pjpd9qdqrnvtadm/luad_snv_data.tsv.gz?dl=0) | [LUAD Indels](https://www.dropbox.com/s/t13t9e995gx2smi/luad_indel_data.tsv.gz?dl=0) | [LUAD haploinsufficiency tables](https://www.dropbox.com/s/kzso2a3rec8n31r/luad_haploins_cnnl_hemi.tsv.gz?dl=0) |
| LUSC     | [LUSC genomic data](https://www.dropbox.com/s/nb3l136b6a2a68w/lusc_genomic_data.tsv.gz?dl=0) | [LUSC SNVs](https://www.dropbox.com/s/purnuj4qptbwwy5/lusc_snv_data.tsv.gz?dl=0) | [LUSC Indels](https://www.dropbox.com/s/0h3iqlts8cvttgr/lusc_indel_data.tsv.gz?dl=0) | [LUSC haploinsufficiency tables](https://www.dropbox.com/s/qdjj5gx59s66dyr/lusc_haploins_cnnl_hemi.tsv.gz?dl=0) |
| MESO     | [MESO genomic data](https://www.dropbox.com/s/zwhq70c0hxhbeh4/meso_genomic_data.tsv.gz?dl=0) | [MESO SNVs](https://www.dropbox.com/s/2c4psbcgkyplpo5/meso_snv_data.tsv.gz?dl=0) | [MESO Indels](https://www.dropbox.com/s/msj59qe1r4vgqwp/meso_indel_data.tsv.gz?dl=0) | [MESO haploinsufficiency tables](https://www.dropbox.com/s/y1hx8ejesfs2jya/meso_haploins_cnnl_hemi.tsv.gz?dl=0) |
| OV       | [OV genomic data](https://www.dropbox.com/s/1pvn9wwi6bputjw/ov_genomic_data.tsv.gz?dl=0)     | [OV SNVs](https://www.dropbox.com/s/z8v7okde6fbx0sw/ov_snv_data.tsv.gz?dl=0)     | [OV Indels](https://www.dropbox.com/s/q5ekvcahcblfgk5/ov_indel_data.tsv.gz?dl=0)     | [OV haploinsufficiency tables](https://www.dropbox.com/s/1hugwr62a3dc8iy/ov_haploins_cnnl_hemi.tsv.gz?dl=0)     |
| PAAD     | [PAAD genomic data](https://www.dropbox.com/s/rs86jh5qw8rxtro/paad_genomic_data.tsv.gz?dl=0) | [PAAD SNVs](https://www.dropbox.com/s/9ar5z4hlayl0pf2/paad_snv_data.tsv.gz?dl=0) | [PAAD Indels](https://www.dropbox.com/s/xclda4e3ecsdxv1/paad_indel_data.tsv.gz?dl=0) | [PAAD haploinsufficiency tables](https://www.dropbox.com/s/uwy4av18egdccn5/paad_haploins_cnnl_hemi.tsv.gz?dl=0) |
| PCPG     | [PCPG genomic data](https://www.dropbox.com/s/elcnuiv8enc4r9p/pcpg_genomic_data.tsv.gz?dl=0) | [PCPG SNVs](https://www.dropbox.com/s/3jj5oa1c3eu04d8/pcpg_snv_data.tsv.gz?dl=0) | [PCPG Indels](https://www.dropbox.com/s/61hqzgwt2ui5ep4/pcpg_indel_data.tsv.gz?dl=0) | [PCPG haploinsufficiency tables](https://www.dropbox.com/s/rjzwqfy0zrzu9hd/pcpg_haploins_cnnl_hemi.tsv.gz?dl=0) |
| PRAD     | [PRAD genomic data](https://www.dropbox.com/s/0w0tuy9g7os4qka/prad_genomic_data.tsv.gz?dl=0) | [PRAD SNVs](https://www.dropbox.com/s/k655py6ed5wdglm/prad_snv_data.tsv.gz?dl=0) | [PRAD Indels](https://www.dropbox.com/s/t21g5vmkbhphlyh/prad_indel_data.tsv.gz?dl=0) | [PRAD haploinsufficiency tables](https://www.dropbox.com/s/6ojqfhrmfq8ilnr/prad_haploins_cnnl_hemi.tsv.gz?dl=0) |
| READ     | [READ genomic data](https://www.dropbox.com/s/rggv9h9enn8k14c/read_genomic_data.tsv.gz?dl=0) | [READ SNVs](https://www.dropbox.com/s/l1q3cb0beeh6f9e/read_snv_data.tsv.gz?dl=0) | [READ Indels](https://www.dropbox.com/s/7ri8dbifc38r212/read_indel_data.tsv.gz?dl=0) | [READ haploinsufficiency tables](https://www.dropbox.com/s/2h7eswp7h38gtzh/read_haploins_cnnl_hemi.tsv.gz?dl=0) |
| SARC     | [SARC genomic data](https://www.dropbox.com/s/eq9mpfezf7n98y3/sarc_genomic_data.tsv.gz?dl=0) | [SARC SNVs](https://www.dropbox.com/s/k1ylqdz36kvwmku/sarc_snv_data.tsv.gz?dl=0) | [SARC Indels](https://www.dropbox.com/s/ishf22cgcoqh6zg/sarc_indel_data.tsv.gz?dl=0) | [SARC haploinsufficiency tables](https://www.dropbox.com/s/mzo3m99nujkjn03/sarc_haploins_cnnl_hemi.tsv.gz?dl=0) |
| SKCM     | [SKCM genomic data](https://www.dropbox.com/s/pr523tfbkpr796d/skcm_genomic_data.tsv.gz?dl=0) | [SKCM SNVs](https://www.dropbox.com/s/vjp713o7hlpwk2r/skcm_snv_data.tsv.gz?dl=0) | [SKCM Indels](https://www.dropbox.com/s/zj5e6ll7vof9q5t/skcm_indel_data.tsv.gz?dl=0) | [SKCM haploinsufficiency tables](https://www.dropbox.com/s/nkk5z1tscc3pqnr/skcm_haploins_cnnl_hemi.tsv.gz?dl=0) |
| STAD     | [STAD genomic data](https://www.dropbox.com/s/k0vc0ovbwwdm03k/stad_genomic_data.tsv.gz?dl=0) | [STAD SNVs](https://www.dropbox.com/s/4yxrcmjr8tt3atd/stad_snv_data.tsv.gz?dl=0) | [STAD Indels](https://www.dropbox.com/s/6k4ssgroz5zdhmn/stad_indel_data.tsv.gz?dl=0) | [STAD haploinsufficiency tables](https://www.dropbox.com/s/auscosecc1g7to9/stad_haploins_cnnl_hemi.tsv.gz?dl=0) |
| TGCT     | [TGCT genomic data](https://www.dropbox.com/s/9g3r9h716mr29r2/tgct_genomic_data.tsv.gz?dl=0) | [TGCT SNVs](https://www.dropbox.com/s/1m6lbh7jpeh0y1z/tgct_snv_data.tsv.gz?dl=0) | [TGCT Indels](https://www.dropbox.com/s/bzoccs3m21489bq/tgct_indel_data.tsv.gz?dl=0) | [TGCT haploinsufficiency tables](https://www.dropbox.com/s/ptg8ip0gvcy7a9v/tgct_haploins_cnnl_hemi.tsv.gz?dl=0) |
| THCA     | [THCA genomic data](https://www.dropbox.com/s/hc74p2zy7xqedxw/thca_genomic_data.tsv.gz?dl=0) | [THCA SNVs](https://www.dropbox.com/s/mej8y1gbgvmycge/thca_snv_data.tsv.gz?dl=0) | [THCA Indels](https://www.dropbox.com/s/56osnho57cid4mp/thca_indel_data.tsv.gz?dl=0) | [THCA haploinsufficiency tables](https://www.dropbox.com/s/z91861lahjd7ur4/thca_haploins_cnnl_hemi.tsv.gz?dl=0) |
| UCEC     | [UCEC genomic data](https://www.dropbox.com/s/whecdo0vftlgv2f/ucec_genomic_data.tsv.gz?dl=0) | [UCEC SNVs](https://www.dropbox.com/s/2bgrfx8z3fprc2e/ucec_snv_data.tsv.gz?dl=0) | [UCEC Indels](https://www.dropbox.com/s/hq36tovv3dz51sa/ucec_indel_data.tsv.gz?dl=0) | [UCEC haploinsufficiency tables](https://www.dropbox.com/s/n6ptafs6faupxl2/ucec_haploins_cnnl_hemi.tsv.gz?dl=0) |
| UVM      | [UVM genomic data](https://www.dropbox.com/s/cypbvc4qopzgc0g/uvm_genomic_data.tsv.gz?dl=0)   | [UVM SNVs](https://www.dropbox.com/s/a9yhha2um805hru/uvm_snv_data.tsv.gz?dl=0)   | [UVM Indels](https://www.dropbox.com/s/k0o46nmz6oihft9/uvm_indel_data.tsv.gz?dl=0)   | [UVM haploinsufficiency tables](https://www.dropbox.com/s/hgmjexry5bfyoxb/uvm_haploins_cnnl_hemi.tsv.gz?dl=0)   |

### Table descriptions

The following sections describe the content of each file linked in the table
above.

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

#### Haploinsufficiency data

The Haploinsufficiency data table contains informations about the
haploinsufficiency of gene expression tested in different genomic contexts
(i.e. hemizygous deletion and copy neutral loss of heterozygosity).

In the table below the field present in the table are described.

| Field               | Description                                                                                                                          |
| --------            | --------                                                                                                                             |
| gene_id             | The symbol of the gene                                                                                                               |
| homo_del            | The number of samples with homozygous deletions of the gene                                                                          |
| hemi_del            | The number of samples with hemizygous deletions of the gene                                                                          |
| cnnl                | The number of samples with copy neutral losses of heterozygosity of the gene                                                         |
| gain_del            | The number of samples with copy gains loss of heterozygosity of the gene                                                             |
| wt                  | The number of samples with normal copy number                                                                                        |
| gain_unb            | The number of samples with an unbalanced gain of the gene                                                                            |
| gain                | The number of samples with a balanced gain of the gene                                                                               |
| amp                 | The number of samples with a balanced amplification of the gene                                                                      |
| amp_del             | The number of samples with copy amplifications loss of heterozygosity of the gene                                                    |
| amp_unb             | The number of samples with a unbalanced amplification of the gene                                                                    |
| nd                  | The number of samples with a undetermied allele specific state of the gene                                                           |
| n_tot               | The total number of considered samples                                                                                               |
| n_snv               | The total number of SNVs in the samples                                                                                              |
| n_double_hit        | The number of events where the SNV is present in LoH genomic context (hemi_del, cnnl, gain_del)                                      |
| mean_fpkm           | The mean expression of the gene in FPKM                                                                                              |
| mean_fpkm_wt        | The mean expression of the gene in FPKM in wt genomic context                                                                        |
| mean_fpkm_hemi_del  | The mean expression of the gene in FPKM in hemi_del genomic context                                                                  |
| mean_fpkm_cnnl      | The mean expression of the gene in FPKM in cnnl genomic context                                                                      |
| mean_fpkm_hetloss   | The mean expression of the gene in FPKM in LoH genomic context (either hemi_del or cnnl)                                             |
| pval_hemi_del_vs_wt | The p-value for the test of haploinsufficiency between samples with hemi_del vs samples with wt genomic context                      |
| pval_cnnl_vs_wt     | The p-value for the test of haploinsufficiency between samples with cnnl vs samples with wt genomic context                          |
| pval_hetloss_vs_wt  | The p-value for the test of haploinsufficiency between samples with LoH (either hemi_del or cnnl) vs samples with wt genomic context |

## Code

The code folder contains the code for the pipeline.

### Pipeline

The pipeline folder contains the complete code of analysis pipeline described
in the paper.

#### Running the pipeline

In order to run the pipeline some resources files and tools are needed. A
complete example folder that includes both the code and the resources can be
downloaded from [here](https://www.dropbox.com/s/ur4wq6zp0xtwbus/pipeline_example.tar.gz?dl=0).
Be aware that the file is big (12GB). The complete instructions needed to
run the pipeline are inside the archive.

## Fundings

This project is funded by (ERC-CoG-2014-648670)
