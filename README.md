![An Image](/images/SPICE_COLORI.png)

# SPICE-pipeline
Framework for allele specific analysis of matched tumor and normal next generation sequencing data

## Pre-computed data

All precomputed TCGA allele-specific data and FAME tables are available [here](https://www.dropbox.com/sh/jzk9roi8c238cc0/AABTTlO7YTneLrjyVUIV8vJia?dl=0)

The table below contains the links to all the output files for each dataset:

| Dataset  | Genomic data                                                                                 | FAME: Hemi-Del vs Hemi-Del                                                                                      | FAME: CN-LOH vs CN-LOH                                                                                  | FAME: (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)                                                                                       |
| -------- | --------                                                                                     | --------                                                                                                        | --------                                                                                                | --------                                                                                                                                 |
| ACC      | [ACC genomic data](https://www.dropbox.com/s/71dcmiv1k7vqo2m/acc_genomic_data.tsv.gz?dl=0)   | [FAME ACC Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/q8sw97cmwlf5eoy/acc_hemi_del_fame_data.tsv.gz?dl=0)   | [FAME ACC CN-LOH vs CN-LOH](https://www.dropbox.com/s/fx0l35esh13esv4/acc_cnnl_fame_data.tsv.gz?dl=0)   | [FAME ACC (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/9mghkoelyoqaijr/acc_hemi_cnnl_fame_data.tsv.gz?dl=0)   |
| BLCA     | [BLCA genomic data](https://www.dropbox.com/s/854v67lzvvedbqp/blca_genomic_data.tsv.gz?dl=0) | [FAME BLCA Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/mvb5siimw86fvqi/blca_hemi_del_fame_data.tsv.gz?dl=0) | [FAME BLCA CN-LOH vs CN-LOH](https://www.dropbox.com/s/tfhtp9f6s39b8dr/blca_cnnl_fame_data.tsv.gz?dl=0) | [FAME BLCA (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/dep6ck00zfwrokc/blca_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| BRCA     | [BRCA genomic data](https://www.dropbox.com/s/qu7xrmmjde7606m/brca_genomic_data.tsv.gz?dl=0) | [FAME BRCA Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/69b2htnxfjisxy2/brca_hemi_del_fame_data.tsv.gz?dl=0) | [FAME BRCA CN-LOH vs CN-LOH](https://www.dropbox.com/s/pebgf2d5b8auysc/brca_cnnl_fame_data.tsv.gz?dl=0) | [FAME BRCA (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/slrbo1jdhfnhg4e/brca_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| COAD     | [COAD genomic data](https://www.dropbox.com/s/vtqf1lndvixispa/coad_genomic_data.tsv.gz?dl=0) | [FAME COAD Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/uvfo1lf3xu3p62x/coad_hemi_del_fame_data.tsv.gz?dl=0) | [FAME COAD CN-LOH vs CN-LOH](https://www.dropbox.com/s/x457f1ypxd9xe7q/coad_cnnl_fame_data.tsv.gz?dl=0) | [FAME COAD (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/5pdmhsd5uqhqt2u/coad_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| GBM      | [GBM genomic data](https://www.dropbox.com/s/0y6sr1dn12ntvuc/gbm_genomic_data.tsv.gz?dl=0)   | [FAME GBM Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/nlefjl67mkzvgbv/gbm_hemi_del_fame_data.tsv.gz?dl=0)   | [FAME GBM CN-LOH vs CN-LOH](https://www.dropbox.com/s/qgdgp49jpr6yubm/gbm_cnnl_fame_data.tsv.gz?dl=0)   | [FAME GBM (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/qympl5c8nsalez8/gbm_hemi_cnnl_fame_data.tsv.gz?dl=0)   |
| HNSC     | [HNSC genomic data](https://www.dropbox.com/s/vuwuiwjg0k70qoy/hnsc_genomic_data.tsv.gz?dl=0) | [FAME HNSC Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/z25v9kxnb35tbvo/hnsc_hemi_del_fame_data.tsv.gz?dl=0) | [FAME HNSC CN-LOH vs CN-LOH](https://www.dropbox.com/s/89fsnkl8zu9wn84/hnsc_cnnl_fame_data.tsv.gz?dl=0) | [FAME HNSC (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/f761cpu6dbwf8m0/hnsc_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| KICH     | [KICH genomic data](https://www.dropbox.com/s/9afh9e3ecp94drz/kich_genomic_data.tsv.gz?dl=0) | [FAME KICH Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/4goxerilm3fsh7y/kich_hemi_del_fame_data.tsv.gz?dl=0) | [FAME KICH CN-LOH vs CN-LOH](https://www.dropbox.com/s/hjfy891aui0j296/kich_cnnl_fame_data.tsv.gz?dl=0) | [FAME KICH (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/811n645km6yp34f/kich_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| KIRC     | [KIRC genomic data](https://www.dropbox.com/s/v7h2h4of3znyin5/kirc_genomic_data.tsv.gz?dl=0) | [FAME KIRC Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/qrh7ys9unkoipfn/kirc_hemi_del_fame_data.tsv.gz?dl=0) | [FAME KIRC CN-LOH vs CN-LOH](https://www.dropbox.com/s/48t655g3ho6cvde/kirc_cnnl_fame_data.tsv.gz?dl=0) | [FAME KIRC (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/f3vqp5vx5ydfjf7/kirc_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| KIRP     | [KIRP genomic data](https://www.dropbox.com/s/lgo88v8d113j9l4/kirp_genomic_data.tsv.gz?dl=0) | [FAME KIRP Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/8a2xzobbh5sh0ie/kirp_hemi_del_fame_data.tsv.gz?dl=0) | [FAME KIRP CN-LOH vs CN-LOH](https://www.dropbox.com/s/28qur1ws0o4t89l/kirp_cnnl_fame_data.tsv.gz?dl=0) | [FAME KIRP (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/c06j3r3jhrkoqlf/kirp_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| LAML     | [LAML genomic data](https://www.dropbox.com/s/vwb6gmmnfne6tyn/laml_genomic_data.tsv.gz?dl=0) | No data available                                                                                               | No data available                                                                                       | No data available                                                                                                                        |
| LGG      | [LGG genomic data](https://www.dropbox.com/s/hcchp4x6ngf3p2y/lgg_genomic_data.tsv.gz?dl=0)   | [FAME LGG Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/44h0632ctr3ihrc/lgg_hemi_del_fame_data.tsv.gz?dl=0)   | [FAME LGG CN-LOH vs CN-LOH](https://www.dropbox.com/s/4qlpz3ybow1fbxt/lgg_cnnl_fame_data.tsv.gz?dl=0)   | [FAME LGG (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/f0dvysqabvb0jea/lgg_hemi_cnnl_fame_data.tsv.gz?dl=0)   |
| LIHC     | [LIHC genomic data](https://www.dropbox.com/s/8708j0s72tvpuqa/lihc_genomic_data.tsv.gz?dl=0) | [FAME LIHC Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/1yttm59n327if40/lihc_hemi_del_fame_data.tsv.gz?dl=0) | [FAME LIHC CN-LOH vs CN-LOH](https://www.dropbox.com/s/3fh6of3pcaxe5t9/lihc_cnnl_fame_data.tsv.gz?dl=0) | [FAME LIHC (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/r3vofxrldt96l0a/lihc_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| LUAD     | [LUAD genomic data](https://www.dropbox.com/s/4l54uwwcxfh31pm/luad_genomic_data.tsv.gz?dl=0) | [FAME LUAD Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/eogq0cdt3g0el5d/luad_hemi_del_fame_data.tsv.gz?dl=0) | [FAME LUAD CN-LOH vs CN-LOH](https://www.dropbox.com/s/4jugmt9sjsygq96/luad_cnnl_fame_data.tsv.gz?dl=0) | [FAME LUAD (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/2rggkfzmda25t0m/luad_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| LUSC     | [LUSC genomic data](https://www.dropbox.com/s/y8c9mlkaw7ch9f3/lusc_genomic_data.tsv.gz?dl=0) | [FAME LUSC Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/4tdu280ghyzhuwq/lusc_hemi_del_fame_data.tsv.gz?dl=0) | [FAME LUSC CN-LOH vs CN-LOH](https://www.dropbox.com/s/8tlfdwkh6e27n3f/lusc_cnnl_fame_data.tsv.gz?dl=0) | [FAME LUSC (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/7hcxr5a3lypppc0/lusc_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| MESO     | [MESO genomic data](https://www.dropbox.com/s/9znywdispg53igs/meso_genomic_data.tsv.gz?dl=0) | [FAME MESO Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/o1bnm0ys5kukjaz/meso_hemi_del_fame_data.tsv.gz?dl=0) | [FAME MESO CN-LOH vs CN-LOH](https://www.dropbox.com/s/wzg3owevnyy7gfs/meso_cnnl_fame_data.tsv.gz?dl=0) | [FAME MESO (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/z45srwbgoxmgqm2/meso_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| OV       | [OV genomic data](https://www.dropbox.com/s/3vuzued5td7q98w/ov_genomic_data.tsv.gz?dl=0)     | [FAME OV Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/up4d8s0a4xs058z/ov_hemi_del_fame_data.tsv.gz?dl=0)     | [FAME OV CN-LOH vs CN-LOH](https://www.dropbox.com/s/mrnawxiw16th2tb/ov_cnnl_fame_data.tsv.gz?dl=0)     | [FAME OV (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/coo8ffi5vwypci7/ov_hemi_cnnl_fame_data.tsv.gz?dl=0)     |
| PAAD     | [PAAD genomic data](https://www.dropbox.com/s/zhl1x36lamlezr6/paad_genomic_data.tsv.gz?dl=0) | [FAME PAAD Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/ctu5x8owf9rmiap/paad_hemi_del_fame_data.tsv.gz?dl=0) | [FAME PAAD CN-LOH vs CN-LOH](https://www.dropbox.com/s/d3541810g93dqmi/paad_cnnl_fame_data.tsv.gz?dl=0) | [FAME PAAD (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/gwzcr4f6mv3nol1/paad_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| PCPG     | [PCPG genomic data](https://www.dropbox.com/s/o4k1z46kmm5yd38/pcpg_genomic_data.tsv.gz?dl=0) | [FAME PCPG Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/illmhoeeix1fx6a/pcpg_hemi_del_fame_data.tsv.gz?dl=0) | [FAME PCPG CN-LOH vs CN-LOH](https://www.dropbox.com/s/2facgop1dcxwasg/pcpg_cnnl_fame_data.tsv.gz?dl=0) | [FAME PCPG (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/nfd4qmkl1urcl4b/pcpg_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| PRAD     | [PRAD genomic data](https://www.dropbox.com/s/2d72kuco8m4qaq8/prad_genomic_data.tsv.gz?dl=0) | [FAME PRAD Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/iu1ka7jlu75sgi8/prad_hemi_del_fame_data.tsv.gz?dl=0) | [FAME PRAD CN-LOH vs CN-LOH](https://www.dropbox.com/s/hketar8bz44ed13/prad_cnnl_fame_data.tsv.gz?dl=0) | [FAME PRAD (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/5qs0grsgojzu4ky/prad_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| READ     | [READ genomic data](https://www.dropbox.com/s/35j4nns9qvrm2l7/read_genomic_data.tsv.gz?dl=0) | [FAME READ Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/6tez8t7l31016yr/read_hemi_del_fame_data.tsv.gz?dl=0) | [FAME READ CN-LOH vs CN-LOH](https://www.dropbox.com/s/5eb0qv8jckqddqu/read_cnnl_fame_data.tsv.gz?dl=0) | [FAME READ (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/5xelio7f19grm22/read_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| SARC     | [SARC genomic data](https://www.dropbox.com/s/c4rm5z4jp1wki99/sarc_genomic_data.tsv.gz?dl=0) | [FAME SARC Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/6j773fp81zj0b8r/sarc_hemi_del_fame_data.tsv.gz?dl=0) | [FAME SARC CN-LOH vs CN-LOH](https://www.dropbox.com/s/2dp5osp4mlnso6r/sarc_cnnl_fame_data.tsv.gz?dl=0) | [FAME SARC (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/85q620tdk7qv9u9/sarc_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| SKCM     | [SKCM genomic data](https://www.dropbox.com/s/frkqol7eg7fdu8o/skcm_genomic_data.tsv.gz?dl=0) | [FAME SKCM Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/idty0735cmnyxqq/skcm_hemi_del_fame_data.tsv.gz?dl=0) | [FAME SKCM CN-LOH vs CN-LOH](https://www.dropbox.com/s/s2ait7sr0bp3gkq/skcm_cnnl_fame_data.tsv.gz?dl=0) | [FAME SKCM (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/gd7yai7i2ss7aj0/skcm_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| STAD     | [STAD genomic data](https://www.dropbox.com/s/fz6y3it6xi6sz3g/stad_genomic_data.tsv.gz?dl=0) | [FAME STAD Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/2938sp0mx9t4fxc/stad_hemi_del_fame_data.tsv.gz?dl=0) | [FAME STAD CN-LOH vs CN-LOH](https://www.dropbox.com/s/vqrk2mwrulu4uls/stad_cnnl_fame_data.tsv.gz?dl=0) | [FAME STAD (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/s6uyfbuowofq37q/stad_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| TGCT     | [TGCT genomic data](https://www.dropbox.com/s/77q19y487xg1d8h/tgct_genomic_data.tsv.gz?dl=0) | [FAME TGCT Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/3tan4zd6z8edtdo/tgct_hemi_del_fame_data.tsv.gz?dl=0) | [FAME TGCT CN-LOH vs CN-LOH](https://www.dropbox.com/s/ldsieevr9n8ml3y/tgct_cnnl_fame_data.tsv.gz?dl=0) | [FAME TGCT (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/jtkvkl3e9gy93iv/tgct_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| THCA     | [THCA genomic data](https://www.dropbox.com/s/rp4aq5dkgywc4h1/thca_genomic_data.tsv.gz?dl=0) | [FAME THCA Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/wnsc8a8lpq202qt/thca_hemi_del_fame_data.tsv.gz?dl=0) | [FAME THCA CN-LOH vs CN-LOH](https://www.dropbox.com/s/xjv8z66ixlmr7rf/thca_cnnl_fame_data.tsv.gz?dl=0) | [FAME THCA (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/epm5caossly67a8/thca_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| UCEC     | [UCEC genomic data](https://www.dropbox.com/s/cg4o0y3mu5phs3b/ucec_genomic_data.tsv.gz?dl=0) | [FAME UCEC Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/zsm1rhukctap0ns/ucec_hemi_del_fame_data.tsv.gz?dl=0) | [FAME UCEC CN-LOH vs CN-LOH](https://www.dropbox.com/s/1n0myixwgm9f561/ucec_cnnl_fame_data.tsv.gz?dl=0) | [FAME UCEC (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/89cffw29h3rns6u/ucec_hemi_cnnl_fame_data.tsv.gz?dl=0) |
| UVM      | [UVM genomic data](https://www.dropbox.com/s/qjbg1rfquyp702d/uvm_genomic_data.tsv.gz?dl=0)   | [FAME UVM Hemi-Del vs Hemi-Del](https://www.dropbox.com/s/tjf7ta53blfy7cr/uvm_hemi_del_fame_data.tsv.gz?dl=0)   | [FAME UVM CN-LOH vs CN-LOH](https://www.dropbox.com/s/crvh9hwnrd1fbud/uvm_cnnl_fame_data.tsv.gz?dl=0)   | [FAME UVM (Hemi-Del or CN-LOH) vs (Hemi-Del or CN-LOH)](https://www.dropbox.com/s/vykat6t1160gfcx/uvm_hemi_cnnl_fame_data.tsv.gz?dl=0)   |

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

#### FAME data

Each FAME output table contains the results of the independence analyses (FAME)
run on each pair of cancer genes vs druggable (as described in the methods
section of the paper). The analyses were run on three different sets of
aberrations:

- Hemizygous deletions (Hemi-Del): one allele is lost and the other retains one
  copy
- Copy neutral Loss of heterozygosity (CN-LOH): one allele is lost while the
  other has a duplication (in non-allele specific analyses this would appear as
  wild type)
- Either Hemi-Del or CN-LOH: the sample is in either of these states.

Each FAME results table contains all the pairs of genes (with at least one
aberration each) filtered to have odds-ratio < 1.

In the table below each field of the results is described:

| Field      | Description                                                                               |
| --------   | --------                                                                                  |
| a1         | The type of aberrations considered for gene 1                                             |
| a2         | The type of aberrations considered for gene 2                                             |
| dataset    | The name of the dataset                                                                   |
| g1         | The symbol for the gene 1                                                                 |
| g2         | The symbol for the gene 2                                                                 |
| p_value    | The p-value for the fisher test                                                           |
| fdr        | The FDR computed on the entire result (method: Benjamini-Hochberg)                        |
| odds_ratio | The odds-ratio statistics for the test                                                    |
| f_min      | The minimum frequency between the aberration of the two genes (f_min = min(f_1, f_2))     |
| n_11       | The number of samples for which gene 1 have aberration a_1 and gene 2 have aberration a_2 |
| n_10       | The number of samples for which gene 1 has aberration a_1 while gene 2 is not             |
| n_01       | The number of samples for which gene 2 has aberration a_2 while gene 1 is not             |
| n_00       | The number of samples for which no gene is aberrant                                       |
| n_tot      | The total number of samples (n_tot = n_11 + n_10 + n_01 + n_00)                           |
| n_1        | The number of samples with aberration a_1 in gene 1 (n_1 = n_11 + n_10)                   |
| n_2        | The number of samples with aberration a_2 in gene 2 (n_2 = n_11 + n_01)                   |
| f_1        | The frequency of samples with aberration a_1 in gene 1 (f_1 = n_1 / n_tot)                |
| f_2        | The frequency of samples with aberration a_2 in gene 2 (f_2 = n_1 / n_tot)                |

## Code

The code folder contains the code for the pipeline and FAME

### Pipeline

The pipeline folder contains the complete code of analysis pipeline described
in the paper.

#### Running the pipeline

In order to run the pipeline some resources files and tools are needed. A
complete example folder that includes both the code and the resources can be
downloaded from [here](https://www.dropbox.com/s/ur4wq6zp0xtwbus/pipeline_example.tar.gz?dl=0).
Be aware that the file is big (12GB). The complete instructions needed to
run the pipeline are inside the archive.

### FaME

FAME is a tool that can test the aberrations of the samples in one or more
datasets in order to find signals of mutual exclusivity or co-occurrence
between all the pairs of genes. FaME can test many different combinations of
aberrations very quickly by using several techniques described in the methods
section of the paper.

#### Running FaME

> :warning: FaME requires large amount of computational resources if run
> on large datasets.

In order to run FaME just copy the genomic data that can be downloaded from the
links reported above to the data/genomic/ folder and run the code. The script
requires the R *librarian* package. The *librarian* package can be installed with
the following command:

```r
install.packages('librarian')
```

In order to run the analysis just run the script `fame_analysis.R`. When
launched the script will automatically install the packages that are needed and
then the analysis will start.

## Fundings

This project is funded by (ERC-CoG-2014-648670)
